import os
import random
import requests
from Bio import Entrez
import time
import google.generativeai as genai
from google.api_core import exceptions as google_exceptions

keyword_groups = {
    "A": [
        "ancient DNA", "semi-nested PCR", "mitochondrial DNA", "DNA extraction",
        "next-generation sequencing", "phylogenetics", "species identification",
        "genetic diversity", "haplotype", "fossil DNA",
        "molecular evolution", "genetic distance", "paleogenomics",
        "bone sample analysis", "DNA degradation"
    ],
    "B": [
        "developmental biology", "gastrulation", "self-organization", "organoids",
        "deep learning", "morphogenesis", "synthetic embryo", "morphogen gradients",
        "3D cell culture", "AI in biology", "computational modeling", "biophysical modeling",
        "embryoid body", "pattern formation", "spatial gene expression"
    ],
    "C": [
        "casein plastic", "biopolymer", "protein-based material", "leather alternative",
        "sustainable materials", "mechanical testing", "tensile strength", "biodegradable plastic",
        "natural polymer", "galalith", "organic materials", "hardness test", "material properties"
    ]
}


# 環境設定
Entrez.email = os.environ.get("EMAIL_ADDRESS")
if not Entrez.email:
    raise ValueError("EMAIL_ADDRESS is not set in environment variables")

SLACK_WEBHOOK_URL = os.environ.get("SLACK_WEBHOOK_URL")
if not SLACK_WEBHOOK_URL:
    raise ValueError("SLACK_WEBHOOK_URL is not set in environment variables")

GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if not GEMINI_API_KEY:
    raise ValueError("GEMINI_API_KEY is not set in environment variables")
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel("gemini-1.5-flash")

def translate_to_japanese(text, max_retries=3):
    for attempt in range(max_retries):
        try:
            prompt = f"Translate the following text to Japanese:\n\n{text}"
            response = model.generate_content(prompt)
            translated_text = response.text.strip()
            return translated_text
        except google_exceptions.ResourceExhausted as e:
            print(f"Resource exhausted on attempt {attempt + 1}: {str(e)}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print("Max retries reached for resource exhausted, using original text.")
                return text
        except google_exceptions.TooManyRequests as e:
            print(f"Rate limit exceeded on attempt {attempt + 1}: {str(e)}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print("Max retries reached for rate limit, using original text.")
                return text
        except Exception as e:
            print(f"Translation failed on attempt {attempt + 1}: {str(e)}")
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                print("Max retries reached, using original text.")
                return text
    return text

def load_processed_pmids(filename="processed_pmids.txt"):
    if os.path.exists(filename):
        with open(filename, "r") as f:
            return set(f.read().splitlines())
    return set()

def save_processed_pmids(pmids, filename="processed_pmids.txt"):
    with open(filename, "a") as f:
        for pmid in pmids:
            f.write(f"{pmid}\n")

def search_pubmed(term, max_results=1, processed_pmids=None, max_retries=3):
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            pmids = [pmid for pmid in record["IdList"] if processed_pmids is None or pmid not in processed_pmids]
            return pmids[:max_results]
        except Exception as e:
            print(f"PubMed search failed for '{term}' on attempt {attempt + 1}: {str(e)}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print(f"Max retries reached for PubMed search, skipping term: {term}")
                # エラーをSlackに通知（失敗してもプロセスを中断しない）
                try:
                    post_to_slack(SLACK_WEBHOOK_URL, f"PubMed search failed for '{term}' after {max_retries} attempts: {str(e)}")
                except Exception as slack_error:
                    print(f"Failed to send PubMed search error to Slack: {str(slack_error)}")
                return []

def fetch_details(pmids):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        papers = []
        for article in records["PubmedArticle"]:
            pmid = article["MedlineCitation"]["PMID"]
            title = article["MedlineCitation"]["Article"]["ArticleTitle"]
            abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No abstract"])[0]
            papers.append({"pmid": pmid, "title": title, "abstract": abstract})
        return papers
    except Exception as e:
        print(f"PubMed fetch failed for PMIDs {pmids}: {str(e)}")
        try:
            post_to_slack(SLACK_WEBHOOK_URL, f"PubMed fetch failed for PMIDs {pmids}: {str(e)}")
        except Exception as slack_error:
            print(f"Failed to send PubMed fetch error to Slack: {str(slack_error)}")
        return []

def post_to_slack(webhook_url, message, thread_ts=None):
    try:
        payload = {"text": message}
        if thread_ts:
            payload["thread_ts"] = thread_ts
        response = requests.post(webhook_url, json=payload, timeout=10)
        print(f"Response status: {response.status_code}, Response text: {response.text}")
        response.raise_for_status()
        try:
            return response.json().get("ts", None)
        except json.JSONDecodeError:
            print("JSON decode failed, assuming successful notification.")
            return None
    except requests.exceptions.HTTPError as e:
        print(f"HTTP Error: {str(e)}, Response: {response.text}")
        raise  # 再スローして上位で処理
    except Exception as e:
        print(f"Slack notification failed: {str(e)}")
        raise  # 再スローして上位で処理

def main():
    try:
        # 処理済みPMIDを読み込み
        processed_pmids = load_processed_pmids()
        new_pmids = set()
        skip_translation = False

        # 各グループから3つのキーワードを選択
        selected_keywords = {}
        for group, keywords in keyword_groups.items():
            selected_keywords[group] = random.sample(keywords, 3)

        papers_with_metadata = []

        # 各グループごとに処理（1論文ずつ）
        for group, keywords in selected_keywords.items():
            combined_term = " AND ".join(keywords)
            pmids = search_pubmed(combined_term, max_results=1, processed_pmids=processed_pmids)
            if not pmids:
                print(f"No papers found for '{combined_term}' in Group {group}, skipping...")
                try:
                    post_to_slack(SLACK_WEBHOOK_URL, f"No papers found for '{combined_term}' in Group {group}")
                except Exception as e:
                    print(f"Failed to send Slack message for no papers: {str(e)}")
                continue
            papers = fetch_details(pmids)
            group_metadata = []
            for paper in papers:
                group_metadata.append({
                    "paper": paper,
                    "group": group,
                    "keywords": keywords
                })
                new_pmids.add(paper["pmid"])
            selected_papers = group_metadata[:1]
            papers_with_metadata.extend(selected_papers)
            time.sleep(2)

        # 新しいPMIDを保存
        if new_pmids:
            save_processed_pmids(new_pmids)

        # Slackメッセージ作成
        if not papers_with_metadata:
            slack_message = "No papers found this week."
            try:
                post_to_slack(SLACK_WEBHOOK_URL, slack_message)
            except Exception as e:
                print(f"Failed to send Slack message for no papers: {str(e)}")
        else:
            thread_ts = None
            try:
                thread_ts = post_to_slack(SLACK_WEBHOOK_URL, "This Week's Papers (3 articles):")
            except Exception as e:
                print(f"Failed to send initial Slack message: {str(e)}")
                thread_ts = None

            for i, meta in enumerate(papers_with_metadata[:3], 1):
                paper = meta["paper"]
                group = meta["group"]
                keywords = ", ".join(meta["keywords"])
                title = paper["title"]
                abstract = paper["abstract"]
                link = f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/"

                # 翻訳をスキップするか判定
                title_ja = translate_to_japanese(title) if not skip_translation else title
                abstract_ja = translate_to_japanese(abstract) if not skip_translation else abstract

                message = (
                    f"**{i}. {title_ja}**\n"
                    f"Group: {group}, Keywords: {keywords}\n"
                    f"要約: {abstract_ja}\n"
                    f"Link: {link}"
                )
                try:
                    post_to_slack(SLACK_WEBHOOK_URL, message, thread_ts)
                except Exception as e:
                    print(f"Failed to send Slack message for paper {i}: {str(e)}")

        # PMIDファイルをリポジトリにコミット
        if new_pmids:
            os.system('git config --global user.email "github-actions@users.noreply.github.com"')
            os.system('git config --global user.name "GitHub Actions"')
            os.system("git add processed_pmids.txt")
            os.system('git commit -m "Update processed PMIDs"')
            os.system("git push origin main || echo 'Git push failed, continuing...'")

    except google_exceptions.ResourceExhausted as e:
        print(f"Resource exhausted globally: {str(e)}")
        try:
            post_to_slack(SLACK_WEBHOOK_URL, "Gemini API quota exceeded, skipping translation.")
        except Exception as e2:
            print(f"Failed to send Slack message for quota error: {str(e2)}")
        skip_translation = True
        # 既存の論文データで処理を継続
        if not thread_ts:
            try:
                thread_ts = post_to_slack(SLACK_WEBHOOK_URL, "This Week's Papers (3 articles):")
            except Exception as e:
                print(f"Failed to send initial Slack message after quota error: {str(e)}")
                thread_ts = None

        for i, meta in enumerate(papers_with_metadata[:3], 1):
            paper = meta["paper"]
            group = meta["group"]
            keywords = ", ".join(meta["keywords"])
            title = paper["title"]
            abstract = paper["abstract"]
            link = f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/"

            title_ja = title
            abstract_ja = abstract

            message = (
                f"**{i}. {title_ja}**\n"
                f"Group: {group}, Keywords: {keywords}\n"
                f"要約: {abstract_ja}\n"
                f"Link: {link}"
            )
            try:
                post_to_slack(SLACK_WEBHOOK_URL, message, thread_ts)
            except Exception as e:
                print(f"Failed to send Slack message for paper {i}: {str(e)}")
    except Exception as e:
        error_msg = f"System error: {str(e)}"
        try:
            post_to_slack(SLACK_WEBHOOK_URL, error_msg)
        except Exception as e2:
            print(f"Failed to send error to Slack: {str(e2)}")
        raise

if __name__ == "__main__":
    main()
