import os
import random
import requests
from Bio import Entrez
import time
import google.generativeai as genai
from google.api_core import exceptions as google_exceptions
import json  # JSONデコード用に追加

# キーワードグループ
keyword_groups = {
    "A": ["semi-nested PCR", "mitochondrial DNA", "amplification", "sequencing", "phylogenetics", "species identification", "genetic diversity", "haplotype"],
    "B": ["bioinformatics", "sequence alignment", "phylogenetic tree", "population genetics", "comparative genomics", "haplotype network", "genome annotation"],
    "C": ["Yamanaka factors", "rejuvenation", "epigenetic reprogramming", "iPS", "cellular senescence", "aging reversal", "transcription factor reprogramming", "epigenetic clocks", "regenerative medicine"]
}

# 環境設定
Entrez.email = "10311kaduken@gmail.com"  # 実際のメールアドレスに変更
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
        except google_exceptions.QuotaExceeded as e:
            print(f"Quota exceeded on attempt {attempt + 1}: {str(e)}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print("Max retries reached for quota exceeded, using original text.")
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

def summarize_abstract(text, max_retries=3):
    for attempt in range(max_retries):
        try:
            prompt = (
                f"以下の日本語の要約を簡潔に（100文字程度）要約してください。重要なポイントを保持してください。\n\n{text}"
            )
            response = model.generate_content(prompt)
            summarized_text = response.text.strip()
            return summarized_text
        except google_exceptions.QuotaExceeded as e:
            print(f"Quota exceeded on attempt {attempt + 1} during summarization: {str(e)}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print("Max retries reached for quota exceeded during summarization, using original text.")
                return text[:100] + "..." if len(text) > 100 else text
        except google_exceptions.TooManyRequests as e:
            print(f"Rate limit exceeded on attempt {attempt + 1} during summarization: {str(e)}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print("Max retries reached for rate limit during summarization, using original text.")
                return text[:100] + "..." if len(text) > 100 else text
        except Exception as e:
            print(f"Summarization failed on attempt {attempt + 1}: {str(e)}")
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                print("Max retries reached for summarization, using original text.")
                return text[:100] + "..." if len(text) > 100 else text
    return text[:100] + "..." if len(text) > 100 else text

def load_processed_pmids(filename="processed_pmids.txt"):
    if os.path.exists(filename):
        with open(filename, "r") as f:
            return set(f.read().splitlines())
    return set()

def save_processed_pmids(pmids, filename="processed_pmids.txt"):
    with open(filename, "a") as f:
        for pmid in pmids:
            f.write(f"{pmid}\n")

def search_pubmed(term, max_results=2, processed_pmids=None):
    try:
        handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        pmids = [pmid for pmid in record["IdList"] if processed_pmids is None or pmid not in processed_pmids]
        return pmids[:max_results]
    except Exception as e:
        post_to_slack(SLACK_WEBHOOK_URL, f"PubMed search failed for '{term}': {str(e)}")
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
            if len(abstract) > 500:
                abstract = abstract[:500] + "... (see link for full text)"
            papers.append({"pmid": pmid, "title": title, "abstract": abstract})
        return papers
    except Exception as e:
        post_to_slack(SLACK_WEBHOOK_URL, f"PubMed fetch failed for PMIDs {pmids}: {str(e)}")
        return []

def post_to_slack(webhook_url, message, thread_ts=None):
    try:
        payload = {"text": message}
        if thread_ts:
            payload["thread_ts"] = thread_ts
        response = requests.post(webhook_url, json=payload, timeout=10)
        print(f"Response status: {response.status_code}, Response text: {response.text}")
        response.raise_for_status()
        # JSONデコードを試行し、失敗時はNoneを返す
        try:
            return response.json().get("ts", None)
        except json.JSONDecodeError:
            print("JSON decode failed, assuming successful notification.")
            return None  # tsを取得しない場合も処理を継続
    except requests.exceptions.HTTPError as e:
        print(f"HTTP Error: {str(e)}, Response: {response.text}")
        raise
    except Exception as e:
        print(f"Slack notification failed: {str(e)}")
        raise

def main():
    try:
        # 処理済みPMIDを読み込み
        processed_pmids = load_processed_pmids()
        new_pmids = set()
        skip_summary = False  # 要約スキップフラグ
        thread_ts = None  # 初期化

        # 各グループから3つのキーワードを選択
        selected_keywords = {}
        for group, keywords in keyword_groups.items():
            selected_keywords[group] = random.sample(keywords, 3)

        papers_with_metadata = []

        # 各グループごとに処理
        for group, keywords in selected_keywords.items():
            combined_term = " AND ".join(keywords)
            pmids = search_pubmed(combined_term, max_results=6, processed_pmids=processed_pmids)
            if not pmids:
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
            selected_papers = group_metadata[:2]
            papers_with_metadata.extend(selected_papers)
            time.sleep(2)  # PubMed APIレート制限対策

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
            try:
                thread_ts = post_to_slack(SLACK_WEBHOOK_URL, "This Week's Papers (6 articles):")
            except Exception as e:
                print(f"Failed to send initial Slack message: {str(e)}")
                thread_ts = None

            for i, meta in enumerate(papers_with_metadata[:6], 1):
                paper = meta["paper"]
                group = meta["group"]
                keywords = ", ".join(meta["keywords"])
                title = paper["title"]
                abstract = paper["abstract"]
                link = f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/"

                title_ja = translate_to_japanese(title)
                abstract_ja = translate_to_japanese(abstract)
                abstract_summary_ja = summarize_abstract(abstract_ja) if not skip_summary else abstract_ja[:100] + "..."

                message = (
                    f"**{i}. {title_ja} ({title})**\n"
                    f"Group: {group}, Keywords: {keywords}\n"
                    f"要約（簡潔）: {abstract_summary_ja}\n"
                    f"要約（詳細）: {abstract_ja}\n"
                    f"原文: {abstract}\n"
                    f"Link: {link}"
                )
                try:
                    post_to_slack(SLACK_WEBHOOK_URL, message, thread_ts)
                except Exception as e:
                    print(f"Failed to send Slack message for paper {i}: {str(e)}")

        # PMIDファイルをリポジトリにコミット
        if new_pmids:
            os.system("git add processed_pmids.txt")
            os.system('git commit -m "Update processed PMIDs"')
            os.system("git push origin main || echo 'Git push failed, continuing...'")

    except google_exceptions.QuotaExceeded as e:
        print(f"Quota exceeded globally: {str(e)}")
        try:
            post_to_slack(SLACK_WEBHOOK_URL, "Gemini API quota exceeded, skipping summarization.")
        except Exception as e2:
            print(f"Failed to send Slack message for quota error: {str(e2)}")
        skip_summary = True
        # 既存の論文データで処理を継続
        if not thread_ts:
            try:
                thread_ts = post_to_slack(SLACK_WEBHOOK_URL, "This Week's Papers (6 articles):")
            except Exception as e:
                print(f"Failed to send initial Slack message after quota error: {str(e)}")
                thread_ts = None

        for i, meta in enumerate(papers_with_metadata[:6], 1):
            paper = meta["paper"]
            group = meta["group"]
            keywords = ", ".join(meta["keywords"])
            title = paper["title"]
            abstract = paper["abstract"]
            link = f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/"

            title_ja = translate_to_japanese(title)
            abstract_ja = translate_to_japanese(abstract)
            abstract_summary_ja = abstract_ja[:100] + "..."

            message = (
                f"**{i}. {title_ja} ({title})**\n"
                f"Group: {group}, Keywords: {keywords}\n"
                f"要約（簡潔）: {abstract_summary_ja}\n"
                f"要約（詳細）: {abstract_ja}\n"
                f"原文: {abstract}\n"
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
