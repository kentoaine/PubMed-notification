import os
import random
import requests
from Bio import Entrez
import time

keyword_groups = {
    "A": ["semi-nested PCR", "mitochondrial DNA", "amplification", "sequencing", "phylogenetics", "species identification", "genetic diversity", "haplotype"],
    "B": ["bioinformatics", "sequence alignment", "phylogenetic tree", "population genetics", "comparative genomics", "haplotype network", "genome annotation"],
    "C": ["Yamanaka factors", "rejuvenation", "epigenetic reprogramming", "iPS", "cellular senescence", "aging reversal", "transcription factor reprogramming", "epigenetic clocks","regenerative medicine"]
}

Entrez.email = "10311kaduken@gmail.com"
SLACK_WEBHOOK_URL = os.environ.get("SLACK_WEBHOOK_URL")
if not SLACK_WEBHOOK_URL:
    raise ValueError("SLACK_WEBHOOK_URL is not set in environment variables")

def search_pubmed(keyword, max_results=2):
    try:
        handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        post_to_slack(SLACK_WEBHOOK_URL, f"PubMed search failed for '{keyword}': {str(e)}")
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
        response = requests.post(webhook_url, json=payload)
        response.raise_for_status()
        return response.json().get("ts")
    except Exception as e:
        print(f"Slack notification failed: {str(e)}")

def main():
    try:
        # 各グループから3つのキーワードを選択
        selected_keywords = {group: random.sample(keywords, 3) for group, keywords in keyword_groups.items()}
        papers_with_metadata = []

        # 各グループごとに処理
        for group, keywords in selected_keywords.items():
            group_papers = []
            group_metadata = []
            
            # 3つのキーワードで検索
            for keyword in keywords:
                pmids = search_pubmed(keyword, max_results=2)
                if pmids:
                    papers = fetch_details(pmids)
                    for paper in papers:
                        group_metadata.append({
                            "paper": paper,
                            "group": group,
                            "keyword": keyword
                        })
                time.sleep(1)

            # グループからランダムに2論文を選択
            if group_metadata:
                selected_papers = random.sample(group_metadata, min(2, len(group_metadata)))
                papers_with_metadata.extend(selected_papers)

        # Slackメッセージ作成
        if not papers_with_metadata:
            slack_message = "No papers found this week."
            post_to_slack(SLACK_WEBHOOK_URL, slack_message)
        else:
            thread_ts = post_to_slack(SLACK_WEBHOOK_URL, "This Week's Papers (up to 6 articles):")
            for i, meta in enumerate(papers_with_metadata[:6], 1):
                paper = meta["paper"]
                group = meta["group"]
                keyword = meta["keyword"]
                title = paper["title"]
                abstract = paper["abstract"]
                link = f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/"
                message = (
                    f"**{i}. {title}**\n"
                    f"Group: {group}, Keyword: {keyword}\n"
                    f"Abstract: {abstract}\n"
                    f"Link: {link}"
                )
                post_to_slack(SLACK_WEBHOOK_URL, message, thread_ts)

    except Exception as e:
        error_msg = f"System error: {str(e)}"
        post_to_slack(SLACK_WEBHOOK_URL, error_msg)
        raise

if __name__ == "__main__":
    main()
