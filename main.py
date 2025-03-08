import random
import requests
from Bio import Entrez
import time

# キーワード群
keyword_groups = {
    "A": ["semi-nested PCR", "mitochondrial DNA", "amplification", "sequencing", "phylogenetics", "species identification", "genetic diversity", "haplotype"],
    "B": ["bioinformatics", "sequence alignment", "phylogenetic tree", "population genetics", "comparative genomics", "haplotype network", "genome annotation"],
    "C": ["Yamanaka factors", "rejuvenation", "epigenetic reprogramming", "iPS", "cellular senescence", "aging reversal", "transcription factor reprogramming", "epigenetic clocks","regenerative medicine"]
}

# API設定
Entrez.email = "10311kaduken@gmail.com"  # メールアドレスを入力
SLACK_WEBHOOK_URL = "https://hooks.slack.com/services/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"  # SlackのWebhook URLを入力

def search_pubmed(keyword, max_results=2):
    handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results, sort="date")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_details(pmids):
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    papers = []
    for article in records["PubmedArticle"]:
        pmid = article["MedlineCitation"]["PMID"]
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No abstract"])[0]
        # 500文字超の場合にカット（Slackの読みやすさのため）
        if len(abstract) > 500:
            abstract = abstract[:500] + "... (see link for full text)"
        papers.append({"pmid": pmid, "title": title, "abstract": abstract})
    return papers

def post_to_slack(webhook_url, message):
    payload = {"text": message}
    requests.post(webhook_url, json=payload)

def main():
    selected_keywords = {group: random.sample(keywords, 3) for group, keywords in keyword_groups.items()}
    all_papers = []
    for group, keywords in selected_keywords.items():
        for keyword in keywords:
            pmids = search_pubmed(keyword, max_results=2)
            papers = fetch_details(pmids)
            all_papers.extend(papers[:2])
            time.sleep(1)  # PubMed APIのレート制限回避

    slack_message = "This Week's Papers (6 articles):\n\n"
    for i, paper in enumerate(all_papers[:6], 1):
        title = paper["title"]
        abstract = paper["abstract"]
        link = f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/"
        slack_message += f"**{i}. {title}**\nAbstract: {abstract}\nLink: {link}\n---\n"

    post_to_slack(SLACK_WEBHOOK_URL, slack_message)

if __name__ == "__main__":
    main()