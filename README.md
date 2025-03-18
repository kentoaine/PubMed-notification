# PubMed notification


## 概要
PubMed-notificationは、PubMedから特定のキーワードに基づいて最新の論文を自動的に検索し、そのタイトルとアブストラクトを日本語に翻訳してSlackに通知するツールです。GitHub Actionsを使用してスケジュール実行が可能で、研究者やチームが最新の論文情報を効率的に把握するのに役立ちます。
## 機能
- PubMedからの論文検索: 3つのキーワードグループ（A, B, C）からランダムに選択したキーワードで論文を検索。

- 日本語翻訳: Gemini APIを使用して、タイトルとアブストラクトを日本語に翻訳。

- Slack通知: 検索結果をSlackに投稿し、チームで共有。

- 重複防止: 処理済みの論文（PMID）を記録し、同じ論文を再通知しない。

- 自動化: GitHub Actionsで毎週日曜22:00 UTC（日本時間月曜7:00）に自動実行。

## 仕様
- 各キーワードグループ（A, B, C）から1論文ずつ、合計3論文を通知。

- 翻訳されたタイトルとアブストラクトのみを表示

- 処理済みのPMIDをリポジトリに保存し、Gitで管理。


## キーワード
|グループ|キーワード|
|----|----|
|A|semi-nested PCR, mitochondrial DNA, amplification, sequencing, phylogenetics, species identification, genetic diversity, haplotype|
|B| bioinformatics, sequence alignment, phylogenetic tree, population genetics, comparative genomics, haplotype network, genome annotation|
|C|Yamanaka factors, rejuvenation, epigenetic reprogramming, iPS, cellular senescence, aging reversal, transcription factor reprogramming, epigenetic clocks,regenerative medicine|
