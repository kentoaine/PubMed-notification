# PubMed notification

このプロジェクトは定義済みのキーワードグループに基づいてPubMedから最新の論文を自動的に取得し、毎週Slackチャンネルに投稿する。スケジューリングにはPythonとGitHub Actionsを使っている。

## 機能
- 3つのグループ(A, B, C)からキーワードを使ってPubMedを検索。
- 各グループ2件（計6件）の論文を日付順に検索。
- タイトル、アブストラクト（500文字以内）、リンクを英語でSlackに投稿。
- 毎週GitHub Actionsで実行（毎週月曜9:00UTC）。
- 失敗時のSlack通知によるエラー処理を含む。

## キーワード
|グループ|キーワード|
|----|----|
|A|semi-nested PCR, mitochondrial DNA, amplification, sequencing, phylogenetics, species identification, genetic diversity, haplotype|
|B| "bioinformatics", "sequence alignment", "phylogenetic tree", "population genetics", "comparative genomics", "haplotype network", "genome annotation"|
|C|"Yamanaka factors", "rejuvenation", "epigenetic reprogramming", "iPS", "cellular senescence", "aging reversal", "transcription factor reprogramming", "epigenetic clocks","regenerative medicine"|
