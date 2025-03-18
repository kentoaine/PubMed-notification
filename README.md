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


## 必要なアカウントとAPI
このプロジェクトを実行するには、以下のアカウントとAPIキーが必要です。
1. GitHubアカウント
- 用途: リポジトリをホストし、GitHub Actionsで自動実行するため。
- 取得方法:
GitHubにアクセスし、アカウントを作成。
新しいリポジトリ（例: PubMed-notification）を作成。
プライベートリポジトリ推奨（機密情報の保護のため）。
Personal Access Token（PAT）
- 用途: GitHub Actionsがリポジトリに書き込み（git push）を行うための認証トークン。
- 取得方法:
GitHubの右上にあるプロフィールアイコン > 「Settings」 > 「Developer settings」 > 「Personal access tokens」 > 「Tokens (classic)」 > 「Generate new token」。
- トークンの権限を設定:
repo（リポジトリへのフルアクセス）を選択。
トークンを生成し、コピー（例: ghp_XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX）。
GitHubリポジトリの「Settings」 > 「Secrets and variables」 > 「Actions」 > 「Secrets」に追加:
Name: PAT
Value: 生成したトークン
2. PubMed API（Entrez E-utilities）
- 用途: PubMedから論文情報を取得するため。
アカウント: 不要（ただし、メールアドレスの登録が必要）。
- 設定:
NCBIのポリシーにより、Entrez E-utilitiesを使用する際はメールアドレスを指定する必要があります。
このメールアドレスは環境変数EMAIL_ADDRESSとして設定します
過剰なリクエストを防ぐため、スクリプト内でリクエスト間にtime.sleep(2)を挿入しています。
3. Gemini API
- 用途: タイトルとアブストラクトを日本語に翻訳するため。
- アカウントとAPIキーの取得方法:
Google Cloud Consoleにアクセスし、プロジェクトを作成。
「APIs & Services」 > 「Library」から「Generative AI API」を検索し、有効化。
「APIs & Services」 > 「Credentials」 > 「Create Credentials」 > 「API Key」を選択してAPIキーを作成。
生成されたAPIキー（例: sk-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX）をコピー。
GitHubリポジトリの「Settings」 > 「Secrets and variables」 > 「Actions」 > 「Secrets」に追加:
Name: GEMINI_API_KEY
Value: 生成したAPIキー
- 注意
Gemini APIには無料枠のリクエスト制限があります（例: 60リクエスト/分）。制限を超えると翻訳がスキップされる可能性があります。
制限が頻繁に発生する場合、有料プランを検討してください。
4. SlackアカウントとWebhook URL
- 用途: 論文情報をSlackチャンネルに通知するため。
- アカウントとWebhook URLの取得方法:
Slackのワークスペースにアクセス。
「Apps」 > 「Incoming Webhooks」を検索してインストール。
通知したいチャンネルを選択し、「Add Incoming Webhooks integration」をクリック。
生成されたWebhook URL（例: https://hooks.slack.com/services/TXXXX/BXXXX/XXXXXXXXXXXXXXXXXXXXXXXX）をコピー。
GitHubリポジトリの「Settings」 > 「Secrets and variables」 > 「Actions」 > 「Secrets」に追加:
Name: SLACK_WEBHOOK_URL
Value: 生成したWebhook URL
- 注意
Webhook URLが無効化された場合（例: チャンネルが削除された場合）、新しいWebhook URLを作成する必要があります。
