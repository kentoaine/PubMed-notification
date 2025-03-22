# PubMed Notification

## 仕組み

`PubMed-notification`は、PubMedから特定のキーワードに基づいて最新の論文を自動的に検索し、そのタイトルとアブストラクトを日本語に翻訳してSlackに通知するツールです。GitHub Actionsを使用して毎週月曜と木曜の朝7:00（日本時間）にスケジュール実行され、研究者やチームが最新の論文情報を効率的に把握するのに役立ちます。

### 処理の流れ
1. **キーワード選択**: 3つのキーワードグループ（A, B, C）からランダムに3つのキーワードを選び、PubMedで検索。
2. **論文検索**: PubMed API（Entrez E-utilities）を使用して論文を検索（各グループ1論文、合計3論文）。
3. **詳細取得**: 検索結果からタイトルとアブストラクトを取得。
4. **翻訳**: Gemini APIを使用してタイトルとアブストラクトを日本語に翻訳。
5. **通知**: 翻訳結果をSlackに投稿。
6. **記録**: 処理済みのPMIDを`processed_pmids.txt`に保存し、GitHubリポジトリにコミット。

### 特徴
- 翻訳された情報のみを表示（原文は非表示）。
- アブストラクトが500文字を超える場合、切り捨てて`...`を付加。
- 一時的なAPIエラーに備え、リトライロジックを実装。
- Slack通知が失敗してもプロセスが継続。

## キーワード

以下は、論文検索に使用するキーワードグループです。各グループからランダムに3つ選択されます。

| グループ | キーワード                                      |
|----------|------------------------------------------------|
| A        | semi-nested PCR, mitochondrial DNA, amplification, sequencing, phylogenetics, species identification, genetic diversity, haplotype |
| B        | bioinformatics, sequence alignment, phylogenetic tree, population genetics, comparative genomics, haplotype network, genome annotation |
| C        | Yamanaka factors, rejuvenation, epigenetic reprogramming, iPS, cellular senescence, aging reversal, transcription factor reprogramming, epigenetic clocks, regenerative medicine |

## 必要要件

このプロジェクトを実行するには、以下のアカウントとAPIキーが必要です。

### 1. GitHubアカウント
- **用途**: リポジトリをホストし、GitHub Actionsで自動実行するため。
- **取得方法**:
  - [GitHub](https://github.com/)にアクセスし、アカウントを作成。
  - 新しいリポジトリ（例: `PubMed-notification`）を作成し、プライベートに設定（機密情報の保護のため）。
- **Personal Access Token（PAT）**:
  - **用途**: GitHub Actionsがリポジトリに書き込み（`git push`）を行うための認証。
  - **取得方法**:
    1. GitHubの「Settings」 > 「Developer settings」 > 「Personal access tokens」 > 「Tokens (classic)」 > 「Generate new token」。
    2. スコープで`repo`（リポジトリへのフルアクセス）を選択。
    3. トークンを生成し、コピー（例: `ghp_XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`）。
    4. GitHubリポジトリの「Settings」 > 「Secrets and variables」 > 「Actions」 > 「Secrets」に追加:
       - **Name**: `PAT`
       - **Value**: 生成したトークン

### 2. PubMed API（Entrez E-utilities）
- **用途**: PubMedから論文情報を取得。
- **アカウント**: 不要（ただしメールアドレスが必要）。
- **設定**:
  - NCBIのポリシーにより、Entrez E-utilitiesを使用する際はメールアドレスを指定。
  - 環境変数`EMAIL_ADDRESS`として設定（詳細は後述）。
  - 過剰なリクエストを防ぐため、スクリプト内で`time.sleep(2)`を挿入。

### 3. Gemini API
- **用途**: タイトルとアブストラクトを日本語に翻訳。
- **アカウントとAPIキーの取得方法**:
  1. [Google Cloud Console](https://console.cloud.google.com/)にアクセスし、プロジェクトを作成。
  2. 「APIs & Services」 > 「Library」から「Generative AI API」を有効化。
  3. 「APIs & Services」 > 「Credentials」 > 「Create Credentials」 > 「API Key」を選択し、APIキーを作成。
  4. 生成されたAPIキー（例: `sk-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`）をコピー。
  5. GitHubリポジトリの「Secrets」に追加:
     - **Name**: `GEMINI_API_KEY`
     - **Value**: APIキー
- **注意**: 無料枠のリクエスト制限（例: 60リクエスト/分）がある。制限を超える場合、翻訳がスキップされ、必要に応じて有料プランを検討。

### 4. SlackアカウントとWebhook URL
- **用途**: 論文情報をSlackチャンネルに通知。
- **取得方法**:
  1. Slackワークスペースにアクセス。
  2. 「Apps」 > 「Incoming Webhooks」をインストール。
  3. 通知したいチャンネルを選択し、「Add Incoming Webhooks integration」をクリック。
  4. 生成されたWebhook URL（例: `https://hooks.slack.com/services/TXXXX/BXXXX/XXXXXXXXXXXXXXXXXXXXXXXX`）をコピー。
  5. GitHubリポジトリの「Secrets」に追加:
     - **Name**: `SLACK_WEBHOOK_URL`
     - **Value**: Webhook URL
- **注意**: URLが無効化された場合（例: チャンネル削除時）、新しいURLを作成。

## セットアップ手順

### 1. リポジトリのクローン
```bash
git clone https://github.com/kentoaine/PubMed-notification.git
cd PubMed-notification
```




### 2. 環境変数の設定
以下の環境変数をGitHub Secretsに設定します。
- EMAIL_ADDRESS: PubMed APIで使用するメールアドレス（例: your-email@yourdomain.com）。

- SLACK_WEBHOOK_URL: Slackに通知を送信するためのWebhook URL（例: https://hooks.slack.com/services/TXXXX/BXXXX/XXXXXXXXXXXXXXXXXXXXXXXX）。

- GEMINI_API_KEY: Gemini APIのAPIキー（例: sk-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX）。

- PAT: GitHubリポジトリへの書き込み権限を持つPersonal Access Token（例: ghp_XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX）。

- 設定手順
1. GitHubリポジトリの「Settings」 > 「Secrets and variables」 > 「Actions」 > 「Secrets」に移動します。

2. 以下のシークレットを追加します:
- EMAIL_ADDRESS


- SLACK_WEBHOOK_URL


- GEMINI_API_KEY


- PAT


3. 「Add secret」または「Update secret」をクリックして保存します。

### 3. GitHub Actionsの設定
GitHub Actionsを使用して、スケジュール実行を設定します。.github/workflows/main.ymlファイルの内容は以下の通りです。
```bash
name: Weekly PubMed Notification

on:
  schedule:
    - cron: "0 22 * * 0"  # 月曜7:00 JST
    - cron: "0 22 * * 3"  # 木曜7:00 JST
  workflow_dispatch:

jobs:
  notify:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
        with:
          token: ${{ secrets.PAT }}
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: "3.9"
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -r requirements.txt -v
      - name: Run script
        env:
          SLACK_WEBHOOK_URL: ${{ secrets.SLACK_WEBHOOK_URL }}
          GEMINI_API_KEY: ${{ secrets.GEMINI_API_KEY }}
          EMAIL_ADDRESS: ${{ secrets.EMAIL_ADDRESS }}
        run: python main.py
```

**スケジュール:** 毎週月曜と木曜の朝7:00（日本時間）に実行。

**手動実行:** 「Actions」タブから「Run workflow」で手動実行が可能です。


