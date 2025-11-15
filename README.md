# triaccel

`triaccel` は、高エネルギー天体物理で扱うトリプレット（k-クリーク）出現頻度を高速にモンテカルロ評価するための Python ライブラリです。C++17 + pybind11 で実装した計算コアを提供し、ヒストグラム/マップ描画や確率集計ユーティリティを含みます。テスト用のランナーや Jupyter ノートブックなしでも、Python だけで完結した解析ループを組めることを目指しています。

## 主な特徴
- Python から直接呼び出せる高速シミュレーター (`triaccel.simulate`)
- 1D/2D ヒスト、Hammer 投影図、トリプレット確率表などの可視化ユーティリティ (`triaccel.viz`)
- `debug=True` 時の詳細ログ出力と自動 PDF 生成
- k-クリーク（既定はトリプレット）列挙。`cluster_size` を変えてペア/クアッド以上にも対応

## インストール
1. リポジトリを取得します。
   ```bash
   git clone https://github.com/i-komae/triaccel.git
   cd triaccel
   ```
2. 仮想環境を推奨します。
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
3. `triaccel` をインストールします。
   ```bash
   pip install .
   ```
   - Python 3.10 以上と、C++17/pybind11/コンパイラ（macOS なら Xcode Command Line Tools）が必要です。
   - NumPy>=2.0 / Matplotlib>=3.7 が依存関係として入ります。
4. 動作確認（任意）: `python tutorial/main.py`
   - `tutorial/` 配下は最小構成のサンプルスクリプトです。好みのパラメータに書き換えて使ってください。

## クイックスタート
```python
import triaccel

sites = [
    {"lat": 39.3, "lon": -112.9, "zmax": 55.0, "n": 32},
    {"lat": -35.2, "lon": -69.2, "zmax": 60.0, "n": 35},
]

res = triaccel.simulate(
    M=1_000_000,
    sites=sites,
    d_deg=5.0,
    bins=(144, 72),
    return_histograms=True,
    cluster_size=3,
)

triaccel.viz.write_triplet_prob_sigma(res, n=30, outdir="results", filename="counts.npy")
triaccel.viz.count_hist(res, log=True, outdir="figs", filename="triplet_counts")
triaccel.viz.hist1d(res, what="events", outdir="figs", filename="events_1d")
triaccel.viz.hammer(res, what="triplets", outdir="figs", filename="triplets_map")
```

### `triaccel.simulate` 主要引数
- `M`: 試行回数 (32bit int 上限 2,147,483,647)。`debug=True` の場合はバックエンドが 1000 にクランプしてから実行します。
- `sites`: 各観測点の緯度 (`lat`), 経度 (`lon`), 最大天頂角 (`zmax`), 検出器数 `n`。
- `d_deg`: ペアをクリークと見なす角距離。
- `bins`: `(n_ra_bins, n_dec_bins)` あるいは単一 int。ヒスト（2D）を返したい場合に必須です。
- `return_histograms`: Hammer 図や 1D プロジェクションを作りたい場合に True。2D ヒストのみを返し、1D 分布は `triaccel.viz.hist1d` が周辺和として生成します。
- `cluster_size`: 2 でペア数、3 でトリプレット、4 以上で DFS による k-クリーク列挙。
- `debug`: True で詳細ログ + 自動 PDF 生成（後述）。

## デバッグモード (`debug=True`)
- `progress=False` になり、代わりに段階ログを stderr に出力。
- `M` は自動で `min(M, 1000)` にクランプされ、`[debug] M clamped to 1000` を 1 度だけ表示。
- 実行ごとに `log/<YYYY-MM-DDTHH-MM-SSZ>/` が作成され、`text/`（各試行のイベント/クリーク TSV）と `fig/`（Hammer PDF）が入ります。
  - `text/sim_00001_events.txt`: `# idx\tsite\tra_deg\tdec_deg\tin_cluster`
  - `text/sim_00001_cliques_k3.txt`: `# v0\tv1\tv2` （k=2 はペア、k>=4 は v0..v{k-1}）
  - `fig/sim_00001_hammer.pdf`: イベントを半径 2.5° の円で描画。赤がクリーク参加、黒が非参加。

## 生成物と保存先
- `figs/`: 通常実行時の図（ヒスト/マップ）。サンプルスクリプトはここに PDF を保存します。
- `results/`: `write_triplet_prob_sigma` で出力した counts / 集計テーブル（ディレクトリ指定可）。
- `log/<timestamp>/text/` & `fig/`: デバッグ時の生ログと自動描画物。

## 可視化ユーティリティ (`triaccel.viz`)
- `count_hist(res, log=True, ...)`: トリプレット数分布。
- `hist1d(res, what="events"|"triplets", ...)`: RA/Dec の 1D ヒスト（2D データを周辺和で圧縮）。
- `hammer(res, ...)`: Hammer 投影図。RA=±180° を跨ぐ bin は自動で反転/分割し、破綻しないよう処理します。
- `write_triplet_prob_sigma(res, n, ...)`: `P(triplets >= n)` と片側ガウス σ をまとめて保存/表示。

## FAQ / トラブルシュート
- `M` が大きすぎる → C++ 側は 32bit int。`debug=False` では Python で `ValueError` を即座に出します。`debug=True` では内部でクランプします。
- macOS でビルドに失敗 → Xcode Command Line Tools, C++17, pybind11 が揃っているか確認してください。`pip install -e .` のログにヒントが出ます。
- Hammer 図で円が崩れる/赤くならない → RA=±π 付近は自動で線分を分割。クリーク色付けは `events` の `in_cluster` と `cliques_k*.txt` を突き合わせ、不一致時は警告＋和集合で赤にします。
- Ctrl+C で止めた → 途中結果は破棄され、`KeyboardInterrupt` が返ります。

## ディレクトリ構成（抜粋）
- `triaccel/` …… ライブラリ本体（pybind11 バインディング + 可視化）
- `tutorial/main.py` …… 最小限のサンプルスクリプト
- `README.md` …… 本ドキュメント
- `pyproject.toml` / `setup.py` …… ビルド設定
- `figs/`, `results/`, `log/` …… 実行時に生成される成果物（必要に応じて作成）

より詳細な API は `triaccel/__init__.py` と `triaccel/viz.py` の docstring を参照してください。
