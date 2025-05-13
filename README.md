# TPE_search

This project uses [uv](https://docs.astral.sh/uv/) for python version management.

## WORK_DIR内のディレクトリ・ファイル構造
{WORK_DIR}
├── run_optuna_py.sh
├── run_vasp6.4.3.sh
├── optuna_ad_search.py
├── INCAR
├── KPOINTS
├── POTCAR (表面の原子+分子の原子の順番(ただし、同種原子をまとめない))
├── .env
├── mole.vasp (.envに合わせる)
└── sur.vasp (.envに合わせる)

## 実行方法
{WORK_DIR}内で`qsub  run_optuna_py.sh`

## surやmole.vaspの条件
- 相対座標
- セルは同じものを使用 (動作が安定するため)
- moleは全原子がセル内に配置されるようにする (重心を求めるため)

## その他の注意事項
- 必要に応じて`optuna_ad_search.py`内の`MASSDICT`に分子内の原子の原子量を追加
    - デフォルト: `MASSDICT = {"H": 1, "C": 12, "N": 14, "O": 16}`
    - スラブモデルの原子の情報は不要


## 仮想環境のアクティベートまでのパス
`/home/pj24001724/ku40000345/venv/TPE_search/.venv/bin/activate`

## ssh接続のためのRSAkey
`/home/pj24001724/ku40000345/rsa_key/busseiken_private`