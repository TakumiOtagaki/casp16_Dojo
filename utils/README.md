# Utils
## add_phosphate_to_rna.py
### 概要
`add_phosphate_to_rna.py`は、RNAの5'末端にリン酸基（P、OP1、OP2）を追加するためのPythonスクリプトです。このスクリプトは、入力PDBファイルを読み込み、指定された残基IDのリボース環のコンフォメーション（C2'-endoまたはC3'-endo）に基づいてリン酸基を適切に追加します。




- --input_pdb: 入力PDBファイルのパス
- --output_pdb: 出力PDBファイルのパス
- --residue_id: 5'末端の残基ID. default では入力 pdb ファイルの最初の INDEX を採用するようになっている
- --help: help


以下のように呼び出す。
```bash
python add_phosphate_to_rna.py --input_pdb input.pdb --output_pdb output_with_phosphate.pdb
```

開始残基番号を指定したいのなら、以下のように変更する。
```bash
python add_phosphate_to_rna.py --input_pdb input.pdb --output_pdb output_with_phosphate.pdb --residue_id 1
```






### 使用例
```sh
python add_phosphate_to_rna.py --input_pdb input.pdb --output_pdb output_with_phosphate.pdb
# このコマンドは、入力PDBファイルから5'末端の残基IDを自動的に検出し、リン酸基を追加した新しいPDBファイルを出力します。

python add_phosphate_to_rna.py --input_pdb input.pdb --output_pdb output_with_phosphate.pdb --residue_id 5
# このコマンドは、指定された残基ID 5 を使用してリン酸基を追加します。
```

### スクリプトの詳細
#### 関数
- read_pdb(file_path):
PDBファイルを読み込み、原子の座標と行を返します。

- determine_ribose_conformation(atoms, resid):
指定された残基のリボース環のコンフォメーション（C2'-endoまたはC3'-endo）を判定します。

- add_phosphate(atoms, conformation, resid):
リボースのコンフォメーションに基づいてリン酸基を追加します。

- write_pdb(output_path, lines, phosphate_coords, resid):
リン酸基を追加した新しいPDBファイルを書き出します。

- check_residue_id(atoms):
残基IDが1であるかどうかを確認し、そうでない場合は警告を出し、最初の残基IDを返します。



### スクリプトの流れ
コマンドライン引数を解析します。
入力PDBファイルを読み込み、原子の座標と行を取得します。
残基IDが指定されていない場合、最初の残基IDを自動的に使用します。
指定された残基のリボース環のコンフォメーションを判定します。
リボースのコンフォメーションに基づいてリン酸基を追加します。
リン酸基を追加した新しいPDBファイルを書き出します。


### 注意点
入力PDBファイルにO5'原子が存在しない場合は、代わりにC5'原子を使用してリン酸基を追加します。この場合、警告メッセージが表示されます。
残基IDが指定されていない場合、スクリプトは自動的に最初の残基IDを使用します。


### エラーメッセージ
- add_phosphate関数で、P、OP1、OP2、O5'の4つの原子を追加するように修正しました。
- write_pdb関数で、リン酸基の原子を5'末端の残基に続けて追加するように修正しました。
- KeyError: "O5'": O5'原子が見つからない場合に発生します。このスクリプトでは、O5'原子が存在しない場合に代わりにC5'原子を使用するように修正されています。

