# ABOUT Utils

## 一発で formatting とリン酸基付与をやってしまう
### `python` 環境整える
以下のコードを実行
```sh
cd $this_repo # edit here
conda create -n rosetta python=3.10.14
pip install -r utils/requirements.txt
```

### 実行
1. `utils/formatting_template.sh`をコピー
2. 25 ~ 34 行目を正しく書き換える
  - `$workdir/unformatted/` 以下に 5 つの unformatted かつ リン酸基が付いていないような pdb ファイルが入っていなくてはならない！それを入力にしてコードが動く
  - 「`$workdir/$CASP_ID.fa`: ターゲットの fasta file」 が存在しなくてはならない
    - fasta file の塩基は必ず小文字でなくてはならない
  - 「`$workdir/$CASP_ID.secstruct`：ターゲットの secondary structure file」が存在しなくてはならない
  - `$workdir/formatted` とか `$workdir/padded` とかが中間生成物としてできるが、無視して良い。
  - `$workdir/merged_for_submission/`に最終生成ファイルができる。

ちなみに `.secstruct` ファイルについては

```.secstruct file
> sequence name
塩基配列（a, t, g, c）
二次構造（ dot bracket）
```
のように3行で構成されている必要がある。

これら全ての条件が満たされていたら、

```sh
bash formatting_template.sh
```
で全て自動で終わる。




## utils/add_residue_to_rna_.py
### 概要
`add_residue_to_rna.py`は、RNAの5'末端にリン酸基（P、OP1、OP2）を追加するためのPythonスクリプトです。このスクリプトは、入力PDBファイルを読み込み、指定された残基IDのリボース環のコンフォメーション（C2'-endoまたはC3'-endo）に基づいてリン酸基を適切に追加します。

- 入力：
    - farfar2 で計算して得られる 5'末端にリン酸を持っていない構造の pdb ファイル
    - farfar2 の入力とした使った塩基配列ファイル（.fasta）
    - farfar2 の入力として扱った 二次構造ファイル（.secstruct）
- 計算
    - 与えられた塩基配列に一塩基先頭に追加する（default: "a" を追加）
    - 与えられた pdb ファイルの原子は固定するという約束（rna_denovo -s ）のもとで構造予測を再度 rna_denovo で行う
    - 二次構造は 与えられた二次構造 ss = "." + ss という形にして、先頭に追加した塩基は塩基対を形成しないものとする
    - 得られた pdb ファイルのうち 必要のない原子（先頭に追加した "a" のうち、欲しかったリン酸基以外の全て）を消去
- output
    - 入力とした pdb に何らかの塩基をつけて再度エネルギー最小化をした後にリン酸以外を消去してできた pdb ファイルのように考えて、coding しています。



### HOW TO USE
```sh
python3 '/large/otgk/casp/casp16/utils/add_residue_to_rna_.py' --help
usage: add_residue_to_rna_.py [-h] [--initial_structure_pdb INITIAL_STRUCTURE_PDB] [--fasta FASTA] [--secondary_structure_file SECONDARY_STRUCTURE_FILE] [--fiveprime_added_out FIVEPRIME_ADDED_OUT]
                              [--adding_residue ADDING_RESIDUE] [--nstruct NSTRUCT] [--rna_denovo_path RNA_DENOVO_PATH] [--rna_extract_path RNA_EXTRACT_PATH] [--rosetta3 ROSETTA3]

adding a residue to the RNA structure and re-running farfar2 with Rosetta, which enables us to predict the RNA tertiary structure with 5 prime phosphate.

options:
  -h, --help            show this help message and exit
  --initial_structure_pdb INITIAL_STRUCTURE_PDB, -pdb INITIAL_STRUCTURE_PDB
                        Path to the initial structure PDB file.
  --fasta FASTA, -f FASTA
                        Path to the sequence file. The sequence should be in 'single' FASTA format.
  --secondary_structure_file SECONDARY_STRUCTURE_FILE, -ss SECONDARY_STRUCTURE_FILE
                        Path to the secondary structure file.format: >filename sequence secondary structure
  --fiveprime_added_out FIVEPRIME_ADDED_OUT, -o FIVEPRIME_ADDED_OUT
                        Path for the output silent file from RNA de novo.
  --adding_residue ADDING_RESIDUE, -r ADDING_RESIDUE
                        Residue to add to the 5' end of the sequence. if you want, you can add more than one residue. However, you should notice all the residues you selected will be attached to the 5' end of the sequence. And
                        you must use lower case.
  --nstruct NSTRUCT, -n NSTRUCT
                        Number of structures to generate.
  --rna_denovo_path RNA_DENOVO_PATH
                        Path to the RNA de novo executable.
  --rna_extract_path RNA_EXTRACT_PATH
                        Path to the RNA extract executable.
  --rosetta3 ROSETTA3   Path to the Rosetta3 directory.

```
-  output directory: 
    - Directory to move output files cannot change from current directory. This is because Rosetta generates output files in the current directory. 

.secstruct format should be as following:
```.secstruct
>sequence_name
uuuugcccuuu
..(.....(()))...
```
- line1: sequence name
- line2: sequence
- line3: dot bracket rna ss

#### EXAMPLE
`utils/examples/` に存在する、
- 'rna.fasta'
- 'rna.secstruct'
- 'rna_initial.pdb'
に対して適用する。

```sh
cd utils
python3 add_residue_to_rna_.py -pdb examples/rna_initial.pdb -f examples/rna.fasta -ss examples/rna.secstruct -o examples/rna.out -r a -n 2 \
 --rna_extract_path $rna_extract_path \
 --rna_denovo_path $rna_denovo_path \
 --rosetta3 $ROSETTA3
```
これで OK.



### INSTALLATION of utils/add_residue_to_rna_,py
#### cloning this repo
```
git clone git@github.com:TakumiOtagaki/casp16_Dojo.git
```

#### installation of rosetta
https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build

### installation of biopython
```
conda install -c conda-forge biopython
```

#### check your env related to rosetta
```sh
echo $ROSETTA3
/path/to/rosetta/source
```
これが空じゃなかったら上手く install できている。


そしてこの $ROSETTA3 が --rosetta3 の引数であって、
```sh
python utils/add_residue_to_rna_.py ... --rosetta3 $ROSETTA3
```
のようにコールしてやる。



## utils/rna_formatter.py
### 概要
このスクリプトはPDBファイルを処理し、以下の変更を行います：

- 先頭の残基を削除（オプションで保持可能）
- 残基番号と原子番号を振り直し
- チェーンIDをアルファベットまたは数字に変更
- 水素原子を削除
- 原子名が1文字を超える場合に警告を出力
- すべての原子の占有率を1.00に設定

## utils/add_residue_to_rna_.py と utils/rna_formatter.py を組み合わせる
```sh
python3 '/large/otgk/casp/casp16/utils/rna_formatter.py' --help
usage: rna_formatter.py [-h] [--input_file INPUT_FILE] [--output_file OUTPUT_FILE] [--keep_first_residue] [--verbose]

Process and format a PDB file. The first residue is removed by default but can be kept with the --keep_first_residue option. The output file is formatted with
renumbered residues and atoms, chain IDs, and occupancy set to 1.00.

options:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE, -i INPUT_FILE
                        Path to the input PDB file.
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        Path to the output formatted PDB file.
  --keep_first_residue  Keep the first residue. By default, the first residue is removed.
  --verbose, -v         Print additional information.
```

実際に `utils/examples/S_000001.pdb` に対して formatting を行って、`utils/examples/S_000001.formatted.pdb`を作成したいとすると、
```sh
python3 'utils/rna_formatter.py' -i 'utils/examples/S_000001.pdb' -o utils/examples/S_000001.formatted.pdb 
```
とすればOK.


また、add_residue_to_rna_.py によって出力した構造以外を入力としたい場合など、
先頭残基を消去したくないのであれば、
`--keep_first_residue` とすることによって先頭の残基（残基ID = 1）のものを消去せずにそのまま放置する。P, OP2, OP1, O5' から始まっていない場合は警告を出す。




## rna_format_validation.py
format の validation を行う。
### HOW TO USE
```sh
python3 '/large/otgk/casp/casp16/utils/rna_format_validation.py' --help
usage: rna_format_validation.py [-h] [--pdb_file PDB_FILE]

Validate a PDB file against specific structural criteria.

options:
  -h, --help            show this help message and exit
  --pdb_file PDB_FILE, -i PDB_FILE
                        Path to the PDB file to validate.
```

### 検証項目
- 水素原子の非存在
  - 構造内の水素原子が削除されていることを確認します。
- 残基IDとインデックスの連続性
  - 残基番号と原子番号が1から始まり、単調増加であること（+1ずつ）を確認します。
- 元素名の長さ
  - 全原子の元素名が1文字であることを確認します。複数文字の元素名には警告を発します。
- RNAのチェーンID
  - RNAのチェーンIDが数値で0から始まることを確認します。
- 占有率の確認
  - 全原子の占有率（Occupancy）が1.00であることを確認します。
- 先頭残基の特定原子の順序
  - 先頭残基の最初の4つの原子がP, OP2, OP1, O5'の順に配置されていることを確認します。

検証して違反している点があれば warning が出る。
warning の個数も最後に表示される。