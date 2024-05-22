# ABOUT Utils
## utils/add_phosphate_to_rna_.py
### 概要
`add_phosphate_to_rna.py`は、RNAの5'末端にリン酸基（P、OP1、OP2）を追加するためのPythonスクリプトです。このスクリプトは、入力PDBファイルを読み込み、指定された残基IDのリボース環のコンフォメーション（C2'-endoまたはC3'-endo）に基づいてリン酸基を適切に追加します。

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
python3 '/large/otgk/casp/casp16/utils/add_phosphate_to_rna_.py' --help
usage: add_phosphate_to_rna_.py [-h] [--initial_structure_pdb INITIAL_STRUCTURE_PDB] [--fasta FASTA] [--secondary_structure_file SECONDARY_STRUCTURE_FILE] [--fiveprime_added_out FIVEPRIME_ADDED_OUT]
                                [--adding_residue ADDING_RESIDUE] [--nstruct NSTRUCT] [--output_dir OUTPUT_DIR]

adding a residue to the RNA structure and re-running farfar2 with Rosetta, which enables us to predict the RNA tertiary structure with 5 prime phosphate.

options:
  -h, --help            show this help message and exit
  --initial_structure_pdb INITIAL_STRUCTURE_PDB
                        Path to the initial structure PDB file.
  --fasta FASTA         Path to the sequence file. The sequence should be in 'single' FASTA format.
  --secondary_structure_file SECONDARY_STRUCTURE_FILE
                        Path to the secondary structure file.format: >filename sequence secondary structure
  --fiveprime_added_out FIVEPRIME_ADDED_OUT
                        Path for the output silent file from RNA de novo.
  --adding_residue ADDING_RESIDUE
                        Residue to add to the 5' end of the sequence. if you want, you can add more than one residue. However, you should notice all the residues you selected will be attached to the 5' end of the
                        sequence. And you must use lower case.
  --nstruct NSTRUCT     Number of structures to generate.

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
python3 add_phosphate_to_rna_.py -pdb examples/rna_initial.pdb -f examples/rna.fasta -ss examples/rna.secstruct -o examples/rna.out -r a -n 2
```



#### errors you can ignore
```sh
/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:388: PDBConstructionWarning: Ignoring unrecognized record '##Begi' at line 1904
  warnings.warn(
/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:388: PDBConstructionWarning: Ignoring unrecognized record 'BINARY' at line 1905
  warnings.warn(
/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:388: PDBConstructionWarning: Ignoring unrecognized record '##End ' at line 1906
  warnings.warn(
/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:388: PDBConstructionWarning: Ignoring unrecognized record 'N_BS 5' at line 1907
  warnings.warn(
/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:388: PDBConstructionWarning: Ignoring unrecognized record 'N_NWC ' at line 1908
  warnings.warn(
/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:388: PDBConstructionWarning: Ignoring unrecognized record 'N_WC 1' at line 1909
  warnings.warn(
/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:388: PDBConstructionWarning: Ignoring unrecognized record 'score ' at line 1910
  warnings.warn(
```
このようなエラーは無視できる。pdb ファイルのコメント行を読み込んでいるだけ。

### INSTALLATION of utils/add_phosphate_to_rna_,py
#### cloning this repo
```
git clone git@github.com:TakumiOtagaki/casp16_Dojo.git
```

#### installation of rosetta
https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build


#### check your env related to rosetta
```sh
echo $ROSETTA3
/path/to/rosetta/source
```
これが空じゃなかったら上手く install できている。

#### change `utils/add_phosphate_to_rna_.py` script
`utils/add_phosphate_to_rna_.py` の16 ~ 18 行目を変更する。

```py
rna_denovo_path = "/path/to/rna_denovo.default.linuxgccrelease"
ROSETTA3 = "/path/to/rosetta/source" # これは上で確認した $ROSETTA3 の中身
```
