#!/bin/bash

# 設定
. ~/.bashrc
module load rosetta

cd /large/otgk/casp/casp16/utils/examples

INPUT_FASTA=/large/otgk/casp/casp16/1_R1203/1.fasta
INPUT_PDB=./R1205_FF2_S000001_with_phosphate.pdb
CUSTOM_TORSIONS=./R1205_FF2_S000001_with_phosphate.torsions
# OUTPUT_PDB=output.pdb

# フラグメントライブラリの作成
echo "TORSIONS extraction"
rna_denovo.default.linuxgccrelease  -sequence "gcccggauagcucagucgguagagcagcgggcacuaugggcgcagugucaauggacgcugacgguacaggccagacaauuauugucugguauagugcccgcggguccaggguucaagucccuguucgggcgcca" \
 -vall_torsions  -o $CUSTOM_TORSIONS -native $INPUT_PDB

# RNAデノボモデリングの実行
# rna_denovo.default.linuxgccrelease \
#     -fasta $INPUT_FASTA \
#     -in:file:silent $INPUT_PDB \
#     -in:file:frag3 $CUSTOM_TORSIONS \
#     -in:file:frag9 $CUSTOM_TORSIONS \
#     -out:file:silent_struct_type binary \
#     -out:file:silent $OUTPUT_PDB
