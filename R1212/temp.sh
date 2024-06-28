# ID=000031
. ~/.bashrc

#$workdir
# ├── formatted <- intermediate2
# ├── merged_for_submission <- final output
# ├── p_added <- intermediate
# ├── R1212.fa <- fastafile in lowercase
# ├── R1212.secstruct <- secondary structure: 1stline: >seqname, 2ndline: sequence, 3rdline: secondary structure
# ├── job.sh
# └── unformatted <- input
#     ├── S_000031.pdb
#     ├── S_000128.pdb
#     ├── S_000321.pdb
#     ├── S_000427.pdb
#     └── S_000781.pdb


# settings
rna_denovo_path=/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease
rna_extract_path=/large/otgk/app/rosetta/v2024.15/source/bin/rna_extract.default.linuxgccrelease
ROSETTA3=/large/otgk/app/rosetta/v2024.15/source
conda activate rosetta # option
repo_dir=/large/otgk/casp/casp16
# python=/home/otgk/.conda/envs/rosetta/bin/python # option

# working dir
workdir=/large/otgk/casp/casp16/R1212 # --plz edit this path--
# もしここに unformatted/*.pdb がなかったらエラー
if [ ! -e $workdir/unformatted ]; then
    echo "Error: $workdir/unformatted does not exist."
    exit 1
fi


cd $workdir 

# casp target ID 
CASP_ID=R1212 # -- plz edit this ID --


# the header of the merged files 




mkdir p_added; mkdir formatted; mkdir merged_for_submission
cd p_added

# ディレクトリ内のファイルをリストアップ
un_formatted_pdbs=($(ls $workdir/unformatted))

# リスト内の各ファイルに対して処理を実行
for pdb_file in "${un_formatted_pdbs[@]}"; do
    # 拡張子を削除
    pdb_base=$(basename "$pdb_file" .pdb) # structure's id.
    echo "Processing file: $pdb_base"

    python $repo_dir/utils/add_residue_to_rna_.py -pdb ${workdir}/unformatted/${pdb_base}.pdb -f ${workdir}/R1212.fa -ss ${workdir}/R1212.secstruct -o ${workdir}/p_added/${pdb_base}.out -r a -n 1
    python $repo_dir/utils/rna_formatter.py -i ${workdir}/p_added/S_000001.pdb -o ${workdir}/formatted/${pdb_base}.formatted.pdb
done
rm tmp.pdb

# merge
python  $repo_dir/utils/merge.py -i ${workdir}/formatted -o ${workdir}/merged_for_submission/$CASP_ID.merged.formatted.pdb -t ${CASP_ID}



