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

# python env
conda activate rosetta # option


# ------------------------------------------------------- plz edit here -------------------------------------------------------
#                                                                                                                             #              
# environment settings
ROSETTA3=/large/otgk/app/rosetta/v2024.15/source 
rna_denovo_path=$ROSETTA3/bin/rna_denovo.default.linuxgccrelease
rna_extract_path=$ROSETTA3/bin/rna_extract.default.linuxgccrelease

# otagaki's repo directory
repo_dir=/large/otgk/casp/casp16
CASP_ID=R1212 

# working dir
workdir=/large/otgk/casp/casp16/R1212 # this dir should contain unformatted/*.pdb

#                                                                                                                             #
# ------------------------------------------------------- plz edit here -------------------------------------------------------

cd $workdir 


# input validation
mkdir p_added; mkdir formatted; mkdir merged_for_submission; cd p_added
if [ ! -e $workdir/unformatted ]; then
    echo "Error: $workdir/unformatted does not exist."
    exit 1
fi



# main process
un_formatted_pdbs=($(ls $workdir/unformatted))
for pdb_file in "${un_formatted_pdbs[@]}"; do
    # 拡張子を削除
    pdb_base=$(basename "$pdb_file" .pdb) # structure's id.
    echo "Processing file: $pdb_base"

    python $repo_dir/utils/add_residue_to_rna_.py -pdb ${workdir}/unformatted/${pdb_base}.pdb \
     -f ${workdir}/R1212.fa -ss ${workdir}/R1212.secstruct \
     -o ${workdir}/p_added/${pdb_base}.out -r a \
     -n 1 \
     --rna_extract_path $rna_extract_path \
     --rna_denovo_path $rna_denovo_path \
     --rosetta3 $ROSETTA3
    python $repo_dir/utils/rna_formatter.py -i ${workdir}/p_added/S_000001.pdb -o ${workdir}/formatted/${pdb_base}.formatted.pdb
done
rm tmp.pdb

# merging
python  $repo_dir/utils/merge_formatted_5structure.py -i ${workdir}/formatted -o ${workdir}/merged_for_submission/$CASP_ID.merged.formatted.pdb -t ${CASP_ID}
echo "Done: formatted & merged pdb is created in $workdir/merged_for_submission/$CASP_ID.merged.formatted.pdb"


