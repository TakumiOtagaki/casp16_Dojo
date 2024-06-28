#! /bin/bash
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
# conda activate rosetta # option
python=/home/otgk/.conda/envs/rosetta/bin/python


# ------------------------------------------------------- plz edit here -------------------------------------------------------
#                                                                                                                             #              
# environment settings
ROSETTA3=/large/otgk/app/rosetta/v2024.15/source 
rna_denovo_path=$ROSETTA3/bin/rna_denovo.default.linuxgccrelease
rna_extract_path=$ROSETTA3/bin/rna_extract.default.linuxgccrelease

# otagaki's repo directory
repo_dir=/large/otgk/casp/casp16
CASP_ID=M1221

# working dir
workdir=/large/otgk/casp/casp16/M1221/unformatted_mixed/1-zdock.S_000461-top1 # this dir should contain unformatted/*.pdb

#                                                                                                                             #
# ------------------------------------------------------- plz edit here -------------------------------------------------------

cd $workdir 


# input validation
mkdir p_added; mkdir formatted; mkdir merged_for_submission; cd p_added
if [ ! -e $workdir/unformatted ]; then
    echo "Error: $workdir/unformatted does not exist."
    exit 1
fi


rna_id_list=(2 3)
# main process
for rna_id in "${rna_id_list[@]}"; do
    echo "rna_id: $rna_id"
    un_formatted_pdb=$workdir/unformatted/rna$rna_id.pdb

    # 拡張子を削除
    pdb_base=$(basename "$pdb_file" .pdb) # structure's id.
    echo "Processing file: $pdb_base"

    $python $repo_dir/utils/add_residue_to_rna_.py \
        -pdb $un_formatted_pdb \
        -f /large/otgk/casp/casp16/M1221/rna$rna_id.fa \
        -ss /large/otgk/casp/casp16/M1221/rna$rna_id.secstruct \
        -o ${workdir}/p_added/${pdb_base}.out \
        -r a \
        -n 1 \
        --rna_extract_path $rna_extract_path \
        --rna_denovo_path $rna_denovo_path \
        --rosetta3 $ROSETTA3
    $python $repo_dir/utils/rna_formatter.py -i ${workdir}/p_added/S_000001.pdb -o ${workdir}/formatted/${pdb_base}.formatted.pdb

    echo "format validation:"
    $python $repo_dir/utils/rna_format_validation.py -i ${workdir}/formatted/rna$rna_id.formatted.pdb
    rm tmp.pdb
done
# un_formatted_pdbs=($(ls $workdir/unformatted))
# for pdb_file in "${un_formatted_pdbs[@]}"; do
#     # 拡張子を削除
#     pdb_base=$(basename "$pdb_file" .pdb) # structure's id.
#     echo "Processing file: $pdb_base"

#     python $repo_dir/utils/add_residue_to_rna_.py -pdb ${workdir}/unformatted/${pdb_base}.pdb \
#      -f ${workdir}/${CASP_ID}.fa -ss ${workdir}/${CASP_ID}.secstruct \
#      -o ${workdir}/p_added/${pdb_base}.out -r a \
#      -n 1 \
#      --rna_extract_path $rna_extract_path \
#      --rna_denovo_path $rna_denovo_path \
#      --rosetta3 $ROSETTA3
#     python $repo_dir/utils/rna_formatter.py -i ${workdir}/p_added/S_000001.pdb -o ${workdir}/formatted/${pdb_base}.formatted.pdb

#     echo "format validation:"
#     python $repo_dir/utils/rna_format_validation.py -i ${workdir}/formatted/${pdb_base}.formatted.pdb
# done
# rm tmp.pdb

# merging
# python  $repo_dir/utils/merge_formatted_5structure.py -i ${workdir}/formatted -o ${workdir}/merged_for_submission/$CASP_ID.merged.formatted.pdb -t ${CASP_ID}
# echo "Done: formatted & merged pdb is created in $workdir/merged_for_submission/$CASP_ID.merged.formatted.pdb"


