. ~/.bashrc
ml rosetta

OUT_FILE=/large/otgk/casp/casp16/1_R1203/for_merging/result_models.out
cd /large/otgk/casp/casp16/1_R1203/for_merging


echo "silent file to pdb"
rna_extract.linuxgccrelease \
 -in:file:silent $OUT_FILE \
 -in:file:silent_struct_type rna 