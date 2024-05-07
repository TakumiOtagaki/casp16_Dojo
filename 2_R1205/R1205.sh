. ~/.bashrc
ml rosetta
# id 
export CASP_TARGET=R1205
cd /large/otgk/casp/casp16/2_${CASP_TARGET}

echo "rna_denovo"
rna_denovo.default.linuxgccrelease \
 -sequence "aaguacccuccaagcccuacagguuggaagagggggcuaucaguccuguaggcagacuc" \
 -secstruct "..((.(((((.....[[[[[[[[......))))).)).......]]]]]]]]......." \
 -nstruct 1000 -out:file:silent ${CASP_TARGET}.out \
 -minimize_rna true \

# .out to .pdb
echo "silent file to pdb"
rna_extract.linuxgccrelease \
 -in:file:silent ${CASP_TARGET}.out \
 -in:file:silent_struct_type rna \