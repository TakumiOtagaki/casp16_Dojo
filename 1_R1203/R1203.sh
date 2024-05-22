# to run this script;
    # bash R1203.sh
# in any directory

. ~/.bashrc
ml rosetta
# id 
export CASP_TARGET=R1203
export problem_ID=1
# move to the directory where the script is located
cd /large/otgk/casp/casp16/${problem_ID}_${CASP_TARGET}

# rna_denovo
echo "rna_denovo"
rna_denovo.default.linuxgccrelease \
 -sequence "gcccggauagcucagucgguagagcagcgggcacuaugggcgcagugucaauggacgcugacgguacaggccagacaauuauugucugguauagugcccgcggguccaggguucaagucccuguucgggcgcca" \
 -secstruct "(((((((..((((........)))).(((((((((((.(.(((((((((....))))))).)).)....(((((((((...)))))))))))))))))))).....(((((.......))))))))))))...." \
 -nstruct 1000 -out:file:silent ${CASP_TARGET}.out \
 -minimize_rna true \

# .out to .pdb
echo "silent file to pdb"
rna_extract.linuxgccrelease \
 -in:file:silent ${CASP_TARGET}.out \
 -in:file:silent_struct_type rna \



