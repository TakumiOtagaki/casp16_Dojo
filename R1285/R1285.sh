# to run this script;
    # bash R1203.sh
# in any directory

. ~/.bashrc
ml rosetta
# id 
export CASP_TARGET=R1285.LL
export problem_ID=2
workdir=/large/otgk/casp/casp16/R1285
cd $workdir

n=5
echo "rna_denovo"
rna_denovo.default.linuxgccrelease \
 -sequence "ggugcaguauucuagucagggaaaugcuuuuugaaggcggggcuaaaaauccgcuaaagggcacaucgaugaaguuccuggugcuggccuuagaaugcccagucuugggcuugugcugggaguuaaaaaagcuggggcacucgcaauggcaugcgacaaaugacccuacuuuuguggaggccaauuauuguauauugagagagauauucaauauacgaaauugggguaaaccugcaaugugguguaaaagcuaugugcaguguagccugccuugagugguauggggagaggagauaaacaagucaaaaauuuuaggccuaaguuuuuguacuauugaacucugaaaccuauguugcaaaagaggcuaagaaagcaucuaacuguugaggaaaacuccuagacuguuuugguaaaaugaggauugcagugcggacuuaguggcaauucaguccugaaaguggcaacacuucagcucggauauuaaagggaaaccgcuauauggcgacguauaguuauucguggggaaagccuacugaaccuaugccguaagauuuacuuauuuuguuaccacauugccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccggugcaguauucuagucagggaaaugcuuuuugaaggcggggcuaaaaauccgcuaaagggcacaucgaugaaguuccuggugcuggccuuagaaugcccagucuugggcuugugcugggaguuaaaaaagcuggggcacucgcaauggcaugcgacaaaugacccuacuuuuguggaggccaauuauuguauauugagagagauauucaauauacgaaauugggguaaaccugcaaugugguguaaaagcuaugugcaguguagccugccuugagugguauggggagaggagauaaacaagucaaaaauuuuaggccuaaguuuuuguacuauugaacucugaaaccuauguugcaaaagaggcuaagaaagcaucuaacuguugaggaaaacuccuagacuguuuugguaaaaugaggauugcagugcggacuuaguggcaauucaguccugaaaguggcaacacuucagcucggauauuaaagggaaaccgcuauauggcgacguauaguuauucguggggaaagccuacugaaccuaugccguaagauuuacuuauuuuguuaccacauugcc" \
 -secstruct "(((.....................((((((((...((((((........))))))))))))))...........[[[[[[[[[[.[[[...{{{....[[((}}}.]]]]].]]]]]]]]]]......[[.)).......(((((......)))))...........((]]..)).....((((((.((((((((((((.......)))))))))))))))))){{...{{{{{{{(((([[[[........)).}}}}}.}}..}}..((....((..((((((....((((.......((((((((((.(((((((((((((((((.(((......((((....(((((...................................((((....)))).................................(((..............)))(((.(((((....)))))))).............((....))((((((((....)))))))).....(((((.....)))))................................[[[[[[[[[[[.............................................................................................................................((((((((...((((((........))))))))))))))...........[[[[[[[[[[.[[[...{{{....[[((}}}.]]]]].]]]]]]]]]].........)).......(((((......)))))...........((....)).....((((((.((((((((((((.......)))))))))))))))))).......((]]]]]]]]]]].....{{{{...[[))..}}}}..]]..........)))))....)))).......)))))))))).))))))))))))))))).)))......))))....))))))..))......))......................((((....)))).................................(((..............)))(((.(((((....)))))))).............((....))((((((((....)))))))).....(((((.....)))))................................]]]]..)).)))" \
 -nstruct ${n} -out:file:silent farfar2_result/${CASP_TARGET}n${n}_nomini.out \
 -minimize_rna false \
 -no_filters true \
 -staged_constraints false \
 -filter_lores_base_pairs false \
 -filter_chain_closure false \
 -output_lores_silent_file false \
 -close_loops false \
 -autofilter false \


# rna_denovo.default.linuxgccrelease -fasta <fasta file with sequence> \
#  -use_legacy_job_distributor false \
#  -nstruct 1 -no_filters true -cycles 0\
#  -staged_constraints false -filter_lores_base_pairs \
#  false -filter_chain_closure false -output_lores_silent_file false \
#  -close_loops false -autofilter false -allow_bulge false \
#  -lores_scorefxn rna_lores_linear_chainbreak_weight 0


# .out to .pdb
# mkdir pdb
# cd pdb
echo "silent file to pdb"
rna_extract.linuxgccrelease \
 -in:file:silent farfar2_result/pdb/${CASP_TARGET}n${n}_nomini.out \
 -in:file:silent_struct_type rna \


# rna_extract.linuxgccrelease \
#  -in:file:silent /large/otgk/casp/casp16/R1285/farfar2_result/R1285.LLn1_nomini.out \
#  -in:file:silent_struct_type rna \
