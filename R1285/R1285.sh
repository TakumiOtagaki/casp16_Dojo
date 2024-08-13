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

echo "rna_denovo"
rna_denovo.default.linuxgccrelease \
 -sequence "ggugcaguauucuagucagggaaaugcuuuuugaaggcggggcuaaaaauccgcuaaagggcacaucgaugaaguuccuggugcuggccuuagaaugcccagucuugggcuugugcugggaguuaaaaaagcuggggcacucgcaauggcaugcgacaaaugacccuacuuuuguggaggccaauuauuguauauugagagagauauucaauauacgaaauugggguaaaccugcaaugugguguaaaagcuaugugcaguguagccugccuugagugguauggggagaggagauaaacaagucaaaaauuuuaggccuaaguuuuuguacuauugaacucugaaaccuauguugcaaaagaggcuaagaaagcaucuaacuguugaggaaaacuccuagacuguuuugguaaaaugaggauugcagugcggacuuaguggcaauucaguccugaaaguggcaacacuucagcucggauauuaaagggaaaccgcuauauggcgacguauaguuauucguggggaaagccuacugaaccuaugccguaagauuuacuuauuuuguuaccacauugccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccggugcaguauucuagucagggaaaugcuuuuugaaggcggggcuaaaaauccgcuaaagggcacaucgaugaaguuccuggugcuggccuuagaaugcccagucuugggcuugugcugggaguuaaaaaagcuggggcacucgcaauggcaugcgacaaaugacccuacuuuuguggaggccaauuauuguauauugagagagauauucaauauacgaaauugggguaaaccugcaaugugguguaaaagcuaugugcaguguagccugccuugagugguauggggagaggagauaaacaagucaaaaauuuuaggccuaaguuuuuguacuauugaacucugaaaccuauguugcaaaagaggcuaagaaagcaucuaacuguugaggaaaacuccuagacuguuuugguaaaaugaggauugcagugcggacuuaguggcaauucaguccugaaaguggcaacacuucagcucggauauuaaagggaaaccgcuauauggcgacguauaguuauucguggggaaagccuacugaaccuaugccguaagauuuacuuauuuuguuaccacauugcc" \
 -secstruct "(((.....................((((((((...((((((........))))))))))))))...........[[[[[[[[[[.[[[...{{{....[[((}}}.]]]]].]]]]]]]]]]......[[.)).......(((((......)))))...........((]]..)).....((((((.((((((((((((.......)))))))))))))))))){{...{{{{{{{(((([[[[........)).}}}}}.}}..}}..((....((..((((((....((((.......((((((((((.(((((((((((((((((.(((......((((....(((((...................................((((....)))).................................(((..............)))(((.(((((....)))))))).............((....))((((((((....)))))))).....(((((.....)))))................................[[[[[[[[[[[.............................................................................................................................((((((((...((((((........))))))))))))))...........[[[[[[[[[[.[[[...{{{....[[((}}}.]]]]].]]]]]]]]]].........)).......(((((......)))))...........((....)).....((((((.((((((((((((.......)))))))))))))))))).......((]]]]]]]]]]].....{{{{...[[))..}}}}..]]..........)))))....)))).......)))))))))).))))))))))))))))).)))......))))....))))))..))......))......................((((....)))).................................(((..............)))(((.(((((....)))))))).............((....))((((((((....)))))))).....(((((.....)))))................................]]]]..)).)))" \
 -nstruct 50 -out:file:silent farfar2_result/${CASP_TARGET}.out \
 -minimize_rna true


# .out to .pdb
# mkdir pdb
# cd pdb
echo "silent file to pdb"
rna_extract.linuxgccrelease \
 -in:file:silent farfar2_result/pdb/${CASP_TARGET}.out \
 -in:file:silent_struct_type rna \
