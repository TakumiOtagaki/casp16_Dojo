# to run:
# python add_phosphate_to_rna.py --input_pdb input.pdb --output_pdb output_with_phosphate.pdb

cd utils/examples
script=../add_phosphate_to_rna.py
input_pdb=./AF3FF2_hybrid_S_000001.pdb
output_with_phosphate=./AF3FF2_hybrid_S_000001_with_phosphate.pdb

python3 $script --input_pdb $input_pdb --output_pdb $output_with_phosphate


input_pdb=./R1205_FF2_S000001.pdb
output_with_phosphate=./R1205_FF2_S000001_with_phosphate.pdb

python3 $script --input_pdb $input_pdb --output_pdb $output_with_phosphate