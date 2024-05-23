# conda activate rosetta
# conda config --add channels https://conda.rosettacommons.org
# conda install pyrosetta


# To run this script;
    # python3 '/large/otgk/casp/casp16/utils/add_phosphate_to_rna_.py' --initial_structure_pdb '/large/otgk/casp/casp16/utils/examples/FF2_R1205_S_000001.pdb' --fasta '/large/otgk/casp/casp16/utils/examples/R1205.fasta' --secondary_structure_file '/large/otgk/casp/casp16/utils/examples/R1205.secstruct' --fiveprime_added_out /large/otgk/casp/casp16/utils/examples/phosphated_R1205_S_000001.out

from Bio.PDB import PDBParser, PDBIO, Select
import os
import subprocess
from modules.myutils import read_ss_file, read_singlefasta
import argparse

# --------------------------------------------- plz rewrite here!! -----------------------------------------------------
rna_denovo_path = "/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease"
ROSETTA3 = "/large/otgk/app/rosetta/v2024.15/source"
# ----------------------------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="adding a residue to the RNA structure and re-running farfar2 with Rosetta, which enables us to predict the RNA tertiary structure with 5 prime phosphate.")
    parser.add_argument("--initial_structure_pdb", "-pdb",
                        help="Path to the initial structure PDB file.")
    parser.add_argument("--fasta", "-f",
                        help="Path to the sequence file. The sequence should be in 'single' FASTA format.")
    parser.add_argument("--secondary_structure_file", "-ss",
                        help="Path to the secondary structure file.format:\n>filename\nsequence\nsecondary structure")
    # parser.add_argument("--output_structure_pdb", help="Path for the output modified structure PDB file.")
    parser.add_argument("--fiveprime_added_out", "-o",
                         help="Path for the output silent file from RNA de novo.")
    # parser.add_argument("--5prime_added.pdb", help="Path for the final structure PDB file.")
    parser.add_argument("--adding_residue", "-r",
                        type=str, default="a", help="Residue to add to the 5' end of the sequence. if you want, you can add more than one residue. However, you should notice all the residues you selected will be attached to the 5' end of the sequence. And you must use lower case.")
    parser.add_argument("--nstruct", "-n", type=int, default=1, help="Number of structures to generate.")
    # parser.add_argument("--output_dir", "-d", type=str, default="./", help="Directory to move output files to. This is because Rosetta generates output files in the current directory. You should execute this script in the directory where you want to store the output files since error may cause when the same name files are already exist in the directory.")
    args = parser.parse_args()

    if not os.path.exists(args.initial_structure_pdb):
        print("Error: initial structure PDB file not found.")
        return
    if not os.path.exists(args.fasta):
        print("Error: sequence file not found.")
        return
    if not os.path.exists(args.secondary_structure_file):
        print("Error: secondary structure file not found.")
        return
    if args.fiveprime_added_out == "" or args.fiveprime_added_out is None:
        print("Error: output silent file not found.")
        return
    # if not os.path.exists(args.output_structure_pdb):
    #     print("Error: output structure PDB file not found.")
    #     return
    # if not os.path.exists(args.output_dir):
    #     print("Error: output directory not found.")
    #     return
    if args.adding_residue != args.adding_residue.lower():
        print("Error: adding residue must be lower case.")
        return
    if args.nstruct < 1:
        print("Error: nstruct must be greater than or equal to 1.")
        return
    
    return args



def run_rna_denovo(modifiedStructurePDB, extendedSequence, modifiedss, nstruct,  rna_denovo_path, finalStructureOut):
    if len(extendedSequence) != len(modifiedss):
        print("Error: length of sequence and secondary structure does not match.")
        return
    
    cmd = f'{rna_denovo_path} \
         -s {modifiedStructurePDB} \
         -sequence "{extendedSequence}" \
         -secstruct "{modifiedss}" \
         -nstruct {nstruct} \
        -out:file:silent {finalStructureOut} \
        -minimize_rna true \
        -show_all_fixes'
        
    print(f"\n\ncmd: {cmd}\n\n\n")
    env = os.environ.copy()
    env["ROSETTA3"] = ROSETTA3
    subprocess.run(cmd, shell=True, env=env)

    # .out to pdb
    cmd = f"rna_extract.linuxgccrelease -in:file:silent {finalStructureOut}  -in:file:silent_struct_type rna"
    subprocess.run(cmd, shell=True, env=env)

class IncrementResidueNumbers(Select):
    """ 残基番号を逆順でインクリメントするクラス """
    def __init__(self, increment):
        self.increment = increment

    def accept_residue(self, residue):
        # 最初に大きな番号に変更して、その後目的の番号に減算
        new_id = residue.id[1] + self.increment
        residue.id = (residue.id[0], new_id, residue.id[2])
        return 1
    
def renumber_residues(input_file, output_file, increment):
    parser = PDBParser()
    structure = parser.get_structure('PDB', input_file)
    
    # 残基を逆順で処理するために、リストを逆転
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            residues.reverse()
            for residue in residues:
                residue.id = (residue.id[0], residue.id[1] + increment, residue.id[2])

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

def main():
    args = parse_args()

    # Data loading
    _, _seq, initial_secondary_structure = read_ss_file(args.secondary_structure_file)
    seq = read_singlefasta(args.fasta)

    # Adding one residue ('a') to the pose sequence
    adding = args.adding_residue
    extended_sequence = adding + seq 
    modified_ss = "." * len(adding) + initial_secondary_structure

    # tmp_pdb = args.output_dir + "tmp.pdb" should deal with "/"
    tmp_pdb  = "tmp.pdb"

    # Append residue to the modified structure and renumber
    renumber_residues(args.initial_structure_pdb, tmp_pdb, len(adding))

    # Run RNA de novo
    run_rna_denovo(tmp_pdb, extended_sequence, modified_ss, args.nstruct, rna_denovo_path, args.fiveprime_added_out)




if __name__ == "__main__":
    main()



