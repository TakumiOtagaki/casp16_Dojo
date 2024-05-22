# conda activate rosetta
# conda config --add channels https://conda.rosettacommons.org
# conda install pyrosetta


# pyrosetta database path is needed;
    # $(pip show pyrosetta)/pyrosetta/database

from pyrosetta import *
from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.constraint_generator import *
from pyrosetta.rosetta.protocols.relax import ClassicRelax

from Bio import PDB
import os
import subprocess

rna_denovo_path = "/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease"
pyrosetta_database = "/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/pyrosetta/database/"

def cif2pdb(cif_file, pdb_file):
    parser = PDB.MMCIFParser()
    structure = parser.get_structure('RNA', cif_file)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def merge_pdbs(output_pdb, *pdbs):
    with open(output_pdb, 'w') as outfile:
        for fname in pdbs:
            with open(fname) as infile:
                for line in infile:
                    if line.startswith('ATOM'):
                        outfile.write(line)
    adjust_atom_indices(output_pdb)

def adjust_atom_indices(pdb_file):
    parser = PDB.PDBParser()
    structure = parser.get_structure('Merged', pdb_file)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)


def phosphate_free_structure(pdb_file):
    pose = pose_from_pdb(pdb_file)
    print(pose.sequence())
    return pose

def extract_i_to_j_nucleotides(pdb_file, i, j, fragment_pdb_file):
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    structure = parser.get_structure('RNA', pdb_file)
    class ResidueSelect(PDB.Select):
        def accept_residue(self, residue):
            # i <= x <= j
            if residue.id[1] >= i and residue.id[1] <= j:
                return 1
            return 0
    io.set_structure(structure)
    io.save(fragment_pdb_file, select=ResidueSelect())
    print(f"Extracted nucleotides from {i} to {j} to {fragment_pdb_file}")

# Read a secondary structure file
# secstruct format:
    # >filename
    # sequence
    # secondary structure
def read_ss_file(ss_file):
    with open(ss_file, "r") as f:
        lines = f.readlines()
        filename = lines[0].strip()
        sequence = lines[1].strip()
        secstruct = lines[2].strip()
    return filename, sequence, secstruct

def data_loading(pfree_pdb_path, pfree_ss_path):
    pfree_chunk = phosphate_free_structure(pfree_pdb_path)
    _, seq, pfree_ss = read_ss_file(pfree_ss_path)
    assert pfree_chunk.sequence() == seq
    print("Data loading is done.")
    return pfree_chunk, pfree_ss

def run_rna_denovo(pfree_pdb, padded_seq, padded_ss, nstruct,  rna_denovo_path, result_models_out_path, padded_pdb):
    if len(padded_seq) != len(padded_ss):
        print("Error: length of sequence and secondary structure does not match.")
        return
    
    cmd = f'{rna_denovo_path} \
         -s {pfree_pdb} \
         -sequence "{padded_seq}" \
         -secstruct "{padded_ss}" \
         -nstruct {nstruct} \
        -out:file:silent {result_models_out_path} \
        -minimize_rna true'
    print(f"\n\ncmd: {cmd}\n\n\n")
    # os.system("source ~/.bashrc")
    # os.system("module load rosetta")
    # os.system(cmd)
    env = os.environ.copy()
    env["ROSETTA3"] = "/large/otgk/app/rosetta/v2024.15/source"
    subprocess.run(cmd, shell=True, env=env)

    # .out to pdb
    # """
    # echo "silent file to pdb"
    # rna_extract.linuxgccrelease \
    # -in:file:silent ${CASP_TARGET}.out \
    # -in:file:silent_struct_type rna \
    # -out:file:silent ${CASP_TARGET}.pdb
    # """
    # cmd = f"rna_extract.linuxgccrelease -in:file:silent {result_models_out_path} \
    #     -in:file:silent_struct_type rna \
    #     -out:file:silent {padded_pdb}"
    # subprocess.run(cmd, shell=True, env=env)

    
def append_residue_to_pose(added_seq, output_tmp_file, input_chunk_pdb):
    # ファイル末尾に
        # TER                                                                             
        # ##Begin comments##
        # BINARY SILENTFILE FULL_MODEL_PARAMETERS  FULL_SEQUENCE gcccggauagcucagucgguagagcagcgggcacuaugggcgcagugucaauggacgcugacgguacaggccagacaauuauugucugguauagugcccgcggguccaggguucaagucccuguucgggcgcca  CONVENTIONAL_RES_CHAIN A:1-134  INPUT_DOMAIN 1-30,98-134  WORKING 1-134
        # ##End comments##
        # ...
    # のようになっているところがあって、FULL_SEQUENCE の右の配列に文字を付け足せばOK

    with open(input_chunk_pdb, "r") as f:
        lines = f.readlines()
        with open(output_tmp_file, "w") as wf:
            print(f"Writing to {output_tmp_file}")
            for line in lines:
                # if line.startswith("##Begin comments##"):
                #     wf.write(line)
                #     wf.write(f"BINARY SILENTFILE FULL_MODEL_PARAMETERS  FULL_SEQUENCE {added_seq}\n")
                #     wf.write("##End comments##\n")
                # else:
                #     wf.write(line)
                if line.startswith("BINARY SILENTFILE FULL_MODEL_PARAMETERS"):
                    wf.write(f"BINARY SILENTFILE FULL_MODEL_PARAMETERS  FULL_SEQUENCE {added_seq}\n")
                else:
                    wf.write(line)
    return 
    

def main():
    # Initialize PyRosetta
    init(extra_options = f"-database {pyrosetta_database}")
    chunk_pdb = "/large/otgk/casp/casp16/utils/examples/R1205_FF2_S000001.pdb"
    pfree_secondary_structure_file = "/large/otgk/casp/casp16/2_R1205/ipknot.secstruct"
    tmp_pdb = "/large/otgk/casp/casp16/utils/examples/test_pfree_R1205_FF2_S000001.tmp.pdb"
    # with tmp pdb, we should added the residue to the head of the sequence


    # data loading
    pfree_chunk, pfree_ss = data_loading(chunk_pdb, pfree_secondary_structure_file)

    # adding one residue to the pose sequence
    adding = "a"
    added_sequence = adding + pfree_chunk.sequence()
    print(added_sequence)
    # added residue must not form any base pair in secondary structure
    padded_ss = "." + pfree_ss

    # change tmp_pose's seq to added_sequence
    append_residue_to_pose(added_sequence, tmp_pdb, chunk_pdb)
    # append_residue_to_pose(pfree_chunk.sequence(), tmp_pdb, chunk_pdb)


    print(added_sequence)

    # run rna_denovo
    nstruct = 3
    result_models_out_path = "/large/otgk/casp/casp16/utils/examples/testR1205_FF2_S000001.out"
    padded_pdb = "/large/otgk/casp/casp16/utils/examples/test_padded_R1205_FF2_S000001.pdb"
    print(f"Running rna_denovo with {nstruct} structures")
    run_rna_denovo(tmp_pdb, added_sequence, padded_ss, nstruct, rna_denovo_path, result_models_out_path, padded_pdb)




if __name__ == "__main__":
    main()


