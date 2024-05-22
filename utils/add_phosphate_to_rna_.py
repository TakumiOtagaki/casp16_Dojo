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

from Bio.PDB import PDBParser, PDBIO, Select
import os
import subprocess

from modules.myutils import read_ss_file


rna_denovo_path = "/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease"
pyrosetta_database = "/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/pyrosetta/database/"

# Read a secondary structure file
# secstruct format:
    # >filename
    # sequence
    # secondary structure
# def read_ss_file(ss_file):
    # in modules/myutils.py


def data_loading(InitialStructurePDB, initialSecondaryStructure_path):
    initialPose = pose_from_pdb(InitialStructurePDB)
    _, seq, initialSecondaryStructure = read_ss_file(initialSecondaryStructure_path)
    assert initialPose.sequence() == seq
    print("Data loading is done.")
    return initialPose, initialSecondaryStructure

def run_rna_denovo(modifiedStructurePDB, extendedSequence, modifiedss, nstruct,  rna_denovo_path, finalStructureOut, finalStructurePDB):
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
    env["ROSETTA3"] = "/large/otgk/app/rosetta/v2024.15/source"
    subprocess.run(cmd, shell=True, env=env)

    # .out to pdb

    cmd = f"rna_extract.linuxgccrelease -in:file:silent {finalStructureOut} \
        -in:file:silent_struct_type rna \
        -out:file:silent {finalStructurePDB}"
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
    # Initialize PyRosetta
    init(extra_options = f"-database {pyrosetta_database}")
    initialStructurePDB = "/large/otgk/casp/casp16/utils/examples/R1205_FF2_S000001.pdb"
    secoudaryStructureFile = "/large/otgk/casp/casp16/2_R1205/ipknot.secstruct"
    modifiedStructurePDB = "/large/otgk/casp/casp16/utils/examples/modifiedR1205_FF2_S000001.pdb"
    finalStructureOut = "/large/otgk/casp/casp16/utils/examples/padded_R1205_FF2_S000001.out"
    finalStructurePDB = "/large/otgk/casp/casp16/utils/examples/padded_R1205_FF2_S000001.pdb"
    # with tmp pdb, we should added the residue to the head of the sequence


    # data loading
    initialPose, initialSecondaryStructure = data_loading(initialStructurePDB, secoudaryStructureFile)

    # adding one residue to the pose sequence
    adding = "a" # adenine consists of 31 atoms in pdb format
    extendedSequence = adding + initialPose.sequence()
    print(extendedSequence)
    # added residue must not form any base pair in secondary structure
    modifiedss = "." + initialSecondaryStructure


    # append residue to the modified structure
    # append_residue_to_pose(extendedSequence, modifiedStructurePDB, initialStructurePDB, num_of_atoms)
    renumber_residues(initialStructurePDB, modifiedStructurePDB, len(adding))


    print(extendedSequence)

    # run rna_denovo
    nstruct = 3
    print(f"Running rna_denovo with {nstruct} structures")
    run_rna_denovo(modifiedStructurePDB, extendedSequence, modifiedss, nstruct, rna_denovo_path, finalStructureOut, finalStructurePDB)




if __name__ == "__main__":
    main()



