from Bio import PDB
import os

def cif2pdb(cif_file, pdb_file):
    parser = PDB.MMCIFParser()
    structure = parser.get_structure('RNA', cif_file)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def extract_and_save_structure(input_pdb, output_pdb, chain_id, start_resi, end_resi):
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    structure = parser.get_structure('RNA', input_pdb)
    class ResidueSelect(PDB.Select):
        def accept_residue(self, residue):
            if residue.id[1] >= start_resi and residue.id[1] <= end_resi and residue.parent.id == chain_id:
                return 1
            return 0
    io.set_structure(structure)
    io.save(output_pdb, select=ResidueSelect())

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

def run_rosetta(pdb_file):
    cmd = f"rna_denovo -s {pdb_file} -nstruct 10 -out:file:silent result_models.out"
    os.system(cmd)

# 抽出する範囲を定義
extract_and_save_structure('alphafold3.pdb', 'alphafold3_extracted.pdb', 'A', 1, 30)
extract_and_save_structure('alphafold3.pdb', 'alphafold3_extracted2.pdb', 'A', 98, 134)
extract_and_save_structure('farfar2.pdb', 'farfar2_extracted.pdb', 'A', 31, 97)

# ファイルをマージしてエネルギー最適化
merge_pdbs('merged_output.pdb', 'alphafold3_extracted.pdb', 'farfar2_extracted.pdb', 'alphafold3_extracted2.pdb')
run_rosetta('merged_output.pdb')