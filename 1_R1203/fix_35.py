from Bio import PDB
import os
import subprocess

"""
extract_and_save_structure: 指定されたPDBファイルから特定の範囲の残基を抽出して新しいPDBファイルに保存します。
merge_pdbs: 抽出された複数のPDBファイルをマージして一つのファイルにします。この関数は、すべての原子レコードを含めるようにファイルを結合します。
adjust_atom_indices: マージしたPDBファイルの原子インデックスを調整します。これは必要に応じてコメントアウトまたは調整してください。
run_rosetta: Rosettaのrna_denovoコマンドを使用してエネルギー最適化と構造予測を行います。このコマンドはRosettaがインストールされ、適切に設定されている環境で実行する必要があります。
"""

farfar2_path = "/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease"

def cif2pdb(cif_file, pdb_file):
    parser = PDB.MMCIFParser()
    structure = parser.get_structure('RNA', cif_file)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)


def extract_and_save_structure(input_pdb, output_pdb, chain_id, ranges):
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    structure = parser.get_structure('RNA', input_pdb)
    
    class ResidueSelect(PDB.Select):
        def accept_residue(self, residue):
            return any(start <= residue.id[1] <= end for start, end in ranges if residue.parent.id == chain_id)

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
    io.save(pdb_file)

def run_rosetta(pdb_file):
    cmd = f'{farfar2_path} -s {pdb_file} \
         -sequence "gcccggauagcucagucgguagagcagcgggcacuaugggcgcagugucaauggacgcugacgguacaggccagacaauuauugucugguauagugcccgcggguccaggguucaagucccuguucgggcgcca"\
         -secstruct "(((((((..((((........)))).(((((((((((.(.(((((((((....))))))).)).)....(((((((((...)))))))))))))))))))).....(((((.......))))))))))))...." \
         -nstruct 10 \
        -out:file:silent /large/otgk/casp/casp16/1_R1203/for_merging/result_models.out \
        -minimize_rna true'
    # os.system("source ~/.bashrc")
    # os.system("module load rosetta")
    # os.system(cmd)
    env = os.environ.copy()
    env["ROSETTA3"] = "/large/otgk/app/rosetta/v2024.15/source"
    subprocess.run(cmd, shell=True)

# 抽出する範囲を定義
alphafold3_cif = "/large/otgk/casp/casp16/1_R1203/for_merging/fold_r1203_s1_model_0.cif"
alphafold3_pdb = "/large/otgk/casp/casp16/1_R1203/for_merging/alphafold3.pdb"
cif2pdb(alphafold3_cif, '/large/otgk/casp/casp16/1_R1203/for_merging/alphafold3.pdb')
print("Extracting and saving structures...")
extract_and_save_structure(alphafold3_pdb, '/large/otgk/casp/casp16/1_R1203/for_merging/alphafold3_terminals.pdb', 'A', [(1, 30), (98, 134)])
# 今回はマージしなくていい。farfar2 の中央部分は固定しないことにする。-s に渡すのは alphafold3 だけ



print("Running Rosetta...")
run_rosetta('alphafold3_terminals.pdb')