import os
from Bio import PDB

def split_pdb_by_chain(input_pdb, output_dir):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb)

    io = PDB.PDBIO()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for model in structure:
        for chain in model:
            chain_id = chain.id
            chain_structure = structure.copy()
            for other_chain in list(chain_structure.get_chains()):
                if other_chain.id != chain_id:
                    chain_structure[0].detach_child(other_chain.id)
            output_file = os.path.join(output_dir, f"chain_{chain_id}.pdb")
            io.set_structure(chain_structure)
            io.save(output_file)
            print(f"Saved chain {chain_id} to {output_file}")

if __name__ == "__main__":
    input_pdb = "/large/otgk/casp/casp16/R0254/multimer/n3_2/raw.pdb"  # 入力PDBファイルのパスを指定
    output_dir = "/large/otgk/casp/casp16/R0254/multimer/n3_2/unformatted"  # 出力ディレクトリのパスを指定
    split_pdb_by_chain(input_pdb, output_dir)