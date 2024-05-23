from Bio import PDB
import argparse
import warnings

from Bio.PDB import PDBParser

def validate_pdb_file(pdb_file_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('PDB', pdb_file_path)
    cnt_of_warnings = 0

    atom_index = 1
    residue_index = 1
    prev_chain_id = None
    valid = True

    expected_atoms = ['P', 'OP2', 'OP1', 'O5\'']
    # 先頭のモデルとチェーンの先頭の4つの原子を取得
    for model in structure:
        for chain in model:
            # 二列目（index）1: P, 2: OP2, 3: OP1, 4: O5'
            first_residue = next(chain.get_residues(), None)
            if first_residue:
                atoms = [atom.name for atom in first_residue.get_atoms()][:4]
                if atoms != expected_atoms:
                    warnings.warn(f"First four atoms in chain {chain.id} are not {expected_atoms}. Found {atoms}.")
                    cnt_of_warnings += 1
                    valid = False


                        


    for model in structure:
        # if chain ID starts from non 0; warning
        if [int(chain.id) for chain in model if chain.id.isdigit()] and min([int(chain.id) for chain in model if chain.id.isdigit()]) != 0:
            warnings.warn("RNA Chain ID should be integer and start from 0.")
            cnt_of_warnings += 1


        for chain in model:
            # Chain ID check for RNA chains starting from 0
            current_chain_id = int(chain.id) if chain.id.isdigit() else -1
            if prev_chain_id is None:
                prev_chain_id = current_chain_id
            else:
                if current_chain_id != prev_chain_id + 1:
                    warnings.warn(f"Chain ID should increment by 1. Found {chain.id} after {prev_chain_id}.")
                    cnt_of_warnings += 1
                    valid = False
            prev_chain_id = current_chain_id

            for residue in chain:
                # Check if residue number is sequentially increasing
                if residue.id[1] != residue_index:
                    warnings.warn(f"Residue numbers not sequential at {residue.id[1]}. Expected {residue_index}.")
                    cnt_of_warnings += 1
                    valid = False
                residue_index += 1

                for atom in residue:
                    # Check if atom index is sequentially increasing
                    if atom.serial_number != atom_index:
                        warnings.warn(f"Atom numbers not sequential at {atom.serial_number}. Expected {atom_index}.")
                        cnt_of_warnings += 1
                        valid = False
                    atom_index += 1

                    # Check for hydrogen atoms
                    if atom.element == 'H':
                        warnings.warn(f"Hydrogen atom found: {atom.name} in residue {residue.resname}.")
                        cnt_of_warnings += 1
                        valid = False

                    # Check for occupancy value
                    if atom.occupancy != 1.00:
                        warnings.warn(f"Occupancy not set to 1.00 for atom {atom.name} in residue {residue.resname}.")
                        cnt_of_warnings += 1
                        valid = False

                    # Check for long atom names
                    if len(atom.element) > 1:
                        warnings.warn(f"Atom name longer than one character: {atom.name} in residue {residue.resname}.")
                        cnt_of_warnings += 1
                        valid = False

                    # Check if B-factor is 0
                    if atom.bfactor != 0:
                        warnings.warn(f"B-factor is not zero for atom {atom.name} in residue {residue.resname}. Found {atom.bfactor}.")
                        cnt_of_warnings += 1
                        valid = False

    return valid, cnt_of_warnings

def main():
    parser = argparse.ArgumentParser(description="Validate a PDB file against specific structural criteria.")
    parser.add_argument("--pdb_file", "-i", help="Path to the PDB file to validate.")
    args = parser.parse_args()

    valid, cnt = validate_pdb_file(args.pdb_file)
    if valid:
        print("PDB file is valid.")
    else:
        print("PDB file has errors; plz check the warnings above.")
        print(f"Total warnings: {cnt}")

    
if __name__ == "__main__":
    main()