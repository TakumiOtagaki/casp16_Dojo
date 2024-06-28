import argparse
from Bio import PDB

def split_chain(pdb_file, output_file):
    parser = PDB.PDBParser(QUIET=True)
    io = PDB.PDBIO()
    
    # Parse the structure
    structure = parser.get_structure('structure', pdb_file)
    
    # Check for multiple models
    if len(structure) > 1:
        raise ValueError("Multiple models found in the PDB file.")
    
    model = structure[0]
    
    # Check for multiple chains
    if len(model) > 1:
        raise ValueError("Multiple chains found in the PDB file.")
    
    chain = next(model.get_chains())
    
    # Get all residues
    residues = list(chain.get_residues())
    num_residues = len(residues)
    
    if num_residues < 2:
        raise ValueError("Not enough residues to split the chain.")
    
    # Determine the split point
    split_point = num_residues // 2
    print(f"Splitting the chain at residue {split_point}")
    
    # Create new chains
    chain_a = PDB.Chain.Chain('A')
    chain_b = PDB.Chain.Chain('B')
    
    # Assign residues to new chains
    for i, residue in enumerate(residues):
        if i == split_point:
            print(f"Splitting residue {residue}")
            print(f"Chain A: {chain_a}")
        if i < split_point:
            chain_a.add(residue)
        else:
            new_residue = residue.copy()
            new_residue.id = (' ', i - split_point + 1, ' ')
            chain_b.add(new_residue)
    
    # Create a new structure to hold the split chains
    new_structure = PDB.Structure.Structure('new_structure')
    new_model = PDB.Model.Model(0)
    new_model.add(chain_a)
    new_model.add(chain_b)
    new_structure.add(new_model)
    
    # Save the new structure to a file
    io.set_structure(new_structure)
    io.save(output_file)
    print(f"New PDB file saved as {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split an RNA chain in a PDB file into two chains.')
    parser.add_argument('--input_pdb', "-i", type=str, help='Input PDB file')
    parser.add_argument('--output_pdb', "-o", type=str, help='Output PDB file')
    
    args = parser.parse_args()
    
    try:
        split_chain(args.input_pdb, args.output_pdb)
    except Exception as e:
        print(f"Error: {e}")
        exit(1)