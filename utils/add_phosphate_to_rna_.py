import numpy as np
import argparse

def read_pdb(file_path):
    """Read PDB file and return atom coordinates and lines."""
    atoms = []
    lines = []
    chain_id = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                atoms.append({
                    'name': line[12:16].strip(),
                    'residue': line[17:20].strip(),
                    'chain': line[21],
                    'resid': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54])
                })
                if chain_id is None:
                    chain_id = line[21]  # Get the chain ID from the first ATOM line
            lines.append(line)
    return atoms, lines, chain_id


def determine_ribose_conformation(atoms, resid):
    """Determine the ribose conformation (C2'-endo or C3'-endo) of a given residue."""
    atom_coords = {atom['name']: np.array([atom['x'], atom['y'], atom['z']]) for atom in atoms if atom['resid'] == resid}
    c1 = atom_coords['C1\'']
    c2 = atom_coords['C2\'']
    c3 = atom_coords['C3\'']
    c4 = atom_coords['C4\'']
    o4 = atom_coords['O4\'']
    
    # Calculate torsion angle for C2'-endo or C3'-endo
    vec_c1c2 = c2 - c1
    vec_c3c2 = c2 - c3
    vec_c3c4 = c4 - c3
    vec_o4c4 = o4 - c4
    
    torsion_angle_c2_endo = np.dot(np.cross(vec_c1c2, vec_c3c2), np.cross(vec_c3c4, vec_o4c4))
    if torsion_angle_c2_endo > 0:
        return 'C2\'-endo'
    else:
        return 'C3\'-endo'

def add_phosphate(atoms, conformation, resid):
    """Add a phosphate group to the 5' end of the RNA, adjusting based on ribose conformation."""
    atom_coords = {atom['name']: np.array([atom['x'], atom['y'], atom['z']]) for atom in atoms if atom['resid'] == resid}
    if 'O5\'' in atom_coords:
        anchor_atom = atom_coords['O5\'']
    else:
        print("Warning: O5' atom not found, using C5' instead.")
        anchor_atom = atom_coords['C5\'']
    
    # Reference: Principles of Nucleic Acid Structure (S. Neidle)
    if conformation == 'C2\'-endo':
        phosphate_coords = {
            'P': anchor_atom + np.array([1.6, -0.2, 0.0]),
            'OP1': anchor_atom + np.array([2.1, 0.5, 0.8]),
            'OP2': anchor_atom + np.array([2.1, 0.5, -0.8]),
            'O5\'': anchor_atom
        }
    else: # C3'-endo
        phosphate_coords = {
            'P': anchor_atom + np.array([1.5, -0.2, 0.0]),
            'OP1': anchor_atom + np.array([2.0, 0.5, 0.8]),
            'OP2': anchor_atom + np.array([2.0, 0.5, -0.8]),
            'O5\'': anchor_atom
        }

    return phosphate_coords

def write_pdb(output_path, lines, phosphate_coords, resid, chain_id):
    """Write the new PDB file with the added phosphate group."""
    new_lines = []
    atom_index = 1
    phosphate_added = False
    temperature_factor = 0.0
    for line in lines:
        if line.startswith('ATOM') and int(line[22:26].strip()) == resid and not phosphate_added:
            # Add phosphate atoms before the first atom of the residue
            new_lines.append(f"ATOM  {atom_index:5d}  P    G {chain_id} {resid:4d}    {phosphate_coords['P'][0]:8.3f}{phosphate_coords['P'][1]:8.3f}{phosphate_coords['P'][2]:8.3f}  1.00 {temperature_factor}           P\n")
            atom_index += 1
            new_lines.append(f"ATOM  {atom_index:5d}  OP1  G {chain_id} {resid:4d}    {phosphate_coords['OP1'][0]:8.3f}{phosphate_coords['OP1'][1]:8.3f}{phosphate_coords['OP1'][2]:8.3f}  1.00 {temperature_factor}           O\n")
            atom_index += 1
            new_lines.append(f"ATOM  {atom_index:5d}  OP2  G {chain_id} {resid:4d}    {phosphate_coords['OP2'][0]:8.3f}{phosphate_coords['OP2'][1]:8.3f}{phosphate_coords['OP2'][2]:8.3f}  1.00 {temperature_factor}           O\n")
            atom_index += 1
            # new_lines.append(f"ATOM  {atom_index:5d}  O5'  G {chain_id} {resid:4d}    {phosphate_coords['O5\''][0]:8.3f}{phosphate_coords['O5\''][1]:8.3f}{phosphate_coords['O5\''][2]:8.3f}  1.00 {temperature_factor}           O\n")
            o5 = phosphate_coords["O5'"]
            new_lines.append(f"ATOM  {atom_index:5d}  O5'  G {chain_id} {resid:4d}    {o5[0]:8.3f}{o5[1]:8.3f}{o5[2]:8.3f}  1.00 {temperature_factor}           O\n")

            atom_index += 1
            phosphate_added = True

        if line.startswith('ATOM'):
            new_lines.append(f"ATOM  {atom_index:5d}{line[11:]}")
            atom_index += 1
        else:
            new_lines.append(line)
    
    with open(output_path, 'w') as file:
        for line in new_lines:
            file.write(line)

def check_residue_id(atoms):
    """Check if the residue ID for the 5' end is 1 and return the first residue ID."""
    residue_ids = [atom['resid'] for atom in atoms]
    first_residue_id = min(residue_ids)
    if first_residue_id != 1:
        print(f"Warning: The first residue ID is {first_residue_id}, not 1. Using {first_residue_id} as the residue ID for the 5' end.")
    return first_residue_id

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add a phosphate group to the 5\' end of an RNA PDB file.')
    parser.add_argument('--input_pdb', required=True, help='Path to the input PDB file')
    parser.add_argument('--output_pdb', required=True, help='Path to the output PDB file')
    parser.add_argument('--residue_id', type=int, default=None, help='Residue ID for the 5\' end (default: None)')
    args = parser.parse_args()

    input_pdb = args.input_pdb
    output_pdb = args.output_pdb

    atoms, lines, chain_ID = read_pdb(input_pdb)

    # If no residue ID is provided, use the first residue ID in the PDB file
    if args.residue_id is None:
        residue_id = check_residue_id(atoms)
    else:
        residue_id = args.residue_id

    conformation = determine_ribose_conformation(atoms, residue_id)
    phosphate_coords = add_phosphate(atoms, conformation, residue_id)
    write_pdb(output_pdb, lines, phosphate_coords, residue_id, chain_ID)

    print(f"Phosphate group added to the 5' end of residue {residue_id} with {conformation} conformation.")
    print(f"Output PDB file saved to {output_pdb}")