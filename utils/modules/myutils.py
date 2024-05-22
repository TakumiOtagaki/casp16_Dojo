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

def read_singlefasta(fasta_file):
    with open(fasta_file, "r") as f:
        lines = f.readlines()
        sequence = lines[1].strip()
    return sequence

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