from Bio import PDB
import argparse
import warnings

def parse_args():
    description = """
    Process and format a PDB file.
    The first residue is removed by default but can be kept with the --keep_first_residue option.
    The output file is formatted with renumbered residues and atoms, chain IDs, and occupancy set to 1.00.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input_file", "-i",  help="Path to the input PDB file.")
    parser.add_argument("--output_file", "-o", help="Path to the output formatted PDB file.")

    parser.add_argument("--keep_first_residue", action="store_false", 
                        help="Keep the first residue. By default, the first residue is removed.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print additional information.")
    args = parser.parse_args()

    if not args.input_file:
        parser.error("Please provide an input file with the -i option.")
    if not args.output_file:
        parser.error("Please provide an output file with the -o option.")
    
    return args

def verbose_print(verbose):
    def printer(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
    return printer

class NonHydrogenSelector(PDB.Select):
    """ 水素以外の原子だけを選択するセレクタ """
    def accept_atom(self, atom):
        return atom.element != 'H'
    
def pdb_formater(input_file, output_file, remove_first_residue=True, verbose=False):
    vprint = verbose_print(verbose)

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB', input_file)
    
    io = PDB.PDBIO()
    
    # Remove the first residue and renumber residues and atoms
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            if remove_first_residue and residues:
                first_residue = residues[0]
                keep_atoms = {'P', 'OP1', 'OP2', 'O5\''}
                for atom in list(first_residue.get_atoms()):
                    if atom.name not in keep_atoms:
                        first_residue.detach_child(atom.id)
                if not first_residue:
                    chain.detach_child(first_residue.id)

    # Renumber residues and atoms, change chain IDs, remove hydrogens, and set occupancy
    atom_number = 1
    residue_number = 1
    chain_counter = 0

    for model in structure:
        # チェーンIDの変更: 現在のチェーンIDが chain_id_map に存在しない場合、新しいチェーンIDを作成します。元のチェーンIDがアルファベットである場合、新しいチェーンIDを 'A', 'B', 'C' などのアルファベットに変更します。そうでない場合、新しいチェーンIDを '0', '1', '2' などの数字に変更します。
        # チェーンマップの作成: chain_id_map を使用して元のチェーンIDと新しいチェーンIDのマッピングを保持し、重複を避けます。
        # チェーンIDの更新: チェーンのIDを新しいIDに変更します。
        # RNA: chain_id = 0, 1, ... (0 origin)
        for chain in model:
            print(f"Processing chain {chain.id} in model {model.id}")
            new_chain_id = str(chain_counter)
            chain.id = new_chain_id
            chain_counter += 1
            

            # 残基IDの更新: 各残基に新しい残基番号を割り当てます。残基番号は1から始まり、順次増加します。
            # 原子の処理: 各残基の各原子について処理を行います。
            # 水素原子の削除: 原子が水素（'H'）である場合、その原子を残基から削除します。
            # 原子番号の更新: 水素以外の原子について、原子番号を1から始めて順次増加させます。
            # Occupancyの設定: 各原子のOccupancyを1.00に設定します。
            # 原子名の長さの警告: 原子名が2文字以上である場合、警告を表示します。
                # 三列目のものではなくて、12 列目のもの
            for residue in chain.get_residues():
                residue.id = (residue.id[0], residue_number, residue.id[2])
                residue_number += 1
                for atom in residue.get_atoms():
                    # if atom.element == 'H':
                    #     vprint(f"Atom {atom.element} in residue {residue.resname} is a hydrogen atom: id={atom.id}, serial_number={atom.serial_number}")
                    #     # residue.detach_child(atom.id)
                    if len(atom.element) > 1:
                        print(f"atom.element: {atom.element}")
                        warnings.warn(f"Atom {atom.name.strip()} in residue {residue.resname} has a long atom name: id={atom.id}, serial_number={atom.serial_number}")
                    else:
                        atom.serial_number = atom_number
                        atom_number += 1
                        atom.occupancy = 1.00

    io.set_structure(structure)
    io.save(output_file, NonHydrogenSelector())

def main():
    args = parse_args()

    pdb_formater(args.input_file, args.output_file, remove_first_residue=args.keep_first_residue, verbose=args.verbose)
    print(f"Formatted PDB file saved to {args.output_file}")

if __name__ == "__main__":
    main()