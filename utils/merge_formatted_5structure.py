import os
import argparse

def merge_pdb_files(input_dir, output_file, casp_target_id):
    pdb_files = [f for f in os.listdir(input_dir) if f.endswith('.pdb')]
    if len(pdb_files) != 5:
        raise ValueError(f"Input directory must contain exactly 5 PDB files: {input_dir}. Found {len(pdb_files)} files.")
    pdb_files.sort()  # 必要に応じてファイル名でソート

    header = f"PFRMAT TS\nTARGET {casp_target_id}\nAUTHOR 0652-1349-3580\nREMARK\nMETHOD newMXfold2 + FARFAR2\n"

    with open(output_file, 'w') as outfile:
        # ヘッダー情報を追加
        outfile.write(header)
        model_index = 1
        
        for pdb_file in pdb_files:
            with open(os.path.join(input_dir, pdb_file), 'r') as infile:
                outfile.write(f"MODEL     {model_index}\n")
                outfile.write("PARENT N/A\n")
                
                for line in infile:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        outfile.write(line)
                
                outfile.write("TER\n")
                outfile.write("ENDMDL\n")
                model_index += 1
            outfile.write("END\n")
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge multiple PDB files into one.')
    parser.add_argument('--input_dir', '-i', required=True, help='Directory containing input PDB files')
    parser.add_argument('--output_file', "-o", required=True, help='Output PDB file')
    parser.add_argument("--casp_target_id", "-t", help="CASP target ID, which will be written in the header")
    args = parser.parse_args()

    merge_pdb_files(args.input_dir, args.output_file, args.casp_target_id)