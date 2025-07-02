import re

def fix_pdb_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as file:
        for line in lines:
            if line.startswith("ATOM"):
                # 正規表現でA+4桁の数値の前にスペースを追加
                line = re.sub(r'A(\d{4})', r'A \1', line)
                # 正規表現で -\d+-\d+ とマイナスの数値がスペースで区切られずに並んでいるところを整える
                line = re.sub(r'(\d+)(-\d+)', r'\1 \2', line)
            file.write(line)

# 使用例
fix_pdb_file('/large/otgk/casp/casp16/R1285/farfar2_result/pdb/nomini_nstruct1_lesscycle/S_000001.pdb',
             '/large/otgk/casp/casp16/R1285/farfar2_result/pdb/nomini_nstruct1_lesscycle/S_000001_fixed.pdb')