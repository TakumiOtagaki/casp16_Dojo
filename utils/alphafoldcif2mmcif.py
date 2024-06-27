from Bio import PDB
import sys

# CIFファイルのパスと出力mmCIFファイルのパスを設定
input_cif = 'input.cif'
output_mmcif = 'output.mmcif'

# CIFパーサーを使用して構造を読み込む
parser = PDB.MMCIFParser()
structure = parser.get_structure('structure', input_cif)

# mmCIFフォーマットで出力
io = PDB.MMCIFIO()
io.set_structure(structure)
io.save(output_mmcif)

print(f"Converted {input_cif} to {output_mmcif}")