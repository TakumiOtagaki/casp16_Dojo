import argparse
import sys
# cat /Users/ootagakitakumi/R1221s2/岩切先生/R1221s2FixPart.pdb  | awk -F "[[:space:]]+" '{print $6}' | uniq > af3_fixpart_resid_list.txt
# for af3.R1221s2_fixed.pdb,

# python3 '/large/otgk/casp/casp16/get_bpseq_from_ss_and_resid_list.py' -ss "((((((((((.[[[[[[[.(((((...))))).{...((.((((((((((((((((......)))))).....(((.(((((((..........))))))(((......)).)((.(((....(((((...(((....)))..))))).....))).....(((((........(...((((((((....))))))))......)..)))..))))).))))))))...(((..........)))...(((((((((((........))).))))))))(((((((((.......))))).))))........]]]]]]].((...(((((.((((.[[[[[[)))).....).))))...))..}....)))))))....)))))))))).]]]]]]" \
# -r /large/otgk/casp/casp16/11_R1221s2/af3_fixpart_resid_list.txt \
# -seq "GUGAUAUUUCGGGUAAUCGCUAUAUUAUAUAGAGGAAAGUCCAUGCUCACACAGUCUGAGAUGAUUGUAGUGUUCGUGCUUGAUGAAACAAUAAAUCAAGGCAUUAAUUUGACGGCAAUGAAAUAUCCUAAGUCUUUCGAUAUGGAUAGAGUAAUUUGAAAGUGCCACAGUGACGUAGCUUUUAUAGAAAUAUAAAAGGUGGAACGCGGUAAACCCCUCGAGUGAGCAAUCCAAAUUUGGUAGGAGCACUUGUUUAACGGAAUUCAACGUAUAAACGAGACACACUUCGCGAAAUGAAGUGGUGUAGACAGAUGGUUAUCACCUGAGUACCAGUGUGACUAGUGCACGUGAUGAGUACGAUGGAACAGAACAUGGCUUAUAGAAAUAUCACUACUAGU" \
# -o /large/otgk/casp/casp16/11_R1221s2/af3_fixpart.bpseq





def parse_args():
    parser = argparse.ArgumentParser(description="Generate BPSEQ file from a secondary structure string and a list of residue IDs.")
    parser.add_argument("--sec_structure", "-ss", type=str, help="Secondary structure string")
    parser.add_argument("--resid_list", "-r", type=str, help="Comma-separated list of residue IDs file. ")
    parser.add_argument("--sequence", "-seq", type=str, help="Sequence string")
    parser.add_argument("--output", "-o", type=str, help="Output BPSEQ file")
    return parser.parse_args()

def parse_layers(ss):
    layers = {1: [], 2: [], 3: []}
    stack = {1: [], 2: [], 3: []}
    
    for index, char in enumerate(ss, start=1):
        if char == ".": continue
        if char in "([{":
            layer_depth = 1 if char == '(' else (2 if char == '{' else 3)
            stack[layer_depth].append((char, index))
        elif char in ")]}":
            layer_depth = 1 if char == ')' else (2 if char == '}' else 3)
            if stack[layer_depth]:
                opening_char, opening_index = stack[layer_depth].pop(-1)
                print(opening_index, index)
                layers[layer_depth].append((opening_index, index))
    
    return layers

def generate_bpseq(seq, sec_structure, res_id_list):
    layers = parse_layers(sec_structure)
    print(layers)
    pairs = {}

    for depth, pair_list in layers.items():
        for i, j in pair_list:
            pairs[i] = j
            pairs[j] = i

    bpseq_data = []
    warnings = []
    for i, base in enumerate(seq, start=1):
        if i not in res_id_list:
            pair_symbol = '.'
        else: # i in res_id_list 
            if i in pairs:
                pair = pairs[i]
                if pair in res_id_list:
                    pair_symbol = str(pair)
                else:
                    warnings.append(f"Warning: Base {i} is paired with {pair}, which is not in res_id_list")
                    pair_symbol = '.'
            else:
                pair_symbol = 'x'
        # if i in pairs: # i が塩基対を形成している場合
        #     pair = pairs[i]
        #     if i in res_id_list and pair in res_id_list: # i が fix 対象だし、pair も fix 対象の場合
        #         pair_symbol = str(pair)
        #     elif i in res_id_list or pair in res_id_list: # i が fix 対象だが、pair が fix 対象でない場合
        #         warnings.append(f"Warning: Base {i} is paired with {pair}, which is not in res_id_list")
        #         pair_symbol = '.' # warning 
        #     else: # i が fix 対象でないし、pair も fix 対象でない場合
        #         print(i, "stopping here")
        #         pair_symbol = '.' # 制約なし
        #         sys.exit()
        # elif i in res_id_list:
        #     pair_symbol = 'x' # i は塩基対をくまないし、i が fix 対象のとき
        # else: # i 塩基対をくまないし、i が fix 対象でないとき
        #     pair_symbol = '.' # 制約なし

        bpseq_data.append((i, base, pair_symbol))

    return bpseq_data, warnings

def write_bpseq(bpseq_data, filename):
    with open(filename, 'w') as f:
        for i, base, pair in bpseq_data:
            f.write(f"{i} {base} {pair}\n")


args = parse_args()

# 入力データ
# sequence = "GCGCUUCGAAAGCGAAGCUCGCUCUUGUUUCGAAUCGGAAGCGACUCGCGGCUCUGCGCGAUUUUUUGAAUGCCGUCUUCGUCGGUGAAUGCGCUCGGCUCGCGGUGAAAGCUCGAACGUGUUUGUCGACUCGUUUGCUUCGACCGUUCGCCCGUUCGCGGAUGCUAGCCGUCCGAAGGCUCGCGGGAUUUUCCUUAGGAUCCGUCUUCGCUCGAUCUCUCGCUCUCGUUUGUGCUCGCUCGCGAUUU"
# sec_structure = "((((((((((.[[[[[[[.(((((...))))).{...((.((((((((((((((((......)))))).....(((.(((((..........))))))(((......)).)((.(((....(((((...(((....)))..))))).....))).....(((((........(...((((((((....))))))))......)..)))..))))).))))))))...(((..........)))...(((((((((((........))).))))))))(((((((((.......))))).))))........]]]]]]].((...(((((.((((.[[[[[[)))).....).))))...))..}....)))))))....)))))))))).]]]]]]"
original_sec_structure = args.sec_structure
sequence = args.sequence
res_id_list = []
with open(args.resid_list, 'r') as f:
    for line in f:
        pos = line.strip()
        if pos == "":
            continue
        # print(pos)
        res_id_list.append(int(pos))
print(res_id_list)


# BPSEQファイルの生成
bpseq_data, warnings = generate_bpseq(sequence, original_sec_structure, res_id_list)
write_bpseq(bpseq_data, args.output)