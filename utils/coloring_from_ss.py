import sys
def parse_layers(input_string):
    layers = {1: [], 2: [], 3: []}
    stack = []

    for index, char in enumerate(input_string, start=1):
        if char in "({[":
            layer_depth = 1 if char == '(' else (2 if char == '{' else 3)
            stack.append((char, index, layer_depth))
        elif char in ")}]":
            opening_char = '(' if char == ')' else ('{' if char == '}' else '[')
            # Find the matching opening character from the stack
            for i in range(len(stack)-1, -1, -1):
                if stack[i][0] == opening_char:
                    _, opening_index, layer_depth = stack.pop(i)
                    layers[layer_depth].append((opening_index, index))
                    break

    return layers


def generate_chimera_commands(layers, colors, model_id=0, chain_id='A'):
    commands = []
    for layer_depth, ranges in layers.items():
        color = colors.get(layer_depth)
        if color:
            for start, end in ranges:
                # Color only the start and end residues
                commands.append(f"color {color} #{model_id}:{start}.{chain_id};")
                commands.append(f"color {color} #{model_id}:{end}.{chain_id};")
    return commands

# input_string = "((((((((((.[[[[[[[.(((((...))))).{...((.((((((((((((((((......)))))).....(((.(((((..........))))))(((......)).)((.(((....(((((...(((....)))..))))).....)))....(((((........(...((((((((....))))))))......)..)))..))))).))))))))...(((..........)))...(((((((((((........))).))))))))(((((((((.......))))).))))........]]]]]]]((...(((((.((((.[[[[[[)))).....).))))...))..}....)))))))....)))))))))).]]]]]]"
input_string = input("Enter a secondary structure string: ")

layers = parse_layers(input_string)
print(layers)
# sys.exit()
colors = {
    1: 'red',    # 1st layer
    2: 'green',  # 2nd layer
    3: 'blue'    # 3rd layer
}


model_id = input("Enter the model ID: ")
chain_id = input("Enter the chain ID: ")
chimera_commands = generate_chimera_commands(layers, colors, model_id, chain_id)

for command in chimera_commands:
    print(command)

"""
AF only
((((((((((.[[[[[[[.(((((...))))).{...((.((((((((((((((((......)))))).....(((.(((((((..........))))))(((......)).)((.(((....(((((...(((....)))..))))).....))).....(((((........(...((((((((....))))))))......)..)))..))))).))))))))...(((..........)))...(((((((((((........))).))))))))(((((((((.......))))).))))........]]]]]]].((...(((((.((((.[[[[[[)))).....).))))...))..}....)))))))....)))))))))).]]]]]]

2A64
..(((((((((((((((((((.(((((((((....))))))))).{...[[.[[[[[(((((((((.((......)))))).....(((((((((((..........))))))((((((------------))))).)((.(((-------------------------------------))).....((.((-----------------------------------------))...))))))))))))))...(((...------..)))...(((((((...(........)...)))))))((((((..........))))))........))))))).((...(.((--------------------)).)...))..}....]]]]]]]..).)))))))))))..
"""