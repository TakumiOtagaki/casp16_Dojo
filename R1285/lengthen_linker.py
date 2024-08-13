import re

def lengthen_linker(sequence, structure, linker_length=100):
    # Find all NNN sequences and their positions
    print("sequence", sequence, len(sequence))
    matches = [(m.start(), m.end()) for m in re.finditer(r'NNN', sequence)]
    print("matches", matches)
    
    new_sequence = sequence
    new_structure = structure
    
    # Offset to account for the increase in length
    offset = 0
    
    for start, end in matches:
        # Calculate new positions considering the offset
        start += offset
        end += offset
        
        # Replace NNN with N*linker_length in the sequence
        new_sequence = new_sequence[:start] + 'C' * linker_length + new_sequence[end:]
        
        # Replace corresponding structure with '.'*linker_length
        new_structure = new_structure[:start] + '.' * linker_length + new_structure[end:]
        
        # Update offset
        offset += linker_length - (end - start)
    
    return new_sequence, new_structure

def read_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
    return lines

def write_file(filepath, lines):
    with open(filepath, 'w') as file:
        file.writelines(lines)

def process_file(input_filepath, output_filepath, linker_length=100):
    lines = read_file(input_filepath)
    
    seq_name = lines[0].strip()
    sequence = lines[1].strip()
    structure = lines[2].strip()
    
    print("seq_name", seq_name, len(seq_name))
    print("sequence", sequence, len(sequence))
    print("structure", structure, len(structure))
    
    new_seq, new_struct = lengthen_linker(sequence, structure, linker_length)
    
    output_lines = [seq_name + "\n", new_seq + "\n", new_struct + "\n"]
    write_file(output_filepath, output_lines)


if __name__ == "__main__":
    input_filepath = '/large/otgk/casp/casp16/R1285/R1285.ShortLinker.secstruct'
    output_filepath = '/large/otgk/casp/casp16/R1285/R1285.LongLinker.secstruct'
    process_file(input_filepath, output_filepath)