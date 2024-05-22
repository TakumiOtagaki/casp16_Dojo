# conda activate rosetta
# conda config --add channels https://conda.rosettacommons.org
# conda install pyrosetta


# pyrosetta database path is needed;
    # $(pip show pyrosetta)/pyrosetta/database

from pyrosetta import *
from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.constraint_generator import *
from pyrosetta.rosetta.protocols.relax import ClassicRelax

from Bio import PDB
import os
import subprocess

rna_denovo_path = "/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease"
pyrosetta_database = "/home/otgk/.conda/envs/rosetta/lib/python3.10/site-packages/pyrosetta/database/"

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


def phosphate_free_structure(pdb_file):
    pose = pose_from_pdb(pdb_file)
    print(pose.sequence())
    return pose

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

# Read a secondary structure file
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

def data_loading(pfree_pdb_path, pfree_ss_path):
    pfree_chunk = phosphate_free_structure(pfree_pdb_path)
    _, seq, pfree_ss = read_ss_file(pfree_ss_path)
    assert pfree_chunk.sequence() == seq
    print("Data loading is done.")
    return pfree_chunk, pfree_ss

def run_rna_denovo(pfree_pdb, padded_seq, padded_ss, nstruct,  rna_denovo_path, result_models_out_path, padded_pdb):
    if len(padded_seq) != len(padded_ss):
        print("Error: length of sequence and secondary structure does not match.")
        return
    
    cmd = f'{rna_denovo_path} \
         -s {pfree_pdb} \
         -sequence "{padded_seq}" \
         -secstruct "{padded_ss}" \
         -nstruct {nstruct} \
        -out:file:silent {result_models_out_path} \
        -minimize_rna true'
    print(f"\n\ncmd: {cmd}\n\n\n")
    # os.system("source ~/.bashrc")
    # os.system("module load rosetta")
    # os.system(cmd)
    env = os.environ.copy()
    env["ROSETTA3"] = "/large/otgk/app/rosetta/v2024.15/source"
    subprocess.run(cmd, shell=True, env=env)

    # .out to pdb
    # """
    # echo "silent file to pdb"
    # rna_extract.linuxgccrelease \
    # -in:file:silent ${CASP_TARGET}.out \
    # -in:file:silent_struct_type rna \
    # -out:file:silent ${CASP_TARGET}.pdb
    # """
    # cmd = f"rna_extract.linuxgccrelease -in:file:silent {result_models_out_path} \
    #     -in:file:silent_struct_type rna \
    #     -out:file:silent {padded_pdb}"
    # subprocess.run(cmd, shell=True, env=env)

    
def append_residue_to_pose(added_seq, output_tmp_file, input_chunk_pdb):
    # ファイル末尾に
        # TER                                                                             
        # ##Begin comments##
        # BINARY SILENTFILE FULL_MODEL_PARAMETERS  FULL_SEQUENCE gcccggauagcucagucgguagagcagcgggcacuaugggcgcagugucaauggacgcugacgguacaggccagacaauuauugucugguauagugcccgcggguccaggguucaagucccuguucgggcgcca  CONVENTIONAL_RES_CHAIN A:1-134  INPUT_DOMAIN 1-30,98-134  WORKING 1-134
        # ##End comments##
        # ...
    # のようになっているところがあって、FULL_SEQUENCE の右の配列に文字を付け足せばOK

    with open(input_chunk_pdb, "r") as f:
        lines = f.readlines()
        with open(output_tmp_file, "w") as wf:
            print(f"Writing to {output_tmp_file}")
            for line in lines:
                # if line.startswith("##Begin comments##"):
                #     wf.write(line)
                #     wf.write(f"BINARY SILENTFILE FULL_MODEL_PARAMETERS  FULL_SEQUENCE {added_seq}\n")
                #     wf.write("##End comments##\n")
                # else:
                #     wf.write(line)
                if line.startswith("BINARY SILENTFILE FULL_MODEL_PARAMETERS"):
                    wf.write(f"BINARY SILENTFILE FULL_MODEL_PARAMETERS  FULL_SEQUENCE {added_seq}\n")
                else:
                    wf.write(line)
    return 
    

def main():
    # Initialize PyRosetta
    init(extra_options = f"-database {pyrosetta_database}")
    chunk_pdb = "/large/otgk/casp/casp16/utils/examples/R1205_FF2_S000001.pdb"
    pfree_secondary_structure_file = "/large/otgk/casp/casp16/2_R1205/ipknot.secstruct"
    tmp_pdb = "/large/otgk/casp/casp16/utils/examples/test_pfree_R1205_FF2_S000001.tmp.pdb"
    # with tmp pdb, we should added the residue to the head of the sequence


    # data loading
    pfree_chunk, pfree_ss = data_loading(chunk_pdb, pfree_secondary_structure_file)

    # adding one residue to the pose sequence
    adding = "a"
    added_sequence = adding + pfree_chunk.sequence()
    print(added_sequence)
    # added residue must not form any base pair in secondary structure
    padded_ss = "." + pfree_ss

    # change tmp_pose's seq to added_sequence
    append_residue_to_pose(added_sequence, tmp_pdb, chunk_pdb)
    # append_residue_to_pose(pfree_chunk.sequence(), tmp_pdb, chunk_pdb)


    print(added_sequence)

    # run rna_denovo
    nstruct = 3
    result_models_out_path = "/large/otgk/casp/casp16/utils/examples/testR1205_FF2_S000001.out"
    padded_pdb = "/large/otgk/casp/casp16/utils/examples/test_padded_R1205_FF2_S000001.pdb"
    print(f"Running rna_denovo with {nstruct} structures")
    run_rna_denovo(tmp_pdb, added_sequence, padded_ss, nstruct, rna_denovo_path, result_models_out_path, padded_pdb)




if __name__ == "__main__":
    main()




"""
ERRORS:

cmd: /large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease          -s /large/otgk/casp/casp16/utils/examples/test_pfree_R1205_FF2_S000001.tmp.pdb          -sequence "aaaguacccuccaagcccuacagguuggaagagggggcuaucaguccuguaggcagacuc"          -secstruct "...((.(((((.....[[[[[[[[......))))).)).......]]]]]]]]......."          -nstruct 3         -out:file:silent /large/otgk/casp/casp16/utils/examples/testR1205_FF2_S000001.out         -minimize_rna true




Basic usage:  /large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease  -fasta <fasta file with sequence>  [ -native <native pdb file> ] 

 Type -help for full slate of options.

********  (C) Copyright Rosetta Commons Member Institutions.  ***************
* Use of Rosetta for commercial purposes may require purchase of a license. *
********  See LICENSE.md or email license@uw.edu for more details. **********
core.init: Checking for fconfig files in pwd and ./rosetta/flags 
core.init: Rosetta version: 2024.15+main.d972b59c53 d972b59c530a12affcbe0eb4a24eedc3ce7d5060 git@github.com:RosettaCommons/rosetta.git 2024-04-02T17:06:29
core.init: Rosetta extras: []
core.init: command: /large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease -s /large/otgk/casp/casp16/utils/examples/test_pfree_R1205_FF2_S000001.tmp.pdb -sequence aaaguacccuccaagcccuacagguuggaagagggggcuaucaguccuguaggcagacuc -secstruct ...((.(((((.....[[[[[[[[......))))).)).......]]]]]]]]....... -nstruct 3 -out:file:silent /large/otgk/casp/casp16/utils/examples/testR1205_FF2_S000001.out -minimize_rna true
basic.random.init_random_generator: 'RNG device' seed mode, using '/dev/urandom', seed=-279197074 seed_offset=0 real_seed=-279197074
basic.random.init_random_generator: RandomGenerator:init: Normal mode, seed=-279197074 RG_type=mt19937
core.init: Resolved executable path: /large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/rna_denovo.default.linuxgccrelease
core.init: Looking for database based on location of executable: /large/otgk/app/rosetta/v2024.15/database/
core.chemical.GlobalResidueTypeSet: Finished initializing fa_standard residue type set.  Created 985 residue types
core.chemical.GlobalResidueTypeSet: Total time to initialize 0.935146 seconds.
core.import_pose.import_pose: File '/large/otgk/casp/casp16/utils/examples/test_pfree_R1205_FF2_S000001.tmp.pdb' automatically determined to be of type PDB
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover:    pdb_seq: aaguacccuccaagcccuacagguuggaagagggggcuaucaguccuguaggcagacuc
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: target_seq: aaaguacccuccaagcccuacagguuggaagagggggcuaucaguccuguaggcagacu
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 3
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target g at 4
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target u at 5
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target a at 6
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target c at 9
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target u at 10
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target c at 12
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 14
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target g at 15
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target c at 18
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target u at 19
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target a at 20
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target c at 21
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 22
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target g at 24
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target u at 26
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target g at 28
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 30
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target g at 31
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 32
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target g at 37
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target c at 38
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target u at 39
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target a at 40
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target u at 41
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target c at 42
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 43
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target g at 44
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target u at 45
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target c at 47
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target u at 48
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target g at 49
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target u at 50
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 51
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target g at 53
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target c at 54
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  g vs target a at 55
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  a vs target g at 56
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target a at 57
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  u vs target c at 58
protocols.rna.denovo.movers.RNA_DeNovoProtocolMover: mismatch in sequence: pdb  c vs target u at 59

ERROR: The sequence in /large/otgk/casp/casp16/utils/examples/test_pfree_R1205_FF2_S000001.tmp.pdb does not match target sequence!!
ERROR:: Exit from: src/protocols/rna/denovo/movers/RNA_DeNovoProtocolMover.cc line: 812

[ ERROR ]: Caught exception:


File: src/protocols/rna/denovo/movers/RNA_DeNovoProtocolMover.cc:812
[ ERROR ] UtilityExitException
ERROR: The sequence in /large/otgk/casp/casp16/utils/examples/test_pfree_R1205_FF2_S000001.tmp.pdb does not match target sequence!!


 ------------------------ Begin developer's backtrace ------------------------- 
BACKTRACE:
/large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/libutility.so(backtrace_string[abi:cxx11](int)+0x5a) [0x7fd73e79135a]
/large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/libutility.so(utility::excn::Exception::Exception(char const*, int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)+0xe0) [0x7fd73e7cf310]
/large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/libutility.so(utility::UtilityExitException::UtilityExitException(char const*, int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)+0x113) [0x7fd73e7964d3]
/large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/libutility.so(utility::exit(char const*, int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, int)+0x3b) [0x7fd73e79611b]
/large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/libprotocols_d.6.so(protocols::rna::denovo::movers::RNA_DeNovoProtocolMover::input_pdb_numbering_setup(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, std::shared_ptr<core::pose::full_model_info::FullModelParameters> const&, utility::vector1<unsigned long, std::allocator<unsigned long> >&, utility::vector1<utility::vector1<unsigned long, std::allocator<unsigned long> >, std::allocator<utility::vector1<unsigned long, std::allocator<unsigned long> > > >&, utility::vector1<unsigned long, std::allocator<unsigned long> >&, utility::vector1<utility::vector1<int, std::allocator<int> >, std::allocator<utility::vector1<int, std::allocator<int> > > >&, utility::vector1<unsigned long, std::allocator<unsigned long> >&, unsigned long&)+0xbc3) [0x7fd740d79143]
/large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/libprotocols_d.6.so(protocols::rna::denovo::movers::RNA_DeNovoProtocolMover::input_numbering_setup(std::tuple<utility::vector1<int, std::allocator<int> >, utility::vector1<char, std::allocator<char> >, utility::vector1<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > const&, std::tuple<utility::vector1<int, std::allocator<int> >, utility::vector1<char, std::allocator<char> >, utility::vector1<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > const&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, std::shared_ptr<core::pose::full_model_info::FullModelParameters> const&, utility::vector1<unsigned long, std::allocator<unsigned long> >&, utility::vector1<utility::vector1<int, std::allocator<int> >, std::allocator<utility::vector1<int, std::allocator<int> > > >&, utility::vector1<unsigned long, std::allocator<unsigned long> >&, utility::vector1<utility::vector1<unsigned long, std::allocator<unsigned long> >, std::allocator<utility::vector1<unsigned long, std::allocator<unsigned long> > > >&, utility::vector1<unsigned long, std::allocator<unsigned long> >&)+0xa6) [0x7fd740d87b46]
/large/otgk/app/rosetta/v2024.15/source/build/src/release/linux/5.15/64/x86/gcc/11/default/libprotocols_d.6.so(protocols::rna::denovo::movers::RNA_DeNovoProtocolMover::de_novo_setup_from_options(utility::options::OptionCollection const&)+0x702) [0x7fd740d88972]
/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease(+0xe998) [0x55f885fc9998]
/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease(+0xff73) [0x55f885fcaf73]
/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease(+0xd1c0) [0x55f885fc81c0]
/lib/x86_64-linux-gnu/libc.so.6(+0x29d90) [0x7fd73e145d90]
/lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80) [0x7fd73e145e40]
/large/otgk/app/rosetta/v2024.15/source/bin/rna_denovo.default.linuxgccrelease(+0xd335) [0x55f885fc8335]
 ------------------------- End developer's backtrace -------------------------- 


AN INTERNAL ERROR HAS OCCURED. PLEASE SEE THE CONTENTS OF ROSETTA_CRASH.log FOR DETAILS.


"""