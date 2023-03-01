'''Generates example of the input file for epic-mace CLI tool'''

#%% Imports

import sys, os
import argparse


#%% Dummy texts

TEXT_SUBS = '''# name: SMILES
# Alk/Ar
H: "[*][H]"
Me: "[*]C"
Et: "[*]CC"
Pr: "[*]CCC"
iPr: "[*]C(C)C"
Bu: "[*]CCCC"
Ph: "[*]c1ccccc1"
# *-oxy
OH: "[*]O"
OMe: "[*]OC"
OEt: "[*]OCC"
OPr: "[*]OCCC"
OiPr: "[*]OC(C)C"
OBu: "[*]OCCCC"
OPh: "[*]Oc1ccccc1"
OAc: "[*]OC(=O)C"
# amino
NH2: "[*]N"
NHMe: "[*]NC"
NMe2: "[*]N(C)C"
# halogens
F: "[*]F"
Cl: "[*]Cl"
Br: "[*]Br"
I: "[*]I"
# acceptors
CN: "[*]C#N"
NO2: "[*]N(=O)=O"
CO2Me: "[*]C(=O)OC"
CO2Et: "[*]C(=O)OCC"
Ac: "[*]C(=O)C"
'''

TEXT_INPUT = '''# example of epic-mace input file

# output directory
out_dir: ./

# structure
name: RhCl_MeCN_bipy
geom: SP
res-structs: 1
# define structure via ligands & CA
ligands:
  - "[*]C1=C[N:4]=C(C=C1)C1=[N:3]C=C([*])C=C1 |$_R1;;;;;;;;;;;_R2;;$,c:3,5,13,t:1,8,10|"
  - "[N:2]#CC"
  - "[Cl-:1]"
CA: "[Rh+]"
## via complex (use either ligands & CA or complex)
#complex: "[Cl-:1][Rh+]1([N:2]#CC)[N:3]2=CC([*])=CC=C2C2=[N:4]1C=C([*])C=C2 |$;;;;;;;;_R2;;;;;;;;_R1;;$,c:8,10,19,t:5,16,C:0.0,2.1,5.4,13.14|"

# stereomer-search
regime: all # all, CA, ligands, none
get-enantiomers: false # true
trans-cycle: no # if no, trans-position for DA-DA donor atoms not allowed
mer-rule: true # false

# conformer-generation
num-confs: 3
rms-thresh: 1.0 # if -1, not rms filtering

# conformer post-processing
num-repr-confs: no # or positive integer
e-rel-max: 25.0 # kJ/mol
drop-close-energy: true # false

# substituents
substituents-file: substituents.yaml # default
R1: # name: SMILES must be defined in substituents file
  - H
  - NMe2
  - OMe
R2:
  - H
  - CN
  - NO2

'''


#%% Functions

class MaceInputError(Exception):
    '''Custom exception for capturing the known errors'''
    def __init__(self, message):
        super().__init__('Input error: ' + message)


def read_args():
    '''Reads CLI arguments'''
    # parser
    parser = argparse.ArgumentParser(
        prog = 'epic-mace-input',
        description = 'Generates example of the input file for epic-mace CLI tool'
    )
    parser.add_argument(
        'path_dir', type = str, default = './', nargs = '?',
        help = 'Directory of epic-mace project'
    )
    # get arguments
    args = parser.parse_args()
    if not os.path.isdir(args.path_dir):
        raise MaceInputError(f'Output directory does not exist: {args.path_dir}')
    
    return args


def main():
    '''Main function'''
    # get arguments
    try:
        args = read_args()
    except MaceInputError as e:
        print(e)
        sys.exit()
    # save files
    with open(os.path.join(args.path_dir, 'substituents.yaml'), 'w') as outf:
        outf.write(TEXT_SUBS)
    with open(os.path.join(args.path_dir, 'mace_input.yaml'), 'w') as outf:
        outf.write(TEXT_INPUT)
    
    return


#%% Main code

if __name__ == '__main__':
    
    main()


