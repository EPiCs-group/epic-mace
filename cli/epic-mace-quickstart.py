'''Generates example of the input file for epic-mace CLI tool'''

#%% Imports

import sys, os
import argparse


#%% Dummy texts

TEXT_SUBS = '''# name: SMILES
# Alk/Ar
H: [*][H]
Me: [*]C
Et: [*]CC
Pr: [*]CCC
iPr: [*]C(C)C
Bu: [*]CCCC
Ph: [*]c1ccccc1
# *-oxy
OH: [*]O
OMe: [*]OC
OEt: [*]OCC
OPr: [*]OCCC
OiPr: [*]OC(C)C
OBu: [*]OCCCC
OPh: [*]Oc1ccccc1
OAc: [*]OC(=O)C
# amino
NH2: [*]N
NHMe: [*]NC
NMe2: [*]N(C)C
# halogens
F: [*]F
Cl: [*]Cl
Br: [*]Br
I: [*]I
# acceptors
CN: [*]C#N
NO2: [*]N(=O)=O
CO2Me: [*]C(=O)OC
CO2Et: [*]C(=O)OCC
Ac: [*]C(=O)C
'''

TEXT_INPUT = '''# example of epic-mace input file

structure:
  name: Pd_bipy
  geom: SP
  # define structure via ligands & CA
  ligands:
  - "[*]C1=C[N:4]=C(C=C1)C1=[N:3]C=C([*])C=C1 |$_R1;;;;;;;;;;;_R2;;$,c:3,5,13,t:1,8,10|"
  - "[Cl-:2]"
  - "[Cl-:1]"
  CA: "[Pd+2]"
  ## via complex (use either ligands & CA or complex)
  #complex: "[Cl-:1][Pd++]1([Cl-:2])[N:3]2=CC([*])=CC=C2C2=[N:4]1C=C([*])C=C2 |$;;;;;;_R2;;;;;;;;_R1;;$,c:6,8,17,t:3,14,C:11.12,3.2,0.0,2.1|"

substituents:
  file: substituents.yaml # default
  R1:
  - H # if name only, SMILES must be defined in substituents file
  - NMe2
  - OMe
  - OAc
  - SMe "[*]SC" # new substituent
  R2:
  - H
  - CN
  - NO2

stereomer-search:
  regime: no # all, CA, ligands
  drop-enantiomers: true # false
  trans-cycle: -1 # if -1, trans-position for DA-DA donor atoms not allowed
  mer-rule: true # false

conformer-generation:
  num-confs: 10
  rms-threshold: 1.0 # if -1, not rms filtering

'''


#%% Functions

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
        raise RuntimeError(f'\nError: output directory does not exist: {args.path_dir}\n')
    
    return args


def main():
    '''Main function'''
    # get arguments
    try:
        args = read_args()
    except RuntimeError as e:
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


