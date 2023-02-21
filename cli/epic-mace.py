'''Generates 3D coordinates for the given complexes'''

#%% Imports

import re, sys, os
import argparse, yaml

import mace


#%% Preparing arguments

def get_parser():
    '''Generates CLI parser'''
    # subs flags are defined in epilog and parsed from unknown arguments
    parser = argparse.ArgumentParser(
        prog = 'epic-mace', # f'epic-mace, v. {mace.__version__}' # add to description
        description = 'CLI tool for stereomer search and 3D coordinates generation of '
                      'octahedral and square-planar metal complexes. '
                      'For more details see https://epic-mace.readthedocs.io/en/latest/',
        epilog = '  --R1, --R2, etc.      lists of substituent names. If substituent name '
                 'is not described in substituents file, SMILES must be provided (sub="SMILES"). '
                 'Must be specified if complex contains the corresponding substituents',
        formatter_class = argparse.RawTextHelpFormatter
    )
    # input file
    inpf = parser.add_argument_group(title = 'Input file')
    inpf.add_argument(
        '--input', type = str,
        help = 'epic-mace input file; if provided, other oparameters will be ignored'
    )
    # structure
    struct = parser.add_argument_group(
        title = 'Complex structure',
        description = 'Parameters describing the complex structure. --complex '
                      'and (--ligands, --CA) are mutually exclusive, but one of '
                      'these parameter groups are required'
    )
    struct.add_argument(
        '--name', type = str,
        help = 'name of the complex; required if --input is not specified'
    )
    struct.add_argument(
        '--geom', type = str, choices = ['OH', 'SP'],
        help = 'molecular geometry of the central atom:\n'
               '  - OH: octahedral\n'
               '  - SP: square-planar'
               'required if --input is not specified'
    )
    struct.add_argument(
        '--complex', type = str,
        help = 'SMILES of the complex; required if --input is not specified'
    )
    struct.add_argument(
        '--ligands', type = str, nargs = '+',
        help = 'list of ligands\' SMILES; required if --input is not specified'
    )
    struct.add_argument(
        '--CA', type = str,
        help = 'SMILES of the central atom; required if --input is not specified'
    )
    # stereomers
    stereo = parser.add_argument_group(
        title = 'Stereomer search',
        description = 'Parameters of the stereomer search'
    )
    stereo.add_argument(
        '--regime', type = str, default = 'all', choices = ['all', 'CA', 'ligands', 'none'],
        help = 'type of the stereomer search:\n'
               '  - all: iterates over all stereocenters\n'
               '  - ligands: iterates over ligand\'s stereocenters only\n'
               '  - CA: changes configuration of central atom only\n'
               '  - none: do not search for other stereomers'
    )
    stereo.add_argument(
        '--get-enantiomers', action = 'store_true',
        help = 'returns both enantiomers for chiral structures'
    )
    stereo.add_argument(
        '--trans-cycle', default = None, type = int,
        help = 'minimal number of bonds between neighboring donor atoms '
               'required for the trans- spatial arrangement of the donor atoms. '
               'If not specified, such arrangement is considered as impossible'
    )
    stereo.add_argument(
        '--mer-rule', action = 'store_true',
        help = 'applies empirical rule forbidding fac- configuration for '
               'the "rigid" DA-DA-DA fragments of the ligand'
    )
    # conformers
    confs = parser.add_argument_group(
        title = '3D generation',
        description = 'Parameters of generation of 3D atomic coordinates'
    )
    confs.add_argument(
        '--num-confs', type = int, default = 5,
        help = 'number of conformers to generate'
    )
    confs.add_argument(
        '--rms-thresh', type = float, default = None,
        help = 'drops one of two conformers if their RMSD is less than this threshold. '
               'If not specified, all generated conformers are returned'
    )
    # substituents
    subs = parser.add_argument_group(
        title = 'Substituents',
        description = 'Info on substituents to modify the core structure'
    )
    subs.add_argument(
        '--substituents-file', type = str, default = './substituents.yaml',
        help = 'YAML-formatted file containing substituents\' name-SMILES mapping; required '
    )
    
    return parser


def read_subs(unknown):
    '''Extracts Rs from unknown arguments'''
    idxs = [int(arg[3:]) for arg in unknown if re.search('^--R\d+$', arg)] # TODO: R, R0
    # set parser
    subs = argparse.ArgumentParser(prog = 'epic-mace') # prog name for errors
    for idx in sorted(idxs):
        subs.add_argument(f'--R{idx}', type = str, nargs = '+')
    # parse Rs & detect unknown at the same time
    subs_args = subs.parse_args(unknown)
    
    return subs_args


def get_args_from_file(path):
    '''Extracts script parameters from an input file'''
    if not os.path.isfile(path):
        raise ValueError('Error: specified input file does not exist')
    with open(path, 'r') as inpf:
        try:
            args = yaml.safe_load(inpf)
        except Exception as e:
            raise ValueError('Error: bad-formatted input file:\n' + str(e))
    
    return args


def args_to_command(args):
    '''Generates terminal command from the dictionary'''
    cmd = []
    for key, val in args.items():
        if type(val) == list:
            cmd += [f'--{key}'] + [str(_) for _ in val]
        elif type(val) == bool:
            if val: cmd.append(f'--{key}')
        else:
            cmd += [f'--{key}', str(val)]
    
    return cmd


def read_arguments():
    '''Reads arguments from command line or input file'''
    # parse known args
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    # read from input file
    if args.input:
        cmd = args_to_command(get_args_from_file(args.input))
        args, unknown = parser.parse_known_args(cmd)
    # parse substituents
    subs_args = read_subs(unknown)
    
    return {**args.__dict__, **subs_args.__dict__}


def check_arguments():
    
    
    pass



#%% Generation

def get_complex_from_smiles():
    pass


def get_complex_from_ligands():
    pass


def get_stereomers():
    pass


def get_conformers():
    pass


def save_conformers():
    pass


def run_mace_for_system():
    pass



#%% Main function

def main():
    pass



#%% Main code

if __name__ == '__main__':
    
    #main()
    args = read_arguments()
    print(args)


