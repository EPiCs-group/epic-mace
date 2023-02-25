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
        prog = 'epic-mace', # f'epic-mace, v. {mace.__version__}' # TODO: add to description
        description = 'CLI tool for stereomer search and 3D coordinates generation of '
                      'octahedral and square-planar metal complexes. '
                      'For more details see https://epic-mace.readthedocs.io/en/latest/',
        epilog = '  --R1, --R2, etc.      lists of substituent names. If substituent name '
                 'is not described in substituents file, SMILES must be provided (sub="SMILES"). '
                 'Must be specified if complex contains the corresponding substituents',
        formatter_class = argparse.RawTextHelpFormatter
    )
    # input/output files
    parser.add_argument(
        'out_dir', type = str, nargs = '?', default = './',
        help = 'directory to store epic-mace output'
    )
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
    struct.add_argument(
        '--res-structs', type = int, default = 1,
        help = 'maximal number of resonance structures'
    )
    # stereomers
    stereo = parser.add_argument_group(
        title = 'Stereomer search',
        description = 'Parameters of the stereomer search'
    )
    stereo.add_argument(
        '--regime', default = 'all', choices = ['all', 'CA', 'ligands', 'none'],
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
    idxs = [int(arg[3:]) for arg in unknown if re.search('^--R\d+$', arg)]
    if 0 in idxs:
        raise ValueError('Input error: --R0 substituent is forbidden; use --R1, etc.')
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
            raise ValueError('Input error: bad-formatted input file:\n' + str(e))
    
    return args


def args_to_command(args):
    '''Generates terminal command from the dictionary'''
    cmd = [args['out_dir']]
    for key, val in args.items():
        if key == 'out_dir':
            continue
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


def get_sub_idxs(smis):
    '''Returns list of substituents' indexes used in a complex (R1, R2, etc)'''
    idxs = []
    for smi in smis:
        mol = mace.MolFromSmiles(smi)
        for a in mol.GetAtoms():
            if not a.GetAtomicNum() and not a.GetAtomMapNum() and a.GetIsotope():
                idxs.append(a.GetIsotope())
    
    return sorted(idxs)


def check_sub(smiles, name):
    '''Checks if SMILES satisfies requirements for substitients (contains
    exactly one single-bonded dummy atom)
    
    Arguments:
        smiles (str): SMILES of the substituent
        name (str): name of the substituent (for the error message)
    '''
    try:
        mol = mace.MolFromSmiles(smiles)
    except Exception: # as e:
        raise ValueError(f'Bad substituent {name}: {smiles}\nCannot read the SMILES string')
    dummies = [_ for _ in mol.GetAtoms() if _.GetSymbol() == '*']
    if len(dummies) != 1:
        raise ValueError(f'Bad substituent {name}: {smiles}\nSubstituent must contain exactly one dummy atom')
    if len(dummies[0].GetNeighbors()) != 1:
        raise ValueError(f'Bad substituent {name}: {smiles}\nSubstituent\'s dummy must be bonded to exactly one atom by single bond')
    if str(dummies[0].GetBonds()[0].GetBondType()) != 'SINGLE':
        raise ValueError(f'Bad substituent {name}: {smiles}\nSubstituent\'s dummy must be bonded to exactly one atom by single bond')
    
    return


def get_subs(args, struct):
    '''Checks substituents and returns them as addend to parameters'''
    # compare Rs in args and stars in struct
    isotopes = set()
    smis = args[struct] if struct == 'ligands' else [args['complex']]
    for smiles in smis:
        mol = mace.MolFromSmiles(smiles)
        for a in mol.GetAtoms():
            if not a.GetAtomicNum() and not a.GetAtomMapNum() and a.GetIsotope():
                isotopes.add(a.GetIsotope())
    Rs_mol = set([f'R{i}' for i in isotopes])
    Rs_inp = set([k for k in args if re.search('^R\d+$', k)])
    if Rs_mol.difference(Rs_inp):
        diff = ', '.join(Rs_mol.difference(Rs_inp))
        msg = f'some substituents are not defined in the input: {{{diff}}}'
        raise ValueError(f'Input error: {msg}')
    if Rs_inp.difference(Rs_mol):
        diff = ', '.join(Rs_inp.difference(Rs_mol))
        msg = f'excessive substituents in the input: {{{diff}}}'
        raise ValueError(f'Input error: {msg}')
    # check if empty
    Rs = sorted(list(Rs_mol))
    if not Rs:
        return {}
    # read a file
    path = args['substituents_file']
    if not os.path.exists(path):
        raise ValueError(f'Input error: substituent file does not exist: {path}')
    with open(path, 'r') as inpf:
        try:
            subs_info = yaml.safe_load(inpf)
        except Exception as e:
            raise ValueError('Input error: bad-formatted substituents file:\n' + str(e))
    # prepare subs
    outp = {}
    for R in Rs:
        subs = args[R]
        addend = []
        if len(subs) != len(set(subs)):
            raise ValueError(f'Input error: repeating substituent in --{R}')
        for sub in subs:
            if sub not in subs_info:
                raise ValueError(f'Input error: missing substituent {sub}, please define it in substituents file')
            smiles = subs_info[sub]
            check_sub(smiles, sub)
            addend.append( (sub, smiles) )
        outp[R] = addend
    
    return outp


def check_arguments(args):
    '''Checks prepared arguments'''
    params = {}
    # check output dir
    if not os.path.isdir(args['out_dir']):
        raise ValueError('Input error: specified output directory does not exist')
    # check required
    for key in ('name', 'geom'):
        if not args[key]:
            raise ValueError(f'Input error: --{key} is not specified')
        params[key] = args[key]
    # check mutually exclusive
    has_complex = bool(args['complex'])
    has_ligands = args['ligands'] and args['CA']
    if not has_complex and not has_ligands:
        raise ValueError('Input error: neither --complex nor --ligands and --CA are not specified')
    if has_complex and has_ligands:
        raise ValueError('Input error: both --complex and (--ligands, --CA) are specified')
    # check SMILES
    if has_complex:
        # --complex
        try:
            smiles = args['complex']
            _ = mace.MolFromSmiles(smiles)
        except Exception: # as e:
            raise ValueError(f'Input error: --complex: unreadable SMILES: {smiles}')
        params['complex'] = args['complex']
    elif has_ligands:
        # --ligands
        try:
            for smiles in args['ligands']:
                _ = mace.MolFromSmiles(smiles)
        except Exception: # as e:
            raise ValueError(f'Input error: --ligands: unreadable SMILES: {smiles}')
        params['ligands'] = args['ligands']
        # --CA
        try:
            smiles = args['CA']
            _ = mace.MolFromSmiles(smiles)
        except Exception: # as e:
            raise ValueError(f'Input error: --CA: unreadable SMILES: {smiles}')
        params['CA'] = args['CA']
    # no need to check
    for key in ('regime', 'get_enantiomers', 'mer_rule'):
        params[key]  = args[key]
    # check numeric params
    if args['trans_cycle'] is not None:
        if args['trans_cycle'] < 1:
            raise ValueError('InputError: --trans-cycle must be a positive integer')
    if args['num_confs'] < 1:
        raise ValueError('InputError: --num-confs must be a positive integer')
    if args['res_structs'] < 1:
        raise ValueError('InputError: --res-structs must be a positive integer')
    if args['rms_thresh'] < 0:
        raise ValueError('InputError: --rms-thresh must be a positive real number')
    for key in ('trans_cycle', 'num_confs', 'rms_thresh'):
        params[key] = args[key]
    # get subs
    struct = 'complex' if 'complex' in params else 'ligands'
    subs = get_subs(args, struct)
    if subs:
        params['subs'] = subs
    
    return params



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

def _main():
    '''Generates 3D coordinates for all stereomers of input complexes'''
    args = read_arguments()
    check_arguments(args)
    print(args)
    
    return


def main():
    '''Main function (error handler for _main)'''
    try:
        _main()
    except Exception as e:
        print(e)
        sys.exit()
    
    return



#%% Main code

if __name__ == '__main__':
    
    _main()

