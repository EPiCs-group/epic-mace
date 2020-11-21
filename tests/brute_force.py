'''
Tests MACE code
'''

#%% Imports

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace


#%% Functions

def read_ligands(path):
    '''
    Reads and systemize CSD ligands info (duct taping pandas
    and minimizing required libs)
    '''
    # read file
    with open(path, 'r') as inpf:
        text = [_.strip() for _ in inpf.readlines()]
        text = [_ for _ in text[1:] if _]
    # parse file
    info = {'SP': [], 'OH': []}
    for line in text:
        name, geom, n, smiles, refcodes = line.split(',')
        n = int(n)
        info[geom].append( (n, name, smiles) )
    
    return info


def generate_3D(systems, params):
    '''
    Runs EPiC MACE for bunch of systems
    '''
    # parameters
    path = params['path']
    subdir = params['subdir']
    CA = params['CA']
    geom = params['geom']
    regime = params['regime']
    minTransCycle = params['mTS']
    numConfs = params['numConfs']
    maxAttempts = params['maxAttempts']
    maxResonanceStructures = params['maxResonanceStructures']
    
    # check subdir
    if not os.path.isdir(f'{path}/{subdir}'):
        os.mkdir(f'{path}/{subdir}')
    
    # prepare output
    with open(f'{path}/{subdir}_temp.csv', 'w') as outf:
        outf.write('complex,success,no3D,error\n')
    
    # generate systems
    for i, (basename, ligands) in enumerate(systems.items()):
        print(f'{i:>4} of {len(systems)}')
        result = {'err': 0, 'no3D': 0, 'success': 0}
        # stereomers
        X = mace.ComplexFromLigands(ligands, CA, geom, maxResonanceStructures)
        Xs = X.GetStereomers(regime = regime, minTransCycle = minTransCycle)
        # 3D
        for i, X in enumerate(Xs):
            try:
                X.AddConformers(numConfs = numConfs, maxAttempts = maxAttempts)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                result['err'] += 1
                continue
            if X.GetNumConformers() == 0:
                result['no3D'] += 1
                continue
            result['success'] += 1
            X.ToXYZ(f'{path}/{subdir}/{basename}_{i}.xyz', -1)
        with open(f'{path}/{subdir}_temp.csv', 'a') as outf:
            outf.write(f'{basename},{result["success"]},{result["no3D"]},{result["err"]}\n')
    
    return



#%% Read ligands

# paths
path_dir = os.path.dirname(__file__)
path_ls = os.path.join(path_dir, 'ligands/cleared_ligands/csd_ligands.csv')
path_xyz = os.path.join(path_dir, 'outputs/brute_xyz')

# ligands
ls = read_ligands(path_ls)
aux_ls = {'H': '[H-:1]', 'Cl': '[Cl-:1]', 'CO': '[C-:1]#[O+]'}



#%% pure SP ligands

# # parameters
# params = {'path'  : path_xyz, 'subdir': 'SP',
#           'CA'    : '[Pd+2]', 'geom'  : 'SP',
#           'regime': 'CA', 'mTS'   : 12,
#           'numConfs': 5, 'maxAttempts': 10, 'maxResonanceStructures': 1}

# # ligands
# systems = {name: [smiles] for n, name, smiles in ls['SP']}

# # run
# generate_3D(systems, params)



#%% pure OH ligands

# parameters
params = {'path'  : path_xyz, 'subdir': 'OH',
          'CA'    : '[Ru+2]', 'geom'  : 'OH',
          'regime': 'CA', 'mTS'   : None,
          'numConfs': 3, 'maxAttempts': 10, 'maxResonanceStructures': 1}

# ligands
systems = {name: [smiles] for n, name, smiles in ls['OH']}

# run
generate_3D(systems, params)



#%% SP simple complexes

# # params
# params = {'path'  : path_xyz,
#           'subdir': 'SP',
#           'CA'    : '[Pd+2]',
#           'geom'  : 'SP',
#           'regime': 'CA',
#           'mTS'   : None}

# # ligands
# systems = {}
# for n in ls['SP']:
#     for name, smiles in ls['SP'][n]:
#         if n == 4:
#             basename = name
#             ligands = [smiles]
#             systems[basename] = ligands
#             continue
#         for aux in aux_ls:
#             basename = f'{name}_{aux}'
#             ligands = [smiles] + [aux_ls[aux]]*(4-n)
#             systems[basename] = ligands






