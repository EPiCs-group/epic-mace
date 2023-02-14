'''
Tests MACE code
'''

#%% Imports

import os, time

from mace import ComplexFromLigands


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



#%% Read ligands

# paths
path_ls = 'ligands/cleared_ligands/csd_ligands.csv'
path_xyz = 'outputs'

# ligands
ls = read_ligands(path_ls)
auxs = {'H': '[H-:1]', 'Cl': '[Cl-:1]', 'CO': '[C-:1]#[O+]'}


#%% Prepare systems

# params
sp_params = {'geom': 'SP', 'CA': '[Pd+2]',
             'regime': 'CA', 'mTS': 12, 'mRS': 1, 'mer': True,
             'numConfs': 5, 'maxAttempts': 10}
oh_params = {'geom': 'OH', 'CA': '[Ru+2]',
             'regime': 'CA', 'mTS': None, 'mRS': 1, 'mer': True,
             'numConfs': 3, 'maxAttempts': 10}
systems = []

# SP
for n, name, smiles in ls['SP']:
    systems.append( (name, n, [smiles] + [auxs['CO']]*(4-n), sp_params, 'SP_CO') )

# OH
for n, name, smiles in ls['OH']:
    systems.append( (name, n, [smiles] + [auxs['CO']]*(6-n), oh_params, 'OH_CO') )


#%% Main code

# check subdirs
subdirs = ['logs'] + list(set([_[4] for _ in systems]))
for subdir in subdirs:
    if not os.path.isdir(f'{path_xyz}/{subdir}'):
        os.mkdir(f'{path_xyz}/{subdir}')

# output initialization
for subdir in set([_[-1] for _ in systems]):
    with open(f'{path_xyz}/logs/{subdir}_temp.csv', 'w') as outf:
        outf.write('complex,n,success,no3D,error,time\n')

# generate systems
for i, (basename, n, ligands, params, subdir) in enumerate(systems):
    print(f'{i:>4} of {len(systems)}')
    # params
    geom = params['geom']
    CA = params['CA']
    regime = params['regime']
    maxResonanceStructures = params['mRS']
    minTransCycle = params['mTS']
    merRule = params['mer']
    numConfs = params['numConfs']
    maxAttempts = params['maxAttempts']
    result = {'err': 0, 'no3D': 0, 'success': 0}
    # stereomers
    t1 = time.time()
    X = ComplexFromLigands(ligands, CA, geom, maxResonanceStructures)
    Xs = X.GetStereomers(regime = regime, minTransCycle = minTransCycle, merRule = merRule)
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
        X.ToXYZ(f'{path_xyz}/{subdir}/{basename}_{i}.xyz', -1)
    T = time.time() - t1
    with open(f'{path_xyz}/logs/{subdir}_temp.csv', 'a') as outf:
        outf.write(f'{basename},{n},{result["success"]},{result["no3D"]},{result["err"]},{T:.1f}\n')


