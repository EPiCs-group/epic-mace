'''
Tests MACE code
'''

#%% Imports

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace

import multiprocessing as mp
from multiprocessing.pool import ThreadPool
import time


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
        outf.write('complex,success,no3D,error,time\n')
    
    # generate systems
    for i, (basename, ligands) in enumerate(systems.items()):
        print(f'{i:>4} of {len(systems)}')
        result = {'err': 0, 'no3D': 0, 'success': 0}
        # stereomers
        t1 = time.time()
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
        T = time.time() - t1
        with open(f'{path}/{subdir}_temp.csv', 'a') as outf:
            outf.write(f'{basename},{result["success"]},{result["no3D"]},{result["err"]},{T:.1f}\n')
    
    return



#%% Read ligands

# paths
path_dir = os.path.dirname(__file__)
path_ls = os.path.join(path_dir, 'ligands/cleared_ligands/csd_ligands.csv')
path_xyz = os.path.join(path_dir, 'outputs/brute_xyz')

# params
timelim = 60

# ligands
ls = read_ligands(path_ls)
auxs = {'H': '[H-:1]', 'Cl': '[Cl-:1]', 'CO': '[C-:1]#[O+]'}


#%% Main code

if __name__ == '__main__':
    
    # params
    sp_params = {'geom': 'SP', 'CA': '[Pd+2]',
                 'regime': 'CA', 'mTS': 12, 'mRS': 1,
                 'numConfs': 5, 'maxAttempts': 10}
    oh_params = {'geom': 'OH', 'CA': '[Ru+2]',
                 'regime': 'CA', 'mTS': None, 'mRS': 1,
                 'numConfs': 3, 'maxAttempts': 10}
    systems = []
    
    # SP
    for n, name, smiles in ls['SP']:
        systems.append( (name, n, [smiles], sp_params, 'SP') )
        systems.append( (name, n, [smiles] + [auxs['H']]*(4-n), sp_params, 'SP_H') )
        systems.append( (name, n, [smiles] + [auxs['Cl']]*(4-n), sp_params, 'SP_Cl') )
        systems.append( (name, n, [smiles] + [auxs['CO']]*(4-n), sp_params, 'SP_CO') )
    
    # OH
    for n, name, smiles in ls['OH']:
        systems.append( (name, n, [smiles], oh_params, 'OH') )
        systems.append( (name, n, [smiles] + [auxs['H']]*(6-n), oh_params, 'OH_H') )
        systems.append( (name, n, [smiles] + [auxs['Cl']]*(6-n), oh_params, 'OH_Cl') )
        systems.append( (name, n, [smiles] + [auxs['CO']]*(6-n), oh_params, 'OH_CO') )
    
    # check subdirs
    subdirs = ['logs'] + list(set([_[4] for _ in systems]))
    for subdir in subdirs:
        if not os.path.isdir(f'{path_xyz}/{subdir}'):
            os.mkdir(f'{path_xyz}/{subdir}')
    
    # generate systems
    #pool = ThreadPool(1)
    for i, (basename, n, ligands, params, subdir) in enumerate(systems):
        if i < 34:
            continue
        print(f'{i:>4} of {len(systems)}')
        # params
        geom = params['geom']
        CA = params['CA']
        regime = params['regime']
        maxResonanceStructures = params['mRS']
        minTransCycle = params['mTS']
        numConfs = params['numConfs']
        maxAttempts = params['maxAttempts']
        result = {'err': 0, 'no3D': 0, 'success': 0, 'timelim': 0}
        # output initialization
        if i == 0:
            with open(f'{path_xyz}/logs/{subdir}_temp.csv', 'w') as outf:
                outf.write('complex,n,success,no3D,error,timelim,time\n')            
        # stereomers
        t1 = time.time()
        X = mace.ComplexFromLigands(ligands, CA, geom, maxResonanceStructures)
        Xs = X.GetStereomers(regime = regime, minTransCycle = minTransCycle)
        # 3D
        for i, X in enumerate(Xs):
            pool = ThreadPool(1)
            try:
                res = pool.apply_async(X.AddConformers,
                                       kwds = {'numConfs': numConfs,
                                               'maxAttempts': maxAttempts})
                res.get(timelim)
                #X.AddConformers(numConfs = numConfs, maxAttempts = maxAttempts)
            except mp.TimeoutError:
                result['timelim'] += 1
                pool.terminate()
                continue
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
            outf.write(f'{basename},{result["success"]},{result["no3D"]},{result["err"]},{result["timelim"]},{T:.1f}\n')



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






