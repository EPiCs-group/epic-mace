'''
Tests MACE code
'''

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace

import pandas as pd
from itertools import product


# paths
path_dir = os.path.dirname(__file__)
path_ls = os.path.join(path_dir, 'inputs/ligands.csv')
path_xyz = os.path.join(path_dir, 'outputs/brute_xyz')
out_name = '_temp.csv'

# ligands
ls = pd.read_csv(path_ls)
ls1 = ls.name[ls.n == 1]
ls2 = ls.name[ls.n == 2]
ls3 = ls.name[ls.n == 3]

# OH B1A4
subdir = 'OH_B1_M4'
CA = '[Ru+2]'
geom = 'OH'
systems = {}
for l1, l2 in product(ls1, ls2):
    basename = f'{geom}_{l2}_{l1}.4'
    ligands = [ls.smiles[ls.name == l2].iloc[0]] + [ls.smiles[ls.name == l1].iloc[0]]*4
    systems[basename] = ligands



# check subdir
if not os.path.isdir(f'{path_xyz}/{subdir}'):
    os.mkdir(f'{path_xyz}/{subdir}')

# generate systems
output = ['complex,success,no3D,error']
for basename, ligands in systems.items():
    # stereomers
    X = mace.ComplexFromLigands(ligands, CA, geom)
    Xs = X.GetStereomers()
    # 3D
    result = {'err': 0, 'no3D': 0, 'success': 0}
    for i, X in enumerate(Xs):
        try:
            X.AddConformers(numConfs = 5, maxAttempts = 30)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            result['err'] += 1
            continue
        if X.GetNumConformers() == 0:
            result['no3D'] += 1
            continue
        result['success'] += 1
        X.ToXYZ(f'{path_xyz}/{subdir}/{basename}_{i}.xyz', -1)
    output.append(f'{basename},{result["success"]},{result["no3D"]},{result["err"]}')

# save output
with open(f'{path_xyz}/{subdir}/{out_name}', 'w') as outf:
    outf.write('\n'.join(output)+'\n')


