'''
Tests MACE code
'''

#%% Imports

import sys, os, shutil, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace


#%% Parameters and data

# paths
path_tests = os.path.join(os.path.dirname(__file__), 'inputs/stereo_tests.json')
path_xyz = os.path.join(os.path.dirname(__file__), 'outputs/stereo_tests')

# load tests
with open(path_tests, 'r') as inpf:
    text = inpf.read()
tests = json.loads(text)

# # testing
# idxs = [10]
# tests = [test for idx, test in enumerate(tests) if idx in idxs]


#%% Run tests

# recreate directory
if os.path.exists(path_xyz):
    shutil.rmtree(path_xyz)
os.mkdir(path_xyz)

# run tests
for idx, test in enumerate(tests):
    print(f'#{idx} {test["name"]}:')
    # initialize complex
    print('... complex initialization: ', end = '')
    try:
        X1 = mace.Complex(test['smiles'], test['geom'],
                          test['maxResonanceStructures'])
        X2 = mace.ComplexFromLigands(test['ligands'], test['CA'], test['geom'],
                                     test['maxResonanceStructures'])
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print('error\n')
        continue
    print('success')
    # check equivalence
    print('... complex equivalence: ', end = '')
    if not X1.IsEqual(X2):
        print('error\n')
        continue
    print('success')
    # make stereomers
    print('... stereomer search: ', end = '')
    try:
        flag = False
        for X in (X1, X2):
            Xs = X.GetStereomers(regime = test['regime'],
                                 minTransCycle = test['minTransCycle'],
                                 dropEnantiomers = False,
                                 merRule = test['merRule'])
            if len(Xs) != test['numEnantiomers']:
                flag = True
            Xs = X.GetStereomers(regime = test['regime'],
                                 minTransCycle = test['minTransCycle'],
                                 dropEnantiomers = True,
                                 merRule = test['merRule'])
            if len(Xs) != test['numStereomers']:
                flag = True
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print('error\n')
        continue
    if flag:
        print('failed\n')
        continue
    print('success')
    # 3D
    print('... 3D generation: ', end = '')
    errors = 0
    no_geom = 0
    success = 0
    for i, X in enumerate(Xs):
        try:
            X.AddConformers(numConfs = 5, rmsThresh = 2.0)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            errors += 1
            continue
        if X.GetNumConformers() == 0:
            no_geom += 1
            continue
        success += 1
        X.ToXYZ(f'{path_xyz}/{idx}_{i}.xyz', confId = -1)
    if not errors and not no_geom:
        print('success\n')
    else:
        print(f'{success} OK, {no_geom} no 3D, {errors} errors\n')



# # RIP good exp test. Nothing personal, you just wanted too much
# ''' {"name": "(A-B*-C)D2F",
#   "smiles": "C1C[NH:3]2CC[P:2](c3ccccc3)(c3ccccc3)[Ru+2]2([C-:4]2=[N+]1C=CN2c1c(C)cc(cc1C)C)([C-:1]#[O+])([H-:6])[H-:5] |c:22,25,C:2.20,5.19,19.21,35.40,36.41,33.38|",
#   "ligands": ["c1ccccc1[P:2](c1ccccc1)CC[NH:3]CC[N+]=1C=CN(c2c(C)cc(C)cc(C)2)[C-:4]1", "[O+]#[C-:1]", "[H-:5]", "[H-:6]"],
#   "CA": "[Ru+2]", "geom": "OH", "maxResonanceStructures": null,
#   "regime": "all", "minTransCycle": null,
#   "numStereomers": 9, "numEnantiomers": 18}'''

# X = mace.ComplexFromLigands(['c1c(C[NH2:1])[n:2]c(CC[NH2:3])cc1', '[Cl-:4]', '[Cl-:5]', '[Cl-:6]'], '[Ru+2]', 'OH')

