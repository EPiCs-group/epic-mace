'''
Tests MACE code
'''

#%% Imports

import os, shutil, json

from mace import Complex, ComplexFromLigands


#%% Parameters and data

# paths
path_tests = 'inputs/stereo_tests.json'
path_xyz = 'outputs/stereo_tests'

# load tests
with open(path_tests, 'r') as inpf:
    text = inpf.read()
tests = json.loads(text)


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
        X1 = Complex(test['smiles'], test['geom'],
                          test['maxResonanceStructures'])
        X2 = ComplexFromLigands(test['ligands'], test['CA'], test['geom'],
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


