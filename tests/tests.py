'''
Tests MACE code
'''

import sys, os, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace


# paths
path_tests = os.path.join(os.path.dirname(__file__), 'inputs/full_tests.json')
path_xyz = os.path.join(os.path.dirname(__file__), 'outputs/full_tests')

# load tests
with open(path_tests, 'r') as inpf:
    text = inpf.read()
tests = json.loads(text)

# testing
idxs = [10]
tests = [test for i, test in enumerate(tests) if i in idxs]

# run tests
for test in tests:
    print(f'#{test["num"]} {test["name"]}:')
    # initialize complex
    print('... complex initialization: ', end = '')
    try:
        X1 = mace.Complex(test['smiles'], test['geom'])
        X2 = mace.ComplexFromLigands(test['ligands'], test['CA'], test['geom'])
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
                                 dropEnantiomers = False)
            if len(Xs) != test['numEnantiomers']:
                flag = True
            Xs = X.GetStereomers(regime = test['regime'],
                                 minTransCycle = test['minTransCycle'],
                                 dropEnantiomers = True)
            if len(Xs) != test['numStereomers']:
                flag = True
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print('error\n')
        continue
    if flag:
        print('failed')
        continue
    print('success')
    # make 3D
    if not test['3D']:
        continue
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
        if not os.path.isdir(f'{path_xyz}/{test["num"]}'):
            os.mkdir(f'{path_xyz}/{test["num"]}')
        X.ToXYZ(f'{path_xyz}/{test["num"]}/{i}_min.xyz', confId = -1)
        X.ToXYZ(f'{path_xyz}/{test["num"]}/{i}_all.xyz', confId = -2)
    if not errors and not no_geom:
        print('success')
    else:
        print(f'{success} OK, {no_geom} no 3D, {errors} errors')
    # empty line
    print()



# Testing

'''

test = tests[6]
X1 = mace.Complex(test['smiles'], test['geom'])
X2 = mace.ComplexFromLigands(test['ligands'], test['CA'], test['geom'])
Xs = X1.GetStereomers(regime = test['regime'],
                      minTransCycle = test['minTransCycle'],
                      dropEnantiomers = False)
Xs = X2.GetStereomers(regime = test['regime'],
                      minTransCycle = test['minTransCycle'],
                      dropEnantiomers = False)
X = Xs[0]
X.AddConformers(numConfs = 20, rmsThresh = 2.0)


X = mace.Complex('[O+]#[C-:2][Ir]([C-:3]#[O+])[C-:1]#[O+] |C:1.1,5.4,3.2|', 'OH')
print(X.AddConformers(numConfs = 20, rmsThresh = 2.0))





'''