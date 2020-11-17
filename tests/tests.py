'''
Tests MACE code
'''

import sys, os, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace


# paths
path_tests = os.path.join(os.path.dirname(__file__), 'inputs/stereo.json')
path_xyz = os.path.join(os.path.dirname(__file__), 'xyz')

# load tests
with open(path_tests, 'r') as inpf:
    text = inpf.read()
stereo_tests = json.loads(text)

# run tests
for test in stereo_tests:
    print(f'#{test["num"]} {test["name"]}: ', end = '')
    # initialize complex
    try:
        X1 = mace.Complex(test['smiles'], test['geom'])
        X2 = mace.ComplexFromLigands(test['ligands'], test['CA'], test['geom'])
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print('complex initialization error')
    # make stereomers
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
        if flag:
            print("failed")
            continue
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print('stereomer search error')
    # make 3D
    if not test['3D']:
        continue
    for X in Xs:
        try:
            X.AddConformers(numConfs = 5, rmsThresh = 2.0)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            print('3D generation error')
        if X.GetNumConformers() == 0:
            print('success but no 3D')
            continue
        if not os.path.isdir(f'{path_xyz}/{test["num"]}'):
            os.mkdir(f'{path_xyz}/{test["num"]}')
        X.ToXYZ(f'{path_xyz}/{test["num"]}/min.xyz', confId = -2)
        X.ToXYZ(f'{path_xyz}/{test["num"]}/all.xyz', confId = -1)
        print('success')


