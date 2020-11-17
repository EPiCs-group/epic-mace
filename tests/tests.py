'''
Tests MACE code
'''

import sys, os, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace


#%% Stereo

# load tests
stereo_tests = json.load(os.path.join(os.path.dirname(__file__), 'inputs/stereo.json'))

# run tests
for test in stereo_tests:
    pass


# # #1: resonance forms, benzoates # TODO: too slow
# smiles = '[Cr+3]12345[O-]C(=O1)c1ccccc1.[O-]2C(=O3)c1ccccc1.[O-]4C(=O5)c1ccccc1 |C:12.14,10.11,3.3,1.0,19.22,21.25|'
# ligands = ['c1ccccc1C(=[O:1])[O-:2]', 'c1ccccc1C(=[O:3])[O-:4]', 'c1ccccc1C(=[O:5])[O-:6]']
# CA = '[Cr+3]'
# geom = 'OH'
# Xs = [mace.Complex(smiles, geom),
#       mace.ComplexFromLigands(ligands, CA, geom)]
# flag = False
# print('Test #1: ', end = '')
# for X in Xs:
#     if len(X.GetStereomers(regime = 'all', dropEnantiomers = True)) != 1:
#         flag = True
#     if len(X.GetStereomers(regime = 'all', dropEnantiomers = False)) != 2:
#         flag = True
# print(f'{"failed" if flag else "success"}\n')


# # #2: resonance forms, aceatates
# smiles = '[Cr+3]12345[O-]C(C)=O1.CC([O-]2)=O3.CC([O-]4)=O5 |C:1.0,4.4,7.7,8.9,11.12,12.14|'
# ligands = ['CC(=[O:1])[O-:2]', 'CC(=[O:3])[O-:4]', 'CC(=[O:5])[O-:6]']
# CA = '[Cr+3]'
# geom = 'OH'
# Xs = [mace.Complex(smiles, geom),
#       mace.ComplexFromLigands(ligands, CA, geom)]
# flag = False
# print('Test #2: ', end = '')
# for X in Xs:
#     if len(X.GetStereomers(regime = 'all', dropEnantiomers = True)) != 1:
#         flag = True
#     if len(X.GetStereomers(regime = 'all', dropEnantiomers = False)) != 2:
#         flag = True
# print(f'{"failed" if flag else "success"}\n')


# #3: resonance forms, symmetric carbene


# 4



#%% 3D generation

# 








