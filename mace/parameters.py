'''
Internal parameters of Complex object
'''

#%% Imports

from itertools import combinations
from copy import deepcopy
from collections import namedtuple

from rdkit.Geometry.rdGeometry import Point3D


#%% Params object

params = namedtuple('ComplexParams', ['FFParams', 'Rcov', 'Syms', 'Geoms',
                                      'Bounds', 'PosVs', 'MinVs', 'EqOrs',
                                      'Nears', 'Angles'])

# force field parameters
params.FFParams = {'X*'     :    2.0,
                   'tether' :  100.0,
                   'kXL'    :  500.0,
                   'kX*'    : 1000.0,
                   'kLA'    :  700.0,
                   'kA*'    :  500.0,
                   'XLA'    :  180.0,
                   'XLO'    :  130.0,
                   'XLA2'   :  120.0,
                   'XLA3'   :  109.5,
                   'kXLA'   :  200.0,
                   'kXLO'   :   80.0,
                   'kALA'   :  200.0,
                   'kZ-LXL' :  200.0,
                   'kE-LXL' :   50.0}

# covalent radii were extracted from CCDC "Elemental Data and Radii" 23.02.2019
# https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx
params.Rcov = [0.23,0.23,1.5,1.28,0.96,0.83,0.68,0.68,0.68,0.64,1.5,1.66,1.41,1.21,1.2,1.05,1.02,
               0.99,1.51,2.03,1.76,1.7,1.6,1.53,1.39,1.61,1.52,1.26,1.24,1.32,1.22,1.22,1.17,
               1.21,1.22,1.21,1.5,2.2,1.95,1.9,1.75,1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.54,
               1.42,1.39,1.39,1.47,1.4,1.5,2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,
               1.94,1.92,1.92,1.89,1.9,1.87,1.87,1.75,1.7,1.62,1.51,1.44,1.41,1.36,1.36,1.32,
               1.45,1.46,1.48,1.4,1.21,1.5,2.6,2.21,2.15,2.06,2,1.96,1.9,1.87,1.8,1.69,1.54,
               1.83,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]
params.Rcov = {i: rcov for i, rcov in enumerate(params.Rcov)}


# SMILES symmetry codes for octahedral geometry @OH1-@OH30
# see OpenSMILES specification for the details
params.Syms = {'OH': { 1: [1,2,3,4,5,6],  2: [1,2,5,4,3,6],  3: [1,2,3,4,6,5],
                       4: [1,2,3,5,4,6],  5: [1,2,3,6,4,5],  6: [1,2,3,5,6,4],
                       7: [1,2,3,6,5,4],  8: [1,2,4,3,5,6],  9: [1,2,4,3,6,5],
                      10: [1,2,5,3,4,6], 11: [1,2,6,3,4,5], 12: [1,2,5,3,6,4],
                      13: [1,2,6,3,5,4], 14: [1,2,4,5,3,6], 15: [1,2,4,6,3,5],
                      16: [1,2,6,4,3,5], 17: [1,2,5,6,3,4], 18: [1,2,6,5,3,4],
                      19: [1,2,4,5,6,3], 20: [1,2,4,6,5,3], 21: [1,2,5,4,6,3],
                      22: [1,2,6,4,5,3], 23: [1,2,5,6,4,3], 24: [1,2,6,5,4,3],
                      25: [1,3,4,5,6,2], 26: [1,3,4,6,5,2], 27: [1,3,5,4,6,2],
                      28: [1,3,6,4,5,2], 29: [1,3,5,6,4,2], 30: [1,3,6,5,4,2]},
               'SP': { 1: [1,2,3,4], 2: [1,2,4,3], 3: [1,3,2,4]}}

# orientations of OH/SP ligands
params.Geoms = {'OH': {'CA': Point3D( 0.0, 0.0, 0.0),
                          1: Point3D( 0.0, 0.0, 2.0),
                          2: Point3D( 2.0, 0.0, 0.0),
                          3: Point3D( 0.0, 2.0, 0.0),
                          4: Point3D(-2.0, 0.0, 0.0),
                          5: Point3D( 0.0,-2.0, 0.0),
                          6: Point3D( 0.0, 0.0,-2.0)},
                'SP': {'CA': Point3D( 0.0, 0.0, 0.0),
                       'X1': Point3D( 0.0, 0.0, 2.0),
                          1: Point3D( 2.0, 0.0, 0.0),
                          2: Point3D( 0.0, 2.0, 0.0),
                          3: Point3D(-2.0, 0.0, 0.0),
                          4: Point3D( 0.0,-2.0, 0.0),
                       'X2': Point3D( 0.0, 0.0,-2.0)}}

# Bounds matrixes
# params for OH and SP
r = 2.0
dr = 0.1
r_max = r + dr # CA-L
r_min = r - dr
r_e_max = 2*r_max # L1..L2, L1-CA-L2 = 180
r_e_min = 2*r_min
r_z_max = 2**0.5 * r_max # L1..L2, L1-CA-L2 = 90
r_z_min = 2**0.5 * r_min
      # CA  # 1 / X1 # 2      # 3      # 4      # 5      # 6 / X2
X = [[  0.0,   r_max,   r_max,   r_max,   r_max,   r_max,   r_max], # CA
     [r_min,     0.0, r_z_max, r_z_max, r_z_max, r_z_max, r_e_max], # 1 / X1
     [r_min, r_z_min,     0.0, r_z_max, r_e_max, r_z_max, r_z_max], # 2
     [r_min, r_z_min, r_z_min,     0.0, r_z_max, r_e_max, r_z_max], # 3
     [r_min, r_z_min, r_e_min, r_z_min,     0.0, r_z_max, r_z_max], # 4
     [r_min, r_z_min, r_z_min, r_e_min, r_z_min,     0.0, r_z_max], # 5
     [r_min, r_e_min, r_z_min, r_z_min, r_z_min, r_z_min,     0.0]] # 6 / X2
# prepare OH
labs = ['CA', 1, 2, 3, 4, 5, 6]
OH = {lab1: {lab2: 0.0 for lab2 in labs} for lab1 in labs}
for (i, lab1), (j, lab2) in combinations(enumerate(labs), r = 2):
    if i > j:
        i, j = j, i
        lab1, lab2 = lab2, lab1
    OH[lab1][lab2] = X[i][j]
    OH[lab2][lab1] = X[j][i]
# prepare SP
labs = ['CA', 'X1', 1, 2, 3, 4, 'X2']
SP = {lab1: {lab2: 0.0 for lab2 in labs} for lab1 in labs}
for (i, lab1), (j, lab2) in combinations(enumerate(labs), r = 2):
    if i > j:
        i, j = j, i
        lab1, lab2 = lab2, lab1
    SP[lab1][lab2] = X[i][j]
    SP[lab2][lab1] = X[j][i]
# add to params
params.Bounds = {'OH': OH, 'SP': SP}

# lists of points corresponding to positive tetrahedra volumes
params.PosVs = {'OH': [['CA',1,2,3],
                       ['CA',1,3,4],
                       ['CA',1,4,5],
                       ['CA',1,5,2],
                       ['CA',6,5,4],
                       ['CA',6,4,3],
                       ['CA',6,3,2],
                       ['CA',6,2,5]],
                'SP': [['CA','X1',1,2],
                       ['CA','X1',2,3],
                       ['CA','X1',3,4],
                       ['CA','X1',4,1],
                       ['CA','X2',4,3],
                       ['CA','X2',3,2],
                       ['CA','X2',2,1],
                       ['CA','X2',1,4]]}

# minimal volumes of polyhedra around CA (to check CA stereo)
params.MinVs = {'OH': 4/3, 'SP': 2/3}

# equivalent orientations of complexes
params.EqOrs = {'OH': {1: [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6],
                       2: [2,3,4,5,1,5,6,3,1,2,6,4,1,3,6,5,1,4,6,2,2,5,4,3],
                       3: [3,4,5,2,5,6,3,1,2,6,4,1,3,6,5,1,4,6,2,1,5,4,3,2],
                       4: [4,5,2,3,6,3,1,5,6,4,1,2,6,5,1,3,6,2,1,4,4,3,2,5],
                       5: [5,2,3,4,3,1,5,6,4,1,2,6,5,1,3,6,2,1,4,6,3,2,5,4],
                       6: [6,6,6,6,4,4,4,4,5,5,5,5,2,2,2,2,3,3,3,3,1,1,1,1]},
                'SP': {1: [1,2,3,4,1,4,3,2],
                       2: [2,3,4,1,4,3,2,1],
                       3: [3,4,1,2,3,2,1,4],
                       4: [4,1,2,3,2,1,4,3]}}
# add enantiomeric orientations for OH (flip 1 and 6 positions)
params.EqOrs['enantOH'] = deepcopy(params.EqOrs['OH'])
params.EqOrs['enantOH'][1] = deepcopy(params.EqOrs['OH'][6])
params.EqOrs['enantOH'][6] = deepcopy(params.EqOrs['OH'][1])

# neigboring ligands positions
params.Nears = {'OH': {1: [2,3,4,5],
                       2: [1,3,6,5],
                       3: [1,2,6,4],
                       4: [1,3,6,5],
                       5: [1,2,6,4],
                       6: [2,3,4,5]},
                'SP': {'X1': [1,2,3,4],
                       1: ['X1',2,'X2',4],
                       2: ['X1',1,'X2',3],
                       3: ['X1',2,'X2',4],
                       4: ['X1',1,'X2',3],
                       'X2': [1,2,3,4]}}

# angles for molecular mechanics
params.Angles = {'OH': {1: {2:  90.0, 3:  90.0, 4:  90.0, 5:  90.0, 6: 180.0},
                        2: {1:  90.0, 3:  90.0, 4: 180.0, 5:  90.0, 6:  90.0},
                        3: {1:  90.0, 2:  90.0, 4:  90.0, 5: 180.0, 6:  90.0},
                        4: {1:  90.0, 2: 180.0, 3:  90.0, 5:  90.0, 6:  90.0},
                        5: {1:  90.0, 2:  90.0, 3: 180.0, 4:  90.0, 6:  90.0},
                        6: {1: 180.0, 2:  90.0, 3:  90.0, 4:  90.0, 5:  90.0}},
                 'SP': {'X1': {1:  90.0, 2:  90.0, 3:  90.0, 4:  90.0, 'X2': 180.0},
                        1: {'X1':  90.0, 2:  90.0, 3: 180.0, 4:  90.0, 'X2':  90.0},
                        2: {'X1':  90.0, 1:  90.0, 3:  90.0, 4: 180.0, 'X2':  90.0},
                        3: {'X1':  90.0, 1: 180.0, 2:  90.0, 4:  90.0, 'X2':  90.0},
                        4: {'X1':  90.0, 1:  90.0, 2: 180.0, 3:  90.0, 'X2':  90.0},
                        'X2': {'X1': 180.0, 1:  90.0, 2:  90.0, 3:  90.0, 4:  90.0}}}


