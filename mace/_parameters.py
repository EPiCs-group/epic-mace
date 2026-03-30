'''Contains symmetric and geometric parameters required for stereomer generation
and 3D embedding
'''

#%% Imports

from itertools import combinations, permutations
from copy import deepcopy
from collections import namedtuple
import numpy as np

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


def _point_to_array(point):
    """Convert Point3D-like objects to NumPy arrays."""
    return np.array([point.x, point.y, point.z], dtype=float)


def _build_bounds(geom):
    """Build a bounds matrix from idealized donor coordinates."""
    labels = ['CA'] + [lab for lab in geom if lab != 'CA']
    bounds = {lab1: {lab2: 0.0 for lab2 in labels} for lab1 in labels}
    r = 2.0
    dr = 0.1
    scale_min = (r - dr) / r
    scale_max = (r + dr) / r
    arrays = {lab: _point_to_array(point) for lab, point in geom.items()}
    for lab in labels:
        if lab == 'CA':
            continue
        dist = np.linalg.norm(arrays[lab] - arrays['CA'])
        bounds['CA'][lab] = dist + dr
        bounds[lab]['CA'] = max(dist - dr, 0.1)
    donors = [lab for lab in labels if lab != 'CA']
    for lab1, lab2 in combinations(donors, r=2):
        dist = np.linalg.norm(arrays[lab1] - arrays[lab2])
        bounds[lab1][lab2] = dist * scale_max
        bounds[lab2][lab1] = dist * scale_min

    return bounds


def _build_angles_and_nears(geom):
    """Build central-atom angle and adjacency tables from idealized coordinates."""
    arrays = {lab: _point_to_array(point) for lab, point in geom.items() if lab != 'CA'}
    labels = [lab for lab in geom if lab != 'CA']
    angles = {lab: {} for lab in labels}
    nears = {lab: [] for lab in labels}
    for lab1, lab2 in combinations(labels, r=2):
        v1 = arrays[lab1]
        v2 = arrays[lab2]
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        angle = float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
        angles[lab1][lab2] = angle
        angles[lab2][lab1] = angle
        if angle < 179.9:
            nears[lab1].append(lab2)
            nears[lab2].append(lab1)

    return angles, nears


def _get_symmetry_permutations(geom):
    """Return proper and improper symmetry permutations for donor positions."""
    labels = sorted(lab for lab in geom if str(lab).isdigit())
    arrays = {lab: _point_to_array(geom[lab]) for lab in labels}
    base = np.column_stack([arrays[lab] for lab in labels])
    gram = base @ base.T
    gram_inv = np.linalg.inv(gram)
    proper, improper = [], []
    for perm in permutations(labels):
        trial = np.column_stack([arrays[lab] for lab in perm])
        rot = trial @ base.T @ gram_inv
        if not np.allclose(rot @ base, trial, atol=1e-6):
            continue
        if not np.allclose(rot.T @ rot, np.eye(3), atol=1e-6):
            continue
        if np.linalg.det(rot) > 0:
            proper.append(perm)
        else:
            improper.append(perm)

    return labels, proper, improper


def _build_eq_ors(labels, perms):
    """Convert symmetry permutations into EqOrs table format."""
    return {
        lab: [perm[idx] for perm in perms]
        for idx, lab in enumerate(labels)
    }


def _build_syms(labels, proper_perms):
    """Generate one representative of each arrangement orbit."""
    rotations = [dict(zip(labels, perm)) for perm in proper_perms]
    seen = set()
    reps = []
    for arrangement in permutations(labels):
        if arrangement in seen:
            continue
        reps.append(arrangement)
        for rotation in rotations:
            rotation_inv = {val: key for key, val in rotation.items()}
            rotated = tuple(arrangement[rotation_inv[lab] - 1] for lab in labels)
            seen.add(rotated)

    return {idx + 1: list(rep) for idx, rep in enumerate(reps)}


def _volume(point1, point2, point3):
    """Signed tetrahedral volume with the central atom placed at the origin."""
    p1 = _point_to_array(point1)
    p2 = _point_to_array(point2)
    p3 = _point_to_array(point3)
    return float(np.dot(p1, np.cross(p2, p3)) / 6.0)


def _register_geometry(name, geom, posvs):
    """Populate all parameter tables for a geometry from idealized coordinates."""
    labels, proper, improper = _get_symmetry_permutations(geom)
    params.Syms[name] = _build_syms(labels, proper)
    params.Geoms[name] = geom
    params.Bounds[name] = _build_bounds(geom)
    params.PosVs[name] = [['CA', *face] for face in posvs]
    ideal_volume = sum(
        _volume(geom[a], geom[b], geom[c])
        for a, b, c in posvs
    )
    params.MinVs[name] = ideal_volume / 8.0
    params.EqOrs[name] = _build_eq_ors(labels, proper)
    if improper:
        params.EqOrs['enant' + name] = _build_eq_ors(labels, improper)
    params.Angles[name], params.Nears[name] = _build_angles_and_nears(geom)


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


# additional coordination polyhedra
_register_geometry(
    'TET',
    {'CA': Point3D( 0.0,  0.0,  0.0),
        1: Point3D( 1.154700538,  1.154700538,  1.154700538),
        2: Point3D( 1.154700538, -1.154700538, -1.154700538),
        3: Point3D(-1.154700538,  1.154700538, -1.154700538),
        4: Point3D(-1.154700538, -1.154700538,  1.154700538)},
    [(1, 2, 3),
     (1, 4, 2),
     (1, 3, 4),
     (2, 4, 3)]
)

_register_geometry(
    'SPY',
    {'CA': Point3D( 0.0,  0.0,  0.0),
        1: Point3D( 0.0,  0.0,  2.0),
        2: Point3D( 2.0,  0.0,  0.0),
        3: Point3D( 0.0,  2.0,  0.0),
        4: Point3D(-2.0,  0.0,  0.0),
        5: Point3D( 0.0, -2.0,  0.0)},
    [(1, 2, 3),
     (1, 3, 4),
     (1, 4, 5),
     (1, 5, 2)]
)

_register_geometry(
    'TBP',
    {'CA': Point3D( 0.0,  0.0,  0.0),
        1: Point3D( 0.0,  0.0,  2.0),
        2: Point3D( 2.0,  0.0,  0.0),
        3: Point3D(-1.0,  1.732050808,  0.0),
        4: Point3D(-1.0, -1.732050808,  0.0),
        5: Point3D( 0.0,  0.0, -2.0)},
    [(1, 2, 3),
     (1, 3, 4),
     (1, 4, 2),
     (5, 3, 2),
     (5, 4, 3),
     (5, 2, 4)]
)

params.Syms['SAN'] = {1: [1, 2]}
params.Geoms['SAN'] = {'CA': Point3D(0.0, 0.0, 0.0),
                          1: Point3D(0.0, 0.0, 2.0),
                          2: Point3D(0.0, 0.0,-2.0)}
params.Bounds['SAN'] = _build_bounds(params.Geoms['SAN'])
params.PosVs['SAN'] = []
params.MinVs['SAN'] = -1.0
params.EqOrs['SAN'] = {1: [1, 2],
                       2: [2, 1]}
params.Nears['SAN'] = {1: [],
                       2: []}
params.Angles['SAN'] = {1: {2: 180.0},
                        2: {1: 180.0}}

# M-to-centroid distances (Å) for haptic ligands, keyed by hapticity (number of
# coordinating atoms). Used in _SetCentralAtomBonds when the donor atom is a
# centroid dummy (* with a map number). These are geometry- and metal-independent
# defaults; the actual distance varies with metal and oxidation state.
#   η²: alkene/alkyne (Pd/Pt ~2.05 Å), σ-H₂ (W ~1.80 Å); approximate midpoint
#   η⁵: Cp (Fe ~1.64 Å, Ru ~1.82 Å, W ~2.03 Å); approximate midpoint
params.HapticDist = {
    2: 1.85,   # η²  (σ-H₂, alkene, alkyne)
    3: 2.00,   # η³  (allyl, open Cp)
    4: 2.00,   # η⁴  (butadiene, 1,5-COD fragment)
    5: 1.90,   # η⁵  (Cp, indenyl)
    6: 1.75,   # η⁶  (benzene, arene)
}

# Centroid-to-haptic-atom distances (Å) used in _SetHapticConstraints.
#   η²: approx half the ligand bond (H–H ≈ 0.41 Å, C=C ≈ 0.69 Å; midpoint used)
#   η⁵: Cp ring radius  (C–C = 1.40 Å, regular pentagon: R = 0.851 × C–C ≈ 1.19 Å)
#   η⁶: arene ring radius (C–C = 1.40 Å, regular hexagon:  R = C–C = 1.40 Å)
params.HapticCentR = {
    2: 0.55,   # η²
    3: 1.20,   # η³
    4: 1.15,   # η⁴
    5: 1.21,   # η⁵
    6: 1.40,   # η⁶
}


