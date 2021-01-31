'''
Supporting functions for the Complex object
'''

#%% Imports

from copy import deepcopy

from rdkit import Chem


#%% Functions

def _CalcTHVolume(conf, idxs):
    '''
    Calculates TH volume with given Point3D objects
    '''
    ps = [conf.GetAtomPosition(idx) for idx in idxs]
    v1 = [ps[1].x-ps[0].x, ps[1].y-ps[0].y, ps[1].z-ps[0].z]
    v2 = [ps[2].x-ps[0].x, ps[2].y-ps[0].y, ps[2].z-ps[0].z]
    v3 = [ps[3].x-ps[0].x, ps[3].y-ps[0].y, ps[3].z-ps[0].z]
    prod = [v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]]
    
    return sum([x*y for x, y in zip(prod, v3)])/6


def _RemoveRs(mol):
    '''
    Transforms R atoms into implicit form
    '''
    mol = deepcopy(mol)
    for a in mol.GetAtoms():
        if a.GetSymbol() in ('*', 'H') and len(a.GetNeighbors()) == 1 and \
           a.GetIsotope() and not a.GetAtomMapNum():
            a.SetAtomicNum(1)
            a.SetIsotope(0)
    
    return Chem.RemoveHs(mol)


