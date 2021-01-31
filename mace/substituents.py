'''
Functions for TH-symmetry-preserved substituents introduction
'''

#%% Imports

from rdkit import Chem


#%% Functions

def _BondByDummies(mol, idx1, idx2):
    '''
    Creates bond A-B from fragments A-* and B-*
    '''
    # prepare mol and atoms
    Chem.SetBondStereoFromDirections(mol)
    a1 = mol.GetAtomWithIdx(idx1).GetNeighbors()[0]
    a2 = mol.GetAtomWithIdx(idx2).GetNeighbors()[0]
    # chiral centers
    CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
    for a, idx in [(a1, idx1), (a2, idx2)]:
        if a.GetChiralTag() in CHIs:
            ns = [_.GetOtherAtom(a).GetIdx() for _ in a.GetBonds()]
            if (len(a.GetNeighbors()) == 4 and not ns.index(idx) % 2) or \
               (len(a.GetNeighbors()) == 3 and ns.index(idx) % 2):
                a.SetChiralTag(CHIs[not CHIs.index(a.GetChiralTag())])
    # E/Z double bonds
    CT = [Chem.rdchem.BondStereo.STEREOCIS, Chem.rdchem.BondStereo.STEREOTRANS]
    bonds = []
    for a, idx, idx_n in [(a1, idx1, a2.GetIdx()), (a2, idx2, a1.GetIdx())]:
        bs = [b for b in a.GetBonds() if b.GetOtherAtom(a).GetIdx() != idx]
        for b in bs:
            if b.GetStereo() in CT:
                i1, i2 = b.GetStereoAtoms()
                if i1 == idx:
                    bonds.append( (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), idx_n, i2, b.GetStereo()) )
                elif i2 == idx1:
                    bonds.append( (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), i1, idx_n, b.GetStereo()) )
    # edit mol
    ed = Chem.EditableMol(mol)
    ed.RemoveBond(a1.GetIdx(), idx1)
    ed.RemoveBond(a2.GetIdx(), idx2)
    ed.AddBond(a1.GetIdx(), a2.GetIdx(), Chem.BondType.SINGLE)
    mol = ed.GetMol()
    # edit bond stereo
    for a1_idx, a2_idx, i, j, bs in bonds:
        bond = mol.GetBondBetweenAtoms(a1_idx, a2_idx)
        bond.SetStereoAtoms(i, j)
        bond.SetStereo(bs)
    Chem.SanitizeMol(mol)
    
    return mol


def _CheckSubs(Rs):
    '''
    Checks substituents format
    '''
    for name, mol in Rs.items():
        if name[0] != 'R' or not name[1:].isdigit():
            raise ValueError(f'Bad substituent\'s name: {name}')
        dummies = [_ for _ in mol.GetAtoms() if _.GetSymbol() == '*']
        smiles = Chem.MolToSmiles(mol)
        if len(dummies) != 1:
            raise ValueError(f'Bad {name} substituent: {smiles}\n Substituent must contain exactly one dummy atom')
        if len(dummies[0].GetNeighbors()) != 1:
            raise ValueError(f'Bad {name} substituent\'s SMILES: {smiles}\n Substituent\'s dummy must be bonded to exactly one atom by single bond')
        if str(dummies[0].GetBonds()[0].GetBondType()) != 'SINGLE':
            raise ValueError(f'Bad {name} substituent\'s SMILES: {smiles}\n Substituent\'s dummy must be bonded to exactly one atom by single bond')
        dummies[0].SetIsotope(int(name[1:]))
    
    return None


def AddSubsToMol(mol, Rs):
    '''
    Updates molecule by adding substituents
    '''
    # check Rs
    _CheckSubs(Rs)
    # get number and type of Rs
    needed_Rs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in (0, 1) and atom.GetIsotope() and not atom.GetAtomMapNum():
            if 'DATIVE' not in [str(_.GetBondType()) for _ in atom.GetBonds()]:
                needed_Rs.append( (atom.GetIdx(), f'R{atom.GetIsotope()}') )
    # check Rs existence
    absent = set([_[1] for _ in needed_Rs]).difference(set(Rs.keys()))
    if absent:
        addend = 's' if len(absent) > 1 else ''
        absent = ','.join(absent)
        raise ValueError(f'No {absent} substituent{addend} found in list of Rs')
    # add Rs
    N = mol.GetNumAtoms()
    drop = []
    for idx1, R in needed_Rs:
        mol = Chem.CombineMols(mol, Rs[R])
        # find dummie's index in sub
        idx2 = None
        for i in range(N, mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() == '*' and atom.GetIsotope() == int(R[1:]):
                idx2 = atom.GetIdx()
                break
        # create bond
        mol = _BondByDummies(mol, idx1, idx2)
        mol.GetAtomWithIdx(idx1).SetIsotope(0)
        mol.GetAtomWithIdx(idx2).SetIsotope(0)
        drop += [idx1, idx2]
    # remove dummies
    ed = Chem.EditableMol(mol)
    for idx in sorted(drop, reverse = True):
        ed.RemoveAtom(idx)
    mol = ed.GetMol()
    Chem.SanitizeMol(mol)
    
    return mol


