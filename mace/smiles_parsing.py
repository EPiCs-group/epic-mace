'''
Functions for SMILES parsing
'''

#%% Imports

import re

from rdkit import Chem


#%% Functions

def _ParseBond(text):
    '''
    Parses bond's text
    '''
    info = {'closures': [],
            'branch_open': 0,
            'branch_close': 0,
            'usual_bond': True}
    # usual bond
    if text and text[-1] == '.':
        info['usual_bond'] = False
        text = text[:-1]
    # closures and branches
    while text:
        # closures
        flag = False
        for regex in [r'[-/\\~:=#$]*(\d)', r'[-/\\~:=#$]*\%(\d\d)']:
            m = re.match(regex, text)
            if m:
                info['closures'].append(int(m.group(1)))
                text = text[m.span()[1]:]
                flag = True
                break
        if flag:
            continue
        # branch closing
        m = re.match(r'\)+', text)
        if m:
            info['branch_close'] = m.span()[1]
            text = text[m.span()[1]:]
            continue
        # branch opening
        m = re.match(r'\([-/\\~:=#$]*', text)
        if m:
            info['branch_open'] += 1
            text = text[m.span()[1]:]
            continue
        # last bond
        if re.fullmatch(r'[-/\\~:=#$]', text):
            text = ''
            continue
        raise ValueError('Bad bond string')
    
    if info['branch_open'] > 1:
        raise ValueError('Bad bond string')
    
    return info


def _GetCXSmilesBondIdxs(smiles):
    '''
    Returns bond indexes in ChemAxon style
    '''
    # split smiles to atoms and bonds
    boundaries = []
    i = 0 # SMILES letter index
    idx = 0 # atom index
    N = len(smiles)
    n_atoms = 0
    while i < N:
        if smiles[i] == '[':
            # bracket atom
            start = i
            while smiles[i] != ']':
                i += 1
            i += 1
            boundaries += [start, i]
            n_atoms += 1
            idx += 1
        elif smiles[i] in ('B','C','N','O','S','P','F','I','b','c','n','o','s','p'):
            # organic subset
            add = 2 if smiles[i:i+2] in ('Cl', 'Br') else 1
            boundaries += [i, i+add]
            i += add
            n_atoms += 1
            idx += 1
        else:
            # if that's bond info just continue
            i += 1
    boundaries += [len(smiles)]
    # analyze bonds
    bonds = []
    closures = {}
    branches = []
    for i in range(1, len(boundaries) - 1, 2):
        a_idx = int((i-1)/2)
        #print(smiles[boundaries[i]:boundaries[i+1]])
        info = _ParseBond(smiles[boundaries[i]:boundaries[i+1]])
        # add closures
        for idx in info['closures']:
            if idx in closures:
                a_idx0 = closures.pop(idx)
                bonds.append( (a_idx0, a_idx) )
            else:
                closures[idx] = a_idx
        # branches
        N = info['branch_close']
        a_idx0 = None
        while N:
            N -= 1
            a_idx0 = branches.pop()
        if info['branch_open']:
            branches.append(a_idx if a_idx0 is None else a_idx0)
        # last bond
        if info['usual_bond'] and a_idx != n_atoms - 1:
            bonds.append( (a_idx if a_idx0 is None else a_idx0, a_idx + 1) )
    
    return bonds


def MolFromCXSmiles(smiles, params = None):
    '''
    Transforms ChemAxon SMILES to mol
    '''
    # split to smiles and info
    ps = smiles.split('|')
    smiles = ps[0].strip()
    info = ps[1]
    # extract indexes of dative bonds
    idxs = re.search(r'C:(\d+\.\d+,*)+', info)
    if idxs:
        idxs = idxs.group(0)[2:].split(',')
        idxs = [[int(_) for _ in db.split('.')] for db in idxs]
    else:
        idxs = []
    # extract R groups
    Rs = {}
    match = re.search('\$(.*?)\$', info)
    if match:
        for idx, char in enumerate(match.group(1).split(';')):
            if char[:2] == '_R':
                Rs[idx] = int(char[2:])
    # get atom indexes
    bonds = _GetCXSmilesBondIdxs(smiles)
    DBs = []
    for a_idx, b_idx in idxs:
        i, j = bonds[b_idx]
        if i != a_idx:
            i, j = j, i
        DBs.append( (i, j) )
    # make mol
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    ps.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params = ps)
    # set R groups
    for idx, num in Rs.items():
        mol.GetAtomWithIdx(idx).SetIsotope(num)
    # make dative bonds
    Chem.SetBondStereoFromDirections(mol)
    CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
    CTs = [Chem.rdchem.BondStereo.STEREOCIS, Chem.rdchem.BondStereo.STEREOTRANS]
    for i, j in DBs:
        # find needed stereo changes
        inverse = []
        doubles = []
        for idx1, idx2 in [(i, j), (j, i)]:
            # atom properties
            a = mol.GetAtomWithIdx(idx1)
            # inversion of chiral ai and aj
            if a.GetChiralTag() in CHIs:
                ns = [_.GetOtherAtom(a).GetIdx() for _ in a.GetBonds()]
                if (len(a.GetNeighbors()) == 4 and not ns.index(idx2) % 2) or \
                   (len(a.GetNeighbors()) == 3 and ns.index(idx2) % 2):
                    a.SetChiralTag(CHIs[not CHIs.index(a.GetChiralTag())])
            # cis/trans double bonds
            bonds = [_ for _ in a.GetBonds() if _.GetOtherAtom(a).GetIdx() != idx2]
            for b in bonds:
                bs = b.GetStereo()
                if bs in CTs:
                    a1 = b.GetBeginAtom()
                    a2 = b.GetEndAtom()
                    n1_idx, n2_idx = b.GetStereoAtoms()
                    if n1_idx not in [_.GetIdx() for _ in a1.GetNeighbors()]:
                        n1_idx, n2_idx = n2_idx, n1_idx
                    doubles.append([a1.GetIdx(), a2.GetIdx(), n1_idx, n2_idx, bs])
        # make bond
        ed = Chem.EditableMol(mol)
        ed.RemoveBond(i, j)
        ed.AddBond(i, j, Chem.BondType.DATIVE)
        mol = ed.GetMol()
        # set stereo
        for idx in inverse:
            a = mol.GetAtomWithIdx(idx)
            a.SetChiralTag(CHIs[not CHIs.index(a.GetChiralTag())])
        # set cis/trans (can't be set before mol edit)
        for a1_idx, a2_idx, n1_idx, n2_idx, bs in doubles:
            bond = mol.GetBondBetweenAtoms(a1_idx, a2_idx)
            bond.SetStereoAtoms(n1_idx, n2_idx)
            bond.SetStereo(bs)
    Chem.SanitizeMol(mol)
    
    return mol


def MolFromSmiles(smiles):
    '''
    Convert SMILES to RDKit Mol object taking into account MACE coding rules
    '''
    # read smiles with explicit Hs
    if '|' in smiles:
        mol = MolFromCXSmiles(smiles)
    else:
        ps = Chem.SmilesParserParams()
        ps.removeHs = False
        mol = Chem.MolFromSmiles(smiles, params = ps)
    # add isotopic label to atom with non-zero atom map number
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum():
            atom.SetIsotope(atom.GetAtomMapNum())
    # # substitute isotopic dummies without dative bonds to Hs
    # for atom in mol.GetAtoms():
    #     if atom.GetAtomicNum() == 0 and atom.GetIsotope():
    #         if 'DATIVE' in [str(_.GetBondType()) for _ in atom.GetBonds()]:
    #             continue
    #         atom.SetAtomicNum(1)
    
    return Chem.RemoveHs(mol)


def MolToSmiles(mol):
    '''
    Anti-RDKit plug
    '''
    
    return Chem.MolToSmiles(mol)


