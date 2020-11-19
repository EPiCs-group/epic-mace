'''
Chemoinformatic support for mononuclear complexes. Generates possible
stereoisomers and 3D coordinates. Supports square planar and octahedral
geometries
'''

#%% Must DOs:

# TODO: isotopes to maps for coord sphere
# TODO: add charge to properties
# TODO: AddBondedLigand: add stereo_dummy flag
# TODO: error texts
# TODO: functions and module descriptions
# TODO: full testing
# TODO: check double bond stereo in AddConstrainedConformer


#%% Optional improvements

# TODO: Na, Ca, Al, etc as CA - fix RDKit problem with dative bonds
# TODO: dummy as CA - substitute for any atom before embedding
# TODO: MolFromCXSmiles: @/@@ support - achiral/chiral ligands (do we need it?)
# TODO: enantiomers: set of unique structures must have same stereo for CA
# TODO: GetEnantiomer() method
# TODO: support of other geometries (tri-/pentagonal pyramid, etc)
# TODO: improved search for stereo centers
# TODO: Open SMILES support


#%% Imports

import os, re, json
from copy import deepcopy
from itertools import product, combinations

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point3D



#%% CXSmiles parsing

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
    # substitute atom mapping by isotopes
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum():
            atom.SetIsotope(atom.GetAtomMapNum())
            atom.SetAtomMapNum(0)
    # # substitute isotopic dummies without dative bonds to Hs
    # for atom in mol.GetAtoms():
    #     if atom.GetAtomicNum() == 0 and atom.GetIsotope():
    #         if 'DATIVE' in [str(_.GetBondType()) for _ in atom.GetBonds()]:
    #             continue
    #         atom.SetAtomicNum(1)
    
    return Chem.RemoveHs(mol)



#%% Substituents

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


def _SubsFromSmiles(RsSmiles):
    '''
    Transforms SMILES of substituents to RDKit Mol objects and checks them
    RsSmiles is dictionary, keys are 'R1', 'R2', etc, and values are
    substituents SMILES
    '''
    Rs = {}
    for name, smiles in RsSmiles.items():
        if name[0] != 'R' or not name[1:].isdigit():
            raise ValueError(f'Bad substituent\'s name: {name}')
        mol = MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f'Bad {name} substituent\'s SMILES: {smiles}')
        dummies = [_ for _ in mol.GetAtoms() if _.GetSymbol() == '*']
        if len(dummies) != 1:
            raise ValueError(f'Bad {name} substituent\'s SMILES: {smiles}\n Substituent must contain exactly one dummy atom')
        if len(dummies[0].GetNeighbors()) != 1:
            raise ValueError(f'Bad {name} substituent\'s SMILES: {smiles}\n Substituent\'s dummy must be bonded to exactly one atom by single bond')
        if str(dummies[0].GetBonds()[0].GetBondType()) != 'SINGLE':
            raise ValueError(f'Bad {name} substituent\'s SMILES: {smiles}\n Substituent\'s dummy must be bonded to exactly one atom by single bond')
        dummies[0].SetIsotope(int(name[1:]))
        Rs[name] = mol
    
    return Rs


def AddSubsToMol(mol, RsSmiles):
    '''
    Updates molecule by adding substituents
    '''
    # prepare Rs
    Rs = _SubsFromSmiles(RsSmiles)
    # get number and type of Rs
    needed_Rs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope():
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



#%% Complex initialization

def ComplexFromMol(mol, geom):
    '''
    Initializes Complex from RDKit Mol. Use it if mol->smiles transform is unwanted,
    e.g. if molecule contains *=DA->CA fragment with stereospecified double bond
    '''
    X = Complex('[*:1]->[*]', geom)
    if not mol or not Chem.MolToSmiles(mol):
        raise ValueError('Bad molecule: None or not RDKit Mol')
    X.mol = mol
    X._CheckMol()
    X._SetComparison()
    X._embedding_prepared = False
    X._ff_prepared = False
    
    return X


def ComplexFromLigands(ligands, CA, geom):
    '''
    Input:
      * ligands: the list of ligands' SMILES. Donor atoms are those ones
                 which have atom mapping
      * CA: the SMILES of central atom
    Output:
      * Complex object
    '''
    # combine ligands
    ligands = [MolFromSmiles(_) for _ in ligands]
    mol = ligands[0]
    for ligand in ligands[1:]:
        mol = Chem.CombineMols(mol, ligand)
    # Add Hs for easy stereo control
    mol = Chem.AddHs(mol)
    Chem.SetBondStereoFromDirections(mol)
    # add CA to Mol object
    CA = Chem.MolFromSmiles(CA)
    if not CA:
        raise ValueError('Error in SMILES of central atom')
    if CA.GetNumAtoms() != 1:
        raise ValueError('CA must be SMILES of one atom')
    CA = CA.GetAtomWithIdx(0)
    ed = Chem.EditableMol(mol)
    idx_CA = ed.AddAtom(CA)
    mol = ed.GetMol()
    # find dummies and make flips if needed
    DAs = [atom for atom in mol.GetAtoms() if atom.GetIsotope()]
    CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
    CTs = [Chem.rdchem.BondStereo.STEREOCIS, Chem.rdchem.BondStereo.STEREOTRANS]
    dummies = {}
    doubles = []
    for DA in DAs:
        idx_DA = DA.GetIdx()
        # is dummy bonded
        idx_dummy = None
        for i, b in enumerate(DA.GetBonds()):
            atom = b.GetOtherAtom(DA)
            if atom.GetSymbol() == '*' and len(atom.GetNeighbors()) == 1:
                idx_dummy = atom.GetIdx()
                break
        if idx_dummy is None:
            continue
        dummies[idx_DA] = idx_dummy
        # chirality
        if DA.GetChiralTag() in CHIs:
            if (len(DA.GetNeighbors()) == 4 and not i % 2) or \
               (len(DA.GetNeighbors()) == 3 and i % 2):
                DA.SetChiralTag(CHIs[not CHIs.index(DA.GetChiralTag())])
        # double bonds
        bs = [b for b in DA.GetBonds() if b.GetOtherAtomIdx(idx_DA) != idx_dummy]
        for b in bs:
            if b.GetStereo() not in CTs:
                continue
            idx1, idx2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            i1, i2 = b.GetStereoAtoms()
            if idx_dummy not in (i1, i2):
                continue
            if idx1 == idx_DA:
                doubles.append( (idx1, idx2, idx_CA, i2, b.GetStereo()) )
            else:
                doubles.append( (idx1, idx2, i1, idx_CA, b.GetStereo()) )
    # create dative bonds and remove bonds with dummies
    ed = Chem.EditableMol(mol)
    for DA in DAs:
        idx_DA = DA.GetIdx()
        if idx_DA in dummies:
            ed.RemoveBond(idx_DA, dummies[idx_DA])
        ed.AddBond(DA.GetIdx(), idx_CA, Chem.BondType.DATIVE)
    # reset double bonds stereo
    if doubles:
        mol = ed.GetMol()
        for a1_idx, a2_idx, i, j, bs in doubles:
            bond = mol.GetBondBetweenAtoms(a1_idx, a2_idx)
            bond.SetStereoAtoms(i, j)
            bond.SetStereo(bs)
        Chem.SanitizeMol(mol)
        ed = Chem.EditableMol(mol)
    # remove dummies
    for idx in sorted(dummies.values(), reverse = True):
        ed.RemoveAtom(idx)
    # final setting
    mol = ed.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    
    return ComplexFromMol(mol, geom)



#%% Read from file

def _ReadXYZ(path):
    '''
    Reads XYZ file of Complex and returns imitation of Complex object
    '''
    # read file
    if not os.path.isfile(path):
        raise ValueError('Bad XYZ path: file does not exist')
    with open(path, 'r') as inpf:
        text = [_.strip() for _ in inpf.readlines()]
    # first check and split to blocks
    if not text[0].isdigit():
        raise ValueError('Bad XYZ file: the first line must contain the number of atoms')
    N = int(text[0])
    if len(text) < N + 2:
        raise ValueError('Bad XYZ file: number of atoms in the coordinates block is less than specified in the first line')
    chunks = [text[i-N-1:i+1] for i in range(N + 1, len(text), N + 2)]
    # check format
    infos = []
    for i, chunk in enumerate(chunks):
        # first line format
        if not chunk[0].isdigit():
            raise ValueError('Bad XYZ file: the first line of a text chunk must contain the number of atoms')
        # second line format
        try:
            info = json.loads(chunk[1])
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            raise ValueError('Bad XYZ file: second line must be a JSON')
        if 'geom' not in info:
            raise ValueError('Bad XYZ file: second line must contain geometry of the central atom')
        if 'smiles' not in info or 'smiles3D' not in info or 'smiles3Dx' not in info:
            raise ValueError('Bad XYZ file: second line must contain SMILES strings of the complex')
        if 'dummies' not in info:
            raise ValueError('Bad XYZ file: second line must contain coordinates of dummies-helpers')
        if len(info['dummies']) % 3:
            raise ValueError('Bad XYZ file: bad number of dummies-helpers\' coordinates')
        if 'E' not in info:
            raise ValueError('Bad XYZ file: second line must contain MM energy of the complex')
        else:
            try:
                float(info['E'])
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                raise ValueError('Bad XYZ file: MM energy must be float')
        info['conf'] = i
        if 'rms' not in info:
            info['rms'] = -1.0
        # parse coordinates
        try:
            atoms = [_.split()[0] for _ in chunk[2:2+N]]
            coords = [Point3D(*[float(_) for _ in line.split()[1:4]]) for line in chunk[2:2+N]]
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            raise ValueError('Bad XYZ file: cannot read coordinates block')
        info['atoms'] = atoms
        info['coords'] = coords
        infos.append(info)
    # check that smiles and atoms is equal everythere
    info0 = infos[0]
    for info1 in infos[1:]:
        if info0['geom'] != info1['geom']:
            raise ValueError('Bad XYZ file: geometry must be the same for all structures')
        if info0['smiles'] != info1['smiles'] or info0['smiles3D'] != info1['smiles3D'] or info0['smiles3Dx'] != info1['smiles3Dx']:
            raise ValueError('Bad XYZ file: SMILES must be the same for all structures')
        if info0['atoms'] != info1['atoms']:
            raise ValueError('Bad XYZ file: atoms in coordinates block must be the same for all structures')
        if len(info0['dummies']) != len(info1['dummies']):
            raise ValueError('Bad XYZ file: number of dummies must be equivalent for all structures')
    # check molecules
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    mol = Chem.MolFromSmiles(info0['smiles'], params = ps)
    mol3D = Chem.MolFromSmiles(info0['smiles3D'], params = ps)
    mol3Dx = Chem.MolFromSmiles(info0['smiles3Dx'], params = ps)
    if not mol or not mol3D or not mol3Dx:
        raise ValueError('Bad XYZ file: not readable SMILES')
    if mol3D.GetNumAtoms() != N:
        raise ValueError('Bad XYZ file: number of atoms in molecule and coordinates block do not match')
    if mol3Dx.GetNumAtoms() - mol3D.GetNumAtoms() != int(len(info0['dummies'])/3):
        raise ValueError('Bad XYZ file: bad number of atoms in mol3Dx')
    # make mol and check mapping
    maps = []
    for atom in mol.GetAtoms():
        maps.append(atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    mol = Chem.RenumberAtoms(mol, tuple(zip(*sorted([(j, i) for i, j in enumerate(maps)])))[1])
    flag = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != atoms[atom.GetIdx()]:
            flag = True
            break
    if flag:
        raise ValueError('Bad XYZ file: atom numbering in SMILES and coordinates block do not match')
    # make mol3D and check mapping
    maps3D = []
    for atom in mol3D.GetAtoms():
        maps3D.append(atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    mol3D = Chem.RenumberAtoms(mol3D, tuple(zip(*sorted([(j, i) for i, j in enumerate(maps3D)])))[1])
    flag = False
    for atom in mol3D.GetAtoms():
        if atom.GetSymbol() != atoms[atom.GetIdx()]:
            flag = True
            break
    if flag:
        raise ValueError('Bad XYZ file: atom numbering in SMILES3D and coordinates block do not match')
    # make mol3Dx and check mapping
    maps3Dx = []
    for atom in mol3Dx.GetAtoms():
        maps3Dx.append(atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    mol3Dx = Chem.RenumberAtoms(mol3Dx, tuple(zip(*sorted([(j, i) for i, j in enumerate(maps3Dx)])))[1])
    flag = False
    for atom in mol3Dx.GetAtoms():
        if atom.GetIdx() < mol3D.GetNumAtoms():
            if atom.GetSymbol() != atoms[atom.GetIdx()]:
                flag = True
                break
        else:
            if atom.GetSymbol() != '*':
                flag = True
                break
    if flag:
        raise ValueError('Bad XYZ file: atom numbering in SMILES3Dx and coordinates block do not match')
    # final checks
    if Chem.MolToSmiles(mol, canonical = False) != Chem.MolToSmiles(Chem.RemoveHs(mol3D), canonical = False):
        raise ValueError('Bad XYZ file: SMILES and SMILES3D do not match')
    
    return mol, mol3D, mol3Dx, infos


def ComplexFromXYZFile(path):
    '''
    Reads XYZ file of Complex and returns imitation of Complex object
    '''
    mol, mol3D, mol3Dx, infos = _ReadXYZ(path)
    # make Complex
    try:
        X = ComplexFromMol(mol, infos[0]['geom'])
    except ValueError as e:
        raise ValueError('Bad XYZ file: unsuccessful SMILES to Complex conversion: ' + e.args[0])
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        raise ValueError('Bad XYZ file: unknown error during SMILES to Complex conversion')
    # check stereo info
    if X.err_init:
        raise ValueError('Bad XYZ file: SMILES does not contain enough stereo information')
    # preset embedding
    X.mol3D = mol3D
    X._SetEmbedding()
    X.mol3Dx = mol3Dx
    # add conformers
    for info in infos:
        # add coords to mol
        conf = AllChem.Conformer()
        for atom in X.mol.GetAtoms():
            conf.SetAtomPosition(atom.GetIdx(), info['coords'][atom.GetIdx()])
        conf.SetDoubleProp('E', info['E'])
        conf.SetDoubleProp('EmbedRMS', info['rms'])
        X.mol.AddConformer(conf, assignId = True)
        # add mol3D and its coords
        conf3D = AllChem.Conformer()
        for atom in X.mol3D.GetAtoms():
            conf3D.SetAtomPosition(atom.GetIdx(), info['coords'][atom.GetIdx()])
        conf3D.SetDoubleProp('E', info['E'])
        conf3D.SetDoubleProp('EmbedRMS', info['rms'])
        X.mol3D.AddConformer(conf3D, assignId = True)
        # add mol3Dx
        conf3Dx = AllChem.Conformer()
        for atom in X.mol3Dx.GetAtoms():
            if atom.GetIdx() < X.mol3D.GetNumAtoms():
                conf3Dx.SetAtomPosition(atom.GetIdx(), info['coords'][atom.GetIdx()])
            else:
                i = atom.GetIdx() - X.mol3D.GetNumAtoms()
                p = Point3D(*info['dummies'][3*i:3*(i+1)])
                conf3Dx.SetAtomPosition(atom.GetIdx(), p)
        conf3Dx.SetDoubleProp('E', info['E'])
        conf3Dx.SetDoubleProp('EmbedRMS', info['rms'])
        X.mol3Dx.AddConformer(conf3Dx, assignId = True)
    
    return X



#%% Support functions

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



#%% Complex object

class Complex():
    '''
    Wrapper around RDKit Mol object for mononuclear complexes
    '''
    
    ##########################
    # Force Field parameters #
    ##########################
    
    # force field parameters
    _FFParams = {'X*'     :    2.0,
                 'tether' :  100.0,
                 'kXL'    :  500.0,
                 'kX*'    : 1000.0,
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
    _Rcov = [0.23,0.23,1.5,1.28,0.96,0.83,0.68,0.68,0.68,0.64,1.5,1.66,1.41,1.21,1.2,1.05,1.02,
             0.99,1.51,2.03,1.76,1.7,1.6,1.53,1.39,1.61,1.52,1.26,1.24,1.32,1.22,1.22,1.17,
             1.21,1.22,1.21,1.5,2.2,1.95,1.9,1.75,1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.54,
             1.42,1.39,1.39,1.47,1.4,1.5,2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,
             1.94,1.92,1.92,1.89,1.9,1.87,1.87,1.75,1.7,1.62,1.51,1.44,1.41,1.36,1.36,1.32,
             1.45,1.46,1.48,1.4,1.21,1.5,2.6,2.21,2.15,2.06,2,1.96,1.9,1.87,1.8,1.69,1.54,
             1.83,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]
    _Rcov = {i: rcov for i, rcov in enumerate(_Rcov)}
    
    
    ############################
    # Basic Complex Parameters #
    ############################
    
    # SMILES symmetry codes for octahedral geometry @OH1-@OH30
    # see OpenSMILES specification for the details
    _Syms = {'OH': { 1: [1,2,3,4,5,6],  2: [1,2,5,4,3,6],  3: [1,2,3,4,6,5],
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
    _Geoms = {'OH': {'CA': Point3D( 0.0, 0.0, 0.0),
                        1: Point3D( 0.0, 0.0, 2.0),
                        2: Point3D( 2.0, 0.0, 0.0),
                        3: Point3D( 0.0, 2.0, 0.0),
                        4: Point3D(-2.0, 0.0, 0.0),
                        5: Point3D( 0.0,-2.0, 0.0),
                        6: Point3D( 0.0, 0.0,-2.0)},
              'SP': {'CA': Point3D( 0.0, 0.0, 0.0),
                        1: Point3D( 2.0, 0.0, 0.0),
                        2: Point3D( 0.0, 2.0, 0.0),
                        3: Point3D(-2.0, 0.0, 0.0),
                        4: Point3D( 0.0,-2.0, 0.0),
                     'X1': Point3D( 0.0, 0.0, 2.0),
                     'X2': Point3D( 0.0, 0.0,-2.0)}}
    
    # lists of points corresponding to positive tetrahedra volumes
    _PosVs = {'OH': [['CA',1,2,3],
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
    _MinVs = {'OH': 4/3, 'SP': 2/3}
    
    # equivalent orientations of complexes
    _EqOrs = {'OH': {1: [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6],
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
    _EqOrs['enantOH'] = deepcopy(_EqOrs['OH'])
    _EqOrs['enantOH'][1] = deepcopy(_EqOrs['OH'][6])
    _EqOrs['enantOH'][6] = deepcopy(_EqOrs['OH'][1])
    
    # neigboring ligands positions
    _Nears = {'OH': {1: [2,3,4,5],
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
    _Angles = {'OH': {1: {2:  90.0, 3:  90.0, 4:  90.0, 5:  90.0, 6: 180.0},
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
    
    
    #############################################
    # Initialization and graph-based comparison #
    #############################################
    
    def _CheckMol(self):
        # find and check dative bonds
        info = [(b.GetIdx(), b.GetBeginAtom(), b.GetEndAtom()) for b in self.mol.GetBonds() \
                if str(b.GetBondType()) == 'DATIVE']
        info.sort(key = lambda x: x[0])
        if len(info) == 0:
            raise ValueError('Bad SMILES: no dative bonds')
        # check number of central atoms
        CAs = set([_[2].GetIdx() for _ in info])
        if len(CAs) > 1:
            raise ValueError('Bad SMILES: there are several acceptors of dative bonds (central atoms)')
        # check CA's bonds
        self._idx_CA = list(CAs)[0]
        CA = self.mol.GetAtomWithIdx(self._idx_CA)
        if len(CA.GetBonds()) > len(info):
            raise ValueError('Bad SMILES: some bonds with central atom are not dative')
        if CA.GetNumImplicitHs() > 0:
            raise ValueError('Bad SMILES: all hydrogens (hydrides) bonded to central atom must be encoded explicitly with isotopic label')
        # check donor atoms labelling
        self._DAs = {_[1].GetIdx(): _[1].GetIsotope() for _ in info}
        labs = list(self._DAs.values())
        if 0 in labs:
            self.err_init = 'Bad SMILES: some donor atoms don\'t have an isotopic label'
            #self._PrintErrorInit()
            return
        elif len(set(labs)) != len(labs):
            self.err_init = 'Bad SMILES: isotopic labels are not unique'
            #self._PrintErrorInit()
            return
        elif max(labs) > max([_ for _ in self._Geoms[self._geom] if str(_).isdigit()]):
            self.err_init = 'Bad SMILES: maximal isotopic label exceeds number of ligands'
            #self._PrintErrorInit()
            return
    
    
    def _SetComparison(self):
        '''
        Prepares set of all possible smiles to compare with other complexes
        '''
        mol_norm = deepcopy(self.mol)
        # modify X<-[C-]=[N+] fragments to X<-[C]-[N] (symmetry problems for carbenes) # GetSubstructMatches does not work properly with kekulization
        Chem.Kekulize(mol_norm, clearAromaticFlags = True)
        for i, j, k in mol_norm.GetSubstructMatches(Chem.MolFromSmarts('[*]<-[C-]=[N+]')):
            mol_norm.GetAtomWithIdx(j).SetFormalCharge(0)
            mol_norm.GetAtomWithIdx(k).SetFormalCharge(0)
            mol_norm.GetBondBetweenAtoms(j, k).SetBondType(Chem.rdchem.BondType.SINGLE)
        mol_norm = Chem.MolFromSmiles(Chem.MolToSmiles(mol_norm, canonical = False))
        # generate all ligand orientations
        mol = deepcopy(mol_norm)
        Chem.SetBondStereoFromDirections(mol)
        # remove dative bonds # HINT: remove after RDKit fix
        drop = []
        CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
        for b in mol.GetBonds():
            if str(b.GetBondType()) == 'DATIVE':
                if b.GetBeginAtom().GetChiralTag() in CHIs:
                    continue
                drop.append( (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) )
        # create bond with dummy
        ed = Chem.EditableMol(mol)
        for i, j in drop:
            ed.RemoveBond(i, j)
        mol = ed.GetMol()
        Chem.SanitizeMol(mol)
        # invert mol
        mol_inv = deepcopy(mol)
        CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
        for atom in mol_inv.GetAtoms():
            # revert all chiral centers
            tag = atom.GetChiralTag()
            if tag in CHIs:
                atom.SetChiralTag(CHIs[not CHIs.index(tag)])
        # generate structure descriptor
        EqOrs = self._EqOrs[self._geom]
        if 'enant' + self._geom in self._EqOrs:
            EqOrsInv = self._EqOrs['enant' + self._geom]
        else:
            EqOrsInv = self._EqOrs[self._geom]
        _ID = []
        _eID = []
        _DAs = {_.GetIdx(): _.GetIsotope() for _ in mol.GetAtoms() if _.GetIsotope()}
        _DAs_inv = {_.GetIdx(): _.GetIsotope() for _ in mol_inv.GetAtoms() if _.GetIsotope()}
        for i in range(len(EqOrs[1])):
            for idx, num in _DAs.items():
                mol.GetAtomWithIdx(idx).SetIsotope(EqOrs[num][i])
            for idx, num in _DAs_inv.items():
                mol_inv.GetAtomWithIdx(idx).SetIsotope(EqOrsInv[num][i])
            # add resonance structures to mols
            for m in Chem.ResonanceMolSupplier(mol):
                _ID.append( Chem.CanonSmiles(Chem.MolToSmiles(m)) )
            # same for inv mol
            for m in Chem.ResonanceMolSupplier(mol_inv):
                _eID.append( Chem.CanonSmiles(Chem.MolToSmiles(m)) )
        self._ID = set(_ID)
        self._eID = set(_eID)
    
    
    def __init__(self, smiles, geom):
        '''
        Generates Complex object from SMILES with stereo info encoded as
        isotopic labels of donor atoms
        '''
        self._smiles_init = smiles
        self._geom = geom
        self.err_init = None
        # check geom
        if self._geom not in self._Geoms:
            raise ValueError(f'Unknown geometry type: {self._geom}')
        # check smiles
        self.mol = MolFromSmiles(self._smiles_init)
        if not self.mol:
            raise ValueError('Bad SMILES: not readable')
        # check mol
        self._CheckMol()
        # other stuff
        self._SetComparison()
        self.mol3D = Chem.AddHs(self.mol)
        self._embedding_prepared = False
        self._ff_prepared = False
    
    
    def _PrintErrorInit(self):
        '''
        Prints warning if complex does not have enough stereo info
        '''
        message = [self.err_init,
                   '',
                   'The initial SMILES contains insufficient or erroneous info',
                   'on the positions of the ligands around the central atom',
                   'encoded with isotopic labels.',
                   'To use 3D generation and other features, generate',
                   'possible stereomers by passing this Complex object',
                   'to GetComplexStereomers function with "CA" or "all" regime.',
                   '', '']
        if self.err_init:
            print('\n'.join(message))
            return True
        
        return False
    
    
    def IsEqual(self, X):
        '''
        Compares two Complex objects for equivalence
        '''
        if self._PrintErrorInit():
            return None
        
        return bool(self._ID.intersection(X._ID))
    
    
    def IsEnantiomeric(self):
        '''
        Checks is the complex enantiomeric
        '''
        if self._PrintErrorInit():
            return None
        
        return not bool(self._ID.intersection(self._eID))
    
    
    def IsEnantiomer(self, X):
        '''
        Checks if two complexes are enantiomers
        '''
        if self._PrintErrorInit():
            return None
        
        return bool(self._ID.intersection(X._eID))
    
    
    #####################
    # Stereomers search #
    #####################
    
    def _FindNeighboringDAs(self, minTransCycle = None):
        '''
        Finds restrictions of multidentate ligands.
        Returns pairs of DAs which must be near each over
        '''
        # get DAs
        DAs = list(self._DAs.keys())
        # get rings
        rings = [list(r) for r in Chem.GetSymmSSSR(self.mol)]
        rings = [r for r in rings if self._idx_CA in r]
        # find restrictions
        restrictions = []
        for r in rings:
            if minTransCycle and len(r) >= minTransCycle:
                continue
            r = [idx for idx in r if idx in DAs]
            if len(r) == 2:
                restrictions.append(r)
        
        return restrictions
    
    
    def GetStereomers(self, regime = 'all', dropEnantiomers = True, minTransCycle = None):
        '''
        Generates all possible stereomers of a complex.
        Saves stereochemistry of existing centers.
        Three regimes are available:
            - "CA": only stereochemistry of central atom;
            - "ligands": only stereochemistry of ligands;
            - "all": both "CA" and "ligands" regimes.
        '''
        if regime not in ('CA', 'ligands', 'all'):
            raise ValueError('regime variable bad value')
        # set Mol object
        mol = deepcopy(self.mol)
        if regime == 'ligands':
            # check numbering
            if self.err_init:
                raise ValueError('Stereo info for the central atom is not specified correctly. Use "CA" or "all" regimes to fix that')
        else:
            # randomly set isotopic numbers
            DAs = list(self._DAs.keys())
            for num, idx in enumerate(DAs):
                mol.GetAtomWithIdx(idx).SetIsotope(num + 1)
        # generate needed stereomers
        if regime == 'CA':
            mols = [mol]
        else:
            idxs = [idx for idx, chi in Chem.FindMolChiralCenters(mol, includeUnassigned = True) if chi == '?']
            # generate all possible combinations of stereocentres
            mols = []
            for chis in product([Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW], repeat = len(idxs)):
                m = deepcopy(mol)
                for idx, chi in zip(idxs, chis):
                    m.GetAtomWithIdx(idx).SetChiralTag(chi)
                mols.append(m)
            # drop similar ones
            smiles = [Chem.MolToSmiles(m) for m in mols]
            drop = []
            for i in range(len(smiles)-1):
                if i in drop:
                    continue
                for j in range(i+1, len(smiles)):
                    if smiles[i] == smiles[j]:
                        drop.append(j)
            mols = [m for i, m in enumerate(mols) if i not in drop]
        # generate all CA isomers
        if regime == 'ligands':
            # transform mols to Complex objects and return them
            stereomers = []
            for m in mols:
                stereomers.append( Complex(Chem.MolToSmiles(m), self._geom) )
            return stereomers
        # find restrictions on DA positions
        pairs = self._FindNeighboringDAs(minTransCycle)
        # generate all possible CA orientations
        stereomers = []
        for m in mols:
            addend = []
            for idx_sym in sorted(list(self._Syms[self._geom].keys())):
                sym = self._Syms[self._geom][idx_sym]
                m1 = deepcopy(m)
                # set new isotopes
                info = {}
                for idx, a_idx in enumerate(DAs):
                    num = sym.index(idx + 1) + 1
                    m1.GetAtomWithIdx(a_idx).SetIsotope(num)
                    info[a_idx] = num
                # check neighboring DAs restriction
                drop = False
                for idx_a, idx_b in pairs:
                    if info[idx_b] not in self._Nears[self._geom][info[idx_a]]:
                        drop = True
                if drop:
                    continue
                addend.append(m1)
            addend = [Complex(Chem.MolToSmiles(m), self._geom) for m in addend]
            # filter uniques
            drop = []
            for i in range(len(addend)-1):
                if i in drop:
                    continue
                for j in range(i+1, len(addend)):
                    if addend[i].IsEqual(addend[j]):
                        drop.append(j)
            addend = [compl for i, compl in enumerate(addend) if i not in drop]
            # add to main
            stereomers += addend
        # final filtering
        drop = []
        for i in range(len(stereomers)-1):
            if i in drop:
                continue
            for j in range(i+1, len(stereomers)):
                if stereomers[i].IsEqual(stereomers[j]):
                    drop.append(j)
                elif dropEnantiomers and stereomers[i].IsEnantiomer(stereomers[j]):
                    drop.append(j)
        stereomers = [compl for i, compl in enumerate(stereomers) if i not in drop]
        
        return stereomers    
    
    
    
    #################
    # 3D Generation #
    #################
    
    def _SetEmbedding(self):
        '''
        Sets cordinates' constraints
        '''
        # find dummies-helpers
        add = [idx for idx in self._Geoms[self._geom] if 'X' in str(idx)]
        must = set([idx for idx in self._Geoms[self._geom] if str(idx).isdigit()])
        have = set([num for num in self._DAs.values()])
        add += list(must.difference(have))
        # basic coordMap
        self._coordMap = {self._idx_CA: self._Geoms[self._geom]['CA']}
        for idx, num in self._DAs.items():
            self._coordMap[idx] = self._Geoms[self._geom][num]
        # add dummies-helpers and their coords
        self._dummies = {}
        ed = Chem.EditableMol(self.mol3D)
        for num in add:
            idx = ed.AddAtom(Chem.Atom(0))
            ed.AddBond(idx, self._idx_CA, Chem.BondType.DATIVE)
            self._dummies[idx] = num
            self._coordMap[idx] = self._Geoms[self._geom][num]
        # final
        self.mol3Dx = ed.GetMol()
        Chem.SanitizeMol(self.mol3Dx)
        self._embedding_prepared = True
    
    
    def _SetCentralAtomAngles(self):
        '''
        Sets restrictions to L->X<-L angles
        '''
        DAs = list(self._DAs.items()) + list(self._dummies.items())
        for i, j in combinations(range(len(DAs)), r = 2):
            a1_idx, a1_num = DAs[i]
            a2_idx, a2_num = DAs[j]
            angle = self._Angles[self._geom][a1_num][a2_num]
            if a1_num in self._Nears[self._geom][a2_num]:
                k = self._FFParams['kZ-LXL']
            else:
                k = self._FFParams['kE-LXL']
            constraint = [a1_idx, self._idx_CA, a2_idx, False, angle, angle, k]
            self._angle_params.append(constraint)
            self._ff.UFFAddAngleConstraint(*constraint)
    
    
    def _SetCentralAtomBonds(self):
        '''
        Sets restrictions to X<-L bonds
        '''
        for idx in self._DAs:
            dist = self._Rcov[self.mol3Dx.GetAtomWithIdx(self._idx_CA).GetAtomicNum()] + \
                   self._Rcov[self.mol3Dx.GetAtomWithIdx(idx).GetAtomicNum()]
            constraint = [idx, self._idx_CA, False, dist, dist, self._FFParams['kXL']]
            self._bond_params.append(constraint)
            self._ff.UFFAddDistanceConstraint(*constraint)
        # dummies-helpers
        for idx, num in self._dummies.items():
            if 'X' in str(num):
                constraint = [idx, self._idx_CA, False, self._FFParams['X*'], self._FFParams['X*'], self._FFParams['kX*']]
                self._bond_params.append(constraint)
                self._ff.UFFAddDistanceConstraint(*constraint)
            else:
                dist = self._Rcov[self.mol3Dx.GetAtomWithIdx(self._idx_CA).GetAtomicNum()] + \
                       self._Rcov[self.mol3Dx.GetAtomWithIdx(idx).GetAtomicNum()]
                constraint = [idx, self._idx_CA, False, dist, dist, self._FFParams['kXL']]
                self._bond_params.append(constraint)
                self._ff.UFFAddDistanceConstraint(*constraint)
    
    
    def _SetDonorAtomsAngles(self):
        '''
        Sets restrictions to X<-L-A angles
        '''
        for DA in self._DAs:
            # get neighbors
            ns = self.mol3Dx.GetAtomWithIdx(DA).GetNeighbors()
            ns = [n.GetIdx() for n in ns if n.GetIdx() != self._idx_CA]
            if not ns or len(ns) > 3:
                continue
            # set angles
            if len(ns) == 1 and str(self.mol3Dx.GetAtomWithIdx(ns[0]).GetHybridization()) == 'SP2':
                angle = self._FFParams['XLO']
                k = self._FFParams['kXLO']
            else:
                angle = {1: self._FFParams['XLA'],
                         2: self._FFParams['XLA2'],
                         3: self._FFParams['XLA3']}[len(ns)]
                k = self._FFParams['kXLA']
            for n in ns:
                constraint = [self._idx_CA, DA, n, False, angle, angle, k]
                self._angle_params.append(constraint)
                self._ff.UFFAddAngleConstraint(*constraint)
    
    
    def _SetDonorAtomsParams(self):
        '''
        Sets DA parameters without CA
        '''
        for DA in self._DAs:
            # get neighbors
            ns = [_.GetIdx() for _ in self.mol3Dx.GetAtomWithIdx(DA).GetNeighbors()]
            ns = [_ for _ in ns if _ != self._idx_CA]
            if not ns:
                continue
            N = len(ns)
            # bonds
            for n in ns:
                d = self._Rcov[self.mol3Dx.GetAtomWithIdx(n).GetAtomicNum()] + \
                    self._Rcov[self.mol3Dx.GetAtomWithIdx(DA).GetAtomicNum()]
                constraint = [DA, n, False, d, d, 700.0]
                self._bond_params.append(constraint)
                self._ff.UFFAddDistanceConstraint(*constraint)
            # X<-L-A angles
            if N < 2 or N > 3:
                continue
            for n1, n2 in combinations(ns, r = 2):
                a = self._FFParams['XLA2'] if N == 2 else self._FFParams['XLA3']
                constraint = [n1, DA, n2, False, a, a, self._FFParams['kALA']]
                self._angle_params.append(constraint)
                self._ff.UFFAddAngleConstraint(*constraint)
            # L-A-B angles
            for n in ns:
                atom = self.mol3Dx.GetAtomWithIdx(n)
                n2s = [_.GetIdx() for _ in atom.GetNeighbors()]
                N = len(n2s)
                if N < 3 or N > 4:
                    continue
                for n21, n22 in combinations(n2s, r = 2):
                    if str(atom.GetHybridization()) == 'SP2':
                        a = self._FFParams['XLA2']
                        constraint = [n21, n, n22, False, a, a, self._FFParams['kALA']]
                        self._angle_params.append(constraint)
                        self._ff.UFFAddAngleConstraint(*constraint)
                    elif str(atom.GetHybridization()) == 'SP3':
                        a = self._FFParams['XLA3']
                        constraint = [n21, n, n22, False, a, a, self._FFParams['kALA']]
                        self._angle_params.append(constraint)
                        self._ff.UFFAddAngleConstraint(*constraint)
    
    
    def _SetForceField(self, confId):
        '''
        Sets force field for the geometry optimization
        '''
        self._ff = AllChem.UFFGetMoleculeForceField(self.mol3Dx, confId = confId)
        if self._ff_prepared:
            # restore them from saved params
            for constraint in self._angle_params:
                self._ff.UFFAddAngleConstraint(*constraint)
            for constraint in self._bond_params:
                self._ff.UFFAddDistanceConstraint(*constraint)
            return
        # set ff parameters
        self._angle_params = []
        self._bond_params = []
        self._SetCentralAtomAngles()
        self._SetCentralAtomBonds()
        self._SetDonorAtomsAngles()
        self._SetDonorAtomsParams()
        self._ff_prepared = True
    
    
    def _CheckStereoCA(self, confId):
        '''
        Checks that generated coordinates corresponds to initial chirality
        '''
        conf = self.mol3Dx.GetConformer(confId)
        DAs = {val: key for key, val in self._DAs.items()}
        DAs['CA'] = self._idx_CA
        for key, val in self._dummies.items():
            DAs[val] = key
        Vs = []
        for idxs in self._PosVs[self._geom]:
            Vs.append(_CalcTHVolume(conf, [DAs[idx] for idx in idxs]))
        
        return sum(Vs) > self._MinVs[self._geom] # sum(Vs) > 0
    
    
    def Optimize(self, confId = 0, maxIts = 1000):
        '''
        Optimizes Complex
        '''
        if self._PrintErrorInit():
            return None
        # optimization
        self._SetForceField(confId)
        self._ff.Initialize()
        flag = self._ff.Minimize(maxIts = maxIts)
        # energy
        E = self._ff.CalcEnergy()
        conf3Dx = self.mol3Dx.GetConformer(confId)
        conf3Dx.SetDoubleProp('E', E)
        # synchronize with mol3D
        conf3D = self.mol3D.GetConformer(confId)
        conf3D.SetDoubleProp('E', E)
        for atom in self.mol3D.GetAtoms():
            idx = atom.GetIdx()
            conf3D.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
        # synchronize with mol
        conf = self.mol.GetConformer(confId)
        conf.SetDoubleProp('E', E)
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            conf.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
        
        return flag
    
    
    def AddConformer(self, clearConfs = True, maxAttempts = 10):
        '''
        Generates complex conformer using constrained embedding
        '''
        if self._PrintErrorInit():
            return None
        flag = -1
        attempt = maxAttempts
        while flag == -1 and attempt > 0:
            attempt -= 1
            # embedding
            if not self._embedding_prepared:
                self._SetEmbedding()
            flag = AllChem.EmbedMolecule(self.mol3Dx, coordMap = self._coordMap,
                                         clearConfs = clearConfs,
                                         enforceChirality = True)
            if flag == -1:
                continue
            # optimization # HINT: do not use self.Optimize as we need to apply self._CheckStereoCA after
            self._SetForceField(flag)
            self._ff.Initialize()
            self._ff.Minimize(maxIts = 1000)
            # check chiral centers from 3D
            if not self._CheckStereoCA(flag):
                self.mol3Dx.RemoveConformer(flag)
                flag = -1
                continue
            # refresh other confs
            if clearConfs:
                self.mol.RemoveAllConformers()
                self.mol3D.RemoveAllConformers()
            # move CA to (0,0,0)
            conf3Dx = self.mol3Dx.GetConformer(flag)
            r = deepcopy(conf3Dx.GetAtomPosition(self._idx_CA))
            for i in range(conf3Dx.GetNumAtoms()):
                conf3Dx.SetAtomPosition(i, conf3Dx.GetAtomPosition(i) - r)
            # energy
            E = self._ff.CalcEnergy()
            conf3Dx.SetDoubleProp('E', E)
            conf3Dx.SetDoubleProp('EmbedRMS', -1)
            # synchronize with mol3D
            conf3D = AllChem.Conformer()
            conf3D.SetDoubleProp('E', E)
            conf3D.SetDoubleProp('EmbedRMS', -1)
            for atom in self.mol3D.GetAtoms():
                idx = atom.GetIdx()
                conf3D.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            conf3D.SetId(flag)
            self.mol3D.AddConformer(conf3D, assignId = True)
            # synchronize with mol
            conf = AllChem.Conformer()
            conf.SetDoubleProp('E', E)
            conf.SetDoubleProp('EmbedRMS', -1)
            for atom in self.mol.GetAtoms():
                idx = atom.GetIdx()
                conf.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            conf.SetId(flag)
            self.mol.AddConformer(conf, assignId = True)
        
        return flag
    
    
    def AddConstrainedConformer(self, core, confId = 0, ignoreHs = True,
                                      clearConfs = True, maxAttempts = 10):
        '''
        Constrained embedding using other complex geometry
        Core complex must be a substructure of complex and
        must contain CA and all DAs (excluding supportive dummies)
        '''
        if len(Chem.GetMolFrags(core.mol)) != 1:
            raise ValueError('Bad core: core must contain exactly one fragment')
        # make mol3Dx and mol3D
        if not self._embedding_prepared:
            self._SetEmbedding()
        # substructure check
        match = self.mol3Dx.GetSubstructMatch(core.mol)
        if not match:
            raise ValueError('Bad core: core is not a substructure of the complex')
        if self._idx_CA not in match:
            raise ValueError('Bad core: core must contain the central atom')
        # prepare core
        if ignoreHs:
            core_mol = deepcopy(core.mol)
        else:
            core_mol = deepcopy(core.mol3D)
            drop = [] # list of Hs to remove
            info = {} # info on number of implicit Hs
            for i, idx in enumerate(match):
                a1 = self.mol3Dx.GetAtomWithIdx(idx)
                a2 = core_mol.GetAtomWithIdx(i)
                n1 = len([_ for _ in a1.GetNeighbors() if _.GetSymbol() == 'H'])
                n2 = len([_ for _ in a2.GetNeighbors() if _.GetSymbol() == 'H'])
                if n1 == n2:
                    continue
                # find H idxs to drop
                drop += [_.GetIdx() for _ in a2.GetNeighbors() if _.GetSymbol() == 'H']
                # explicit atoms
                info[a2.GetIdx()] = (n2, a2.GetNumRadicalElectrons())
            # remove Hs
            ed = Chem.EditableMol(core_mol)
            for idx in sorted(drop, reverse = True):
                ed.RemoveAtom(idx)
            core_mol = ed.GetMol()
            # set implicit Hs and radicals
            for idx, (nh, nr) in info.items():
                core_mol.GetAtomWithIdx(idx).SetNumExplicitHs(nh)
                core_mol.GetAtomWithIdx(idx).SetNumRadicalElectrons(nr)
            Chem.SanitizeMol(core_mol)
            match = self.mol3Dx.GetSubstructMatch(core_mol)
            if not match:
                raise RuntimeError('Something went wrong: core is not substructure of the molecule (ignoreHs = False)')
        # prepare embed params
        coordMap = {}
        coreConf = core_mol.GetConformer(confId)
        for i, idxI in enumerate(match):
            coordMap[idxI] = coreConf.GetAtomPosition(i)
        # check if some DAs missed in new coordMap
        add = [idx for idx in self._coordMap if idx not in match]
        if add:
            # make dummy mol
            dummy = Chem.MolFromSmiles('.'.join(['[*]']*len(self._coordMap))) # CA
            dummyMap = {}
            conf = Chem.Conformer()
            for i, (idx, point) in enumerate(self._coordMap.items()):
                dummyMap[i] = idx
                conf.SetAtomPosition(i, point)
            dummy.AddConformer(conf)
            # orient dummy mol over core
            algMap = [(key, match.index(val)) for key, val in dummyMap.items() if val in match]
            AllChem.AlignMol(dummy, core_mol, atomMap = algMap, maxIters = 200)
            # renew coordmap
            dummyMap = {val: key for key, val in dummyMap.items()}
            for idx in add:
                coordMap[idx] = dummy.GetConformer().GetAtomPosition(dummyMap[idx])
        # embedding
        flag = -1
        attempt = maxAttempts
        while flag == -1 and attempt > 0:
            attempt -= 1
            flag = AllChem.EmbedMolecule(self.mol3Dx, coordMap = coordMap,
                                         clearConfs = clearConfs,
                                         enforceChirality = True)
            if flag == -1:
                continue
            # set ff
            self._SetForceField(flag)
            # reorient core
            algMap = [(j, i) for i, j in enumerate(match)]
            AllChem.AlignMol(self.mol3Dx, core_mol, atomMap = algMap, maxIters = 200)
            # add tethers
            conf = core_mol.GetConformer(confId)
            for i in range(core_mol.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                pIdx = self._ff.AddExtraPoint(p.x, p.y, p.z, fixed = True) - 1
                self._ff.UFFAddDistanceConstraint(pIdx, match[i], False, 0, 0, 100.)
            # optimize
            self._ff.Initialize()
            self._ff.Minimize(maxIts = 1000)
            rms = AllChem.AlignMol(self.mol3Dx, core_mol, atomMap = algMap, maxIters = 200)
            # check chiral centers from 3D
            if not self._CheckStereoCA(flag):
                self.mol3Dx.RemoveConformer(flag)
                flag = -1
                continue
            # refresh other confs
            if clearConfs:
                self.mol.RemoveAllConformers()
                self.mol3D.RemoveAllConformers()
            # energy
            E = self._ff.CalcEnergy()
            conf3Dx = self.mol3Dx.GetConformer(flag)
            conf3Dx.SetDoubleProp('E', E)
            conf3Dx.SetDoubleProp('EmbedRMS', rms)
            # synchronize with mol3D
            conf3D = AllChem.Conformer()
            conf3D.SetDoubleProp('E', E)
            conf3D.SetDoubleProp('EmbedRMS', rms)
            for atom in self.mol3D.GetAtoms():
                idx = atom.GetIdx()
                conf3D.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            self.mol3D.AddConformer(conf3D, assignId = True)
            # synchronize with mol
            conf = AllChem.Conformer()
            conf.SetDoubleProp('E', E)
            conf.SetDoubleProp('EmbedRMS', rms)
            for atom in self.mol.GetAtoms():
                idx = atom.GetIdx()
                conf.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            self.mol.AddConformer(conf, assignId = True)        
        
        return flag
    
    
    def AddConstrainedConformerFromXYZ(self, pathCore, ignoreHs = True,
                                             clearConfs = True, maxAttempts = 10):
        '''
        Constrained embedding of the complex using geometry stored to XYZ file
        Central atom and donor atoms must be in complete accord
        '''
        core = ComplexFromXYZFile(pathCore)
        
        return self.AddConstrainedConformer(core, confId = 0, ignoreHs = ignoreHs,
                                            clearConfs = clearConfs, maxAttempts = maxAttempts)
    
    
    def GetNumConformers(self):
        '''
        Returns number of conformers
        '''
        
        return self.mol.GetNumConformers()
    
    
    def RemoveConformer(self, confId):
        '''
        Removes conformer with given idx
        '''
        self.mol.RemoveConformer(confId)
        self.mol3D.RemoveConformer(confId)
        self.mol3Dx.RemoveConformer(confId)
    
    
    def RemoveAllConformers(self):
        '''
        Removes all conformers
        '''
        self.mol.RemoveAllConformers()
        self.mol3D.RemoveAllConformers()
        self.mol3Dx.RemoveAllConformers()
    
    
    def AddConformers(self, numConfs = 10, clearConfs = True, maxAttempts = 10, rmsThresh = -1):
        '''
        Generates several conformers
        '''
        if self._PrintErrorInit():
            return None
        # generate 3D
        flags = []
        for i in range(numConfs):
            clearConfsIter = False if flags else clearConfs
            flag = self.AddConformer(clearConfs = clearConfsIter,
                                     maxAttempts = maxAttempts)
            # check flag and rms
            if flag == -1:
                continue
            if rmsThresh == -1:
                flags.append(flag)
                continue
            # check rms with previous conformers
            remove_conf = False
            for cid in flags:
                rms = AllChem.GetConformerRMS(self.mol3D, cid, flag)
                #print(rms)
                if rms < rmsThresh:
                    remove_conf = True
                    break
            if remove_conf:
                self.RemoveConformer(flag)
            else:
                flags.append(flag)
        
        return flags
    
    
    def AddConstrainedConformers(self, core, confId = 0, ignoreHs = True, numConfs = 10,
                                 clearConfs = True, maxAttempts = 10, rmsThresh = -1):
        '''
        Generates several conformers
        '''
        if self._PrintErrorInit():
            return None
        # generate 3D
        flags = []
        for i in range(numConfs):
            clearConfsIter = False if flags else clearConfs
            flag = self.AddConstrainedConformer(core, confId = confId, ignoreHs = ignoreHs,
                                                clearConfs = clearConfsIter,
                                                maxAttempts = maxAttempts)
            # check flag and rms
            if flag == -1:
                continue
            if rmsThresh == -1:
                flags.append(flag)
                continue
            # check rms with previous conformers
            remove_conf = False
            for cid in flags:
                rms = AllChem.GetConformerRMS(self.mol3D, cid, flag)
                #print(rms)
                if rms < rmsThresh:
                    remove_conf = True
                    break
            if remove_conf:
                self.RemoveConformer(flag)
            else:
                flags.append(flag)
        
        return flags
    
    
    def AddConstrainedConformersFromXYZ(self, pathCore, ignoreHs = True, numConfs = 10,
                                        clearConfs = True, maxAttempts = 10, rmsThresh = -1):
        '''
        Generates several conformers
        '''
        if self._PrintErrorInit():
            return None
        core = ComplexFromXYZFile(pathCore)
        # generate 3D
        flags = []
        for i in range(numConfs):
            clearConfsIter = False if flags else clearConfs
            flag = self.AddConstrainedConformer(core, confId = 0, ignoreHs = ignoreHs,
                                                clearConfs = clearConfsIter,
                                                maxAttempts = maxAttempts)
            # check flag and rms
            if flag == -1:
                continue
            if rmsThresh == -1:
                flags.append(flag)
                continue
            # check rms with previous conformers
            remove_conf = False
            for cid in flags:
                rms = AllChem.GetConformerRMS(self.mol3D, cid, flag)
                #print(rms)
                if rms < rmsThresh:
                    remove_conf = True
                    break
            if remove_conf:
                self.RemoveConformer(flag)
            else:
                flags.append(flag)
        
        return flags
    
    
    ######################
    # MolSimplify helper #
    ######################
    
    def GetBondedLigand(self, num):
        '''
        Extracts ligand with available 3D from the embedded complex
        num is an isotopic number of DA owned by target ligand
        '''
        if self._PrintErrorInit():
            return None
        N = self.mol3D.GetNumConformers()
        if not N:
            raise ValueError('Complex has no conformers')
        # is num in DAs
        idx_DA = None
        for DA, isotope in self._DAs.items():
            if isotope == num:
                idx_DA = DA
        if idx_DA is None:
            raise ValueError(f'Complex has donor atom with {num} order number')
        # remove dative bonds
        ed = Chem.EditableMol(self.mol3D)
        for DA in self._DAs:
            ed.RemoveBond(DA, self._idx_CA)
        mol = ed.GetMol()
        # remove unneeded fragments
        frags = Chem.GetMolFrags(mol)
        drop = []
        for frag in frags:
            if idx_DA not in frag:
                drop += list(frag)
        drop = sorted(drop, reverse = True)
        ed = Chem.EditableMol(mol)
        for idx in drop:
            ed.RemoveAtom(idx)
        mol = ed.GetMol()
        Chem.SanitizeMol(mol)
        
        return mol
    
    
    ##########
    # Output #
    ##########
    
    def _ConfToXYZ(self, confId):
        '''
        Generates text of XYZ file of conformer
        '''
        # coordinates
        xyz = []
        conf = self.mol3Dx.GetConformer(confId) # not mol3D as it uses in AlignMol
        for atom in self.mol3D.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol == '*':
                symbol = 'X'
            pos = conf.GetAtomPosition(atom.GetIdx())
            line = f'{symbol:2} {pos.x:>-10.4f} {pos.y:>-10.4f} {pos.z:>-10.4f}'
            xyz.append(line)
        # conf params
        E = conf.GetDoubleProp('E')
        rms = conf.GetDoubleProp('EmbedRMS')
        # mol smiles
        mol = deepcopy(self.mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        smiles = Chem.MolToSmiles(mol, canonical = False)
        # mol3D smiles
        mol3D = deepcopy(self.mol3D)
        for atom in mol3D.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        smiles3D = Chem.MolToSmiles(mol3D, canonical = False)
        # mol3Dx smiles
        mol3Dx = deepcopy(self.mol3Dx)
        for atom in mol3Dx.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        smiles3Dx = Chem.MolToSmiles(mol3Dx, canonical = False)
        # dummies' coords
        dummies = []
        conf3Dx = mol3Dx.GetConformer(confId)
        for i in range(mol3D.GetNumAtoms(), mol3Dx.GetNumAtoms()):
            if i < mol3D.GetNumAtoms():
                continue
            p = conf3Dx.GetAtomPosition(i)
            dummies += [p.x, p.y, p.z]
        # make text
        info = {'conf': confId, 'E': float(f'{E:.2f}'), 'rms': float(f'{rms:.4f}'),
                'geom': self._geom, 'smiles': smiles, 'smiles3D': smiles3D,
                'smiles3Dx': smiles3Dx, 'dummies': dummies}
        text = [str(len(xyz)), json.dumps(info)] + xyz
        
        return '\n'.join(text)+'\n'
    
    
    def ToXYZBlock(self, confId = -2):
        '''
        Returns XYZ as text block
        If confId is -1 or "min", conformer with the lowest energy will be saved
        If confId is -2 or "all", all conformers will be saved
        '''
        if self._PrintErrorInit():
            return None
        N = self.mol3D.GetNumConformers()
        if not N:
            raise ValueError('Bad conformer ID: complex has no conformers')
        # prepare conf idxs
        if confId in (-2, 'all'):
            conf_idxs = list(range(N))
        elif confId in (-1, 'min'):
            minE, minI = float('inf'), None
            for i in range(N):
                E = self.mol3D.GetConformer(i).GetDoubleProp('E')
                if E < minE:
                    minE = E
                    minI = i
            conf_idxs = [minI]
        else:
            conf_idxs = [confId]
        text = ''
        for conf_id in conf_idxs:
            text += self._ConfToXYZ(conf_id)
        
        return text
    
    
    def ToXYZ(self, path, confId = -2):
        '''
        Saves found conformers in XYZ format
        If confId is -1 or "min", conformer with the lowest energy will be saved
        If confId is -2 or "all", all conformers will be saved
        '''
        text = self.ToXYZBlock(confId)
        with open(path, 'w') as outf:
            outf.write(text)



