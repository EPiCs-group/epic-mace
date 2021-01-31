'''
Functions for SMILES parsing
'''

#%% Imports

import json, os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point3D

from .complex_init_mols import ComplexFromMol


#%% Functions

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
    # make mol
    maps = []
    for atom in mol.GetAtoms():
        maps.append(atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    mol = Chem.RenumberAtoms(mol, tuple(zip(*sorted([(j, i) for i, j in enumerate(maps)])))[1])
    # check mol mapping
    flag = False
    for atom in mol.GetAtoms():
        sym = '*' if atoms[atom.GetIdx()] == 'X' else atoms[atom.GetIdx()]
        if atom.GetSymbol() != sym:
            flag = True
            break
    if flag:
        raise ValueError('Bad XYZ file: atom numbering in SMILES and coordinates block do not match')
    # restore DAs in mol
    for a in mol.GetAtoms():
        if not a.GetIsotope():
            continue
        if 'DATIVE' not in [str(b.GetBondType()) for b in a.GetBonds()]:
            continue
        a.SetAtomMapNum(a.GetIsotope())
    # make mol3D
    maps3D = []
    for atom in mol3D.GetAtoms():
        maps3D.append(atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    mol3D = Chem.RenumberAtoms(mol3D, tuple(zip(*sorted([(j, i) for i, j in enumerate(maps3D)])))[1])
    # check mol3D mapping
    flag = False
    for atom in mol3D.GetAtoms():
        sym = '*' if atoms[atom.GetIdx()] == 'X' else atoms[atom.GetIdx()]
        if atom.GetSymbol() != sym:
            flag = True
            break
    if flag:
        raise ValueError('Bad XYZ file: atom numbering in SMILES3D and coordinates block do not match')
    # restore DAs in mol3D
    for a in mol3D.GetAtoms():
        if not a.GetIsotope():
            continue
        if 'DATIVE' not in [str(b.GetBondType()) for b in a.GetBonds()]:
            continue
        a.SetAtomMapNum(a.GetIsotope())
    # make mol3Dx
    maps3Dx = []
    for atom in mol3Dx.GetAtoms():
        maps3Dx.append(atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
    mol3Dx = Chem.RenumberAtoms(mol3Dx, tuple(zip(*sorted([(j, i) for i, j in enumerate(maps3Dx)])))[1])
    # check mol3Dx mapping
    flag = False
    for atom in mol3Dx.GetAtoms():
        if atom.GetIdx() < mol3D.GetNumAtoms():
            sym = '*' if atoms[atom.GetIdx()] == 'X' else atoms[atom.GetIdx()]
            if atom.GetSymbol() != sym:
                flag = True
                break
        else:
            if atom.GetSymbol() != '*':
                flag = True
                break
    if flag:
        raise ValueError('Bad XYZ file: atom numbering in SMILES3Dx and coordinates block do not match')
    # restore DAs in mol3Dx
    for a in mol3Dx.GetAtoms():
        if not a.GetIsotope():
            continue
        if 'DATIVE' not in [str(b.GetBondType()) for b in a.GetBonds()]:
            continue
        a.SetAtomMapNum(a.GetIsotope())
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
    X.mol3Dx = mol3Dx
    X._SetEmbedding()
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


