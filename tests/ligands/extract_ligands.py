'''
Prepares ligands from CSD complexes encoded as SMILES
'''

#%% Imports

import os

from rdkit import Chem


#%% Functions

def read_smi(path, metals, CNs):
    '''
    Extracts SMILES of mono-metal complexes from ConQuest output
    '''
    # read text
    with open(path, 'r') as inpf:
        text = [_.strip() for _ in inpf.readlines()]
        text = [_ for _ in text if _]
    # set SMILES parsing
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    ps.sanitize = False
    # split lines
    ligand_elems = ['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'I']
    complexes = []
    for line in text:
        smiles, refcode = line.split('\t')
        for smi in smiles.split('.'):
            mol = Chem.MolFromSmiles(smi, params = ps)
            ms = [a for a in mol.GetAtoms() if a.GetSymbol() in metals]
            CA = ms[0] if ms else None
            n = len(ms)
            nd = len([a for a in mol.GetAtoms() if a.GetSymbol() not in ligand_elems])
            if n == 1 and nd == 1 and len(CA.GetNeighbors()) in CNs:
                complexes.append( (refcode, smi) )
    
    return complexes


def extract_ligands(smiles):
    '''
    Extracts ligands from metal complex
    '''
    # read mol
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    ps.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params = ps)
    # find CA
    ligand_elems = ['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'I']
    CA = [a for a in mol.GetAtoms() if a.GetSymbol() not in ligand_elems][0]
    # label DAs
    for a in CA.GetNeighbors():
        a.SetAtomMapNum(1)
    # remove bonds
    ed_mol = Chem.EditableMol(mol)
    bonds = [(_.GetBeginAtomIdx(), _.GetEndAtomIdx()) for _ in CA.GetBonds()]
    for i, j in bonds:
        ed_mol.RemoveBond(i, j)
    # get fragments
    ligands = []
    CA.SetAtomicNum(0)
    for idxs in Chem.GetMolFrags(ed_mol.GetMol()):
        idxs = set([*idxs, CA.GetIdx()])
        ed_frag = Chem.EditableMol(mol)
        drop = sorted([i for i in range(mol.GetNumAtoms()) if i not in idxs], reverse = True)
        for idx in drop:
            ed_frag.RemoveAtom(idx)
        ligand = ed_frag.GetMol()
        if len([a for a in ligand.GetAtoms() if a.GetAtomMapNum()]) > 1:
            ligands.append(Chem.MolToSmiles(ligand))
    
    return ligands


def get_valence(atom):
    '''
    Gets valence for atom in non-sanitized RDKit Mol
    '''
    bonds = {'SINGLE': 0, 'DOUBLE': 0, 'TRIPLE': 0, 'AROMATIC': 0}
    idx_CA = [_.GetIdx() for _ in atom.GetNeighbors() if _.GetSymbol() == '*'][0]
    val = 0
    for b in atom.GetBonds():
        if idx_CA in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()):
            val += 1 # skipped dative bond
            continue
        bt = str(b.GetBondType())
        if bt not in bonds:
            print(bt)
            return None
        bonds[bt] += 1
    if bonds['AROMATIC'] == 3:
        print('wow')
        return 4
    elif bonds['AROMATIC'] == 2:
        val += 3
    val += bonds['SINGLE'] + 2*bonds['DOUBLE'] + 3*bonds['TRIPLE'] + atom.GetTotalNumHs()
    
    return val


def pretify_ligand(smiles):
    '''
    Makes reasonable H-count and charge for each donor atom
    '''
    # read mol
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    ps.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params = ps)
    # find CA and DAs
    CA = [a for a in mol.GetAtoms() if a.GetSymbol() == '*'][0]
    DAs = [a for a in mol.GetAtoms() if a.GetAtomMapNum()]
    # check and modify DAs' properties
    for DA in DAs:
        # DA info
        sym = DA.GetSymbol()
        nR = DA.GetNumRadicalElectrons()
        val = get_valence(DA)
        # carbene
        if sym == 'C' and val == 3 and nR == 1:
            DA.SetNumRadicalElectrons(2)
        # carbene or carbanion
        elif sym == 'C' and val == 4:
            bonds = DA.GetBonds()
            bts = [str(b.GetBondType()) for b in bonds]
            if 'DOUBLE' in bts:
                idx = bts.index('DOUBLE')
                b = bonds[idx]
                if b.GetOtherAtom(DA).GetSymbol() == 'N':
                    b.SetBondType(Chem.BondType.SINGLE)
                    DA.SetNumRadicalElectrons(2)
                elif b.GetOtherAtom(DA).GetSymbol() == 'C':
                    DA.SetNumRadicalElectrons(0)
                    DA.SetFormalCharge(-1)
            elif 'AROMATIC' in bts:
                ns = [a.GetSymbol() for a in DA.GetNeighbors()]
                if 'N' not in ns:
                    DA.SetNumRadicalElectrons(0)
                    DA.SetNumExplicitHs(0)
                    DA.SetFormalCharge(-1)
        # nitrogen
        elif sym == 'N' and val == 3:
            DA.SetFormalCharge(-1)
        # oxygen and sulfur anions
        elif sym in ('O', 'S') and val == 2:
            DA.SetFormalCharge(-1)
    # break bonds
    ed = Chem.EditableMol(mol)
    for DA in DAs:
        ed.RemoveBond(CA.GetIdx(), DA.GetIdx())
    ed.RemoveAtom(CA.GetIdx())    
    ligand = ed.GetMol()
    
    return Chem.MolToSmiles(ligand)


def is_good_ligand(smiles):
    '''
    Drops "bad" ligands
    '''
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False
    # check atoms
    for a in mol.GetAtoms():
        # info
        isDA = bool(a.GetAtomMapNum())
        sym = a.GetSymbol()
        nH = a.GetTotalNumHs()
        val = a.GetTotalValence()
        Z = a.GetFormalCharge()
        nR = a.GetNumRadicalElectrons()
        # radicals (not carbenes)
        if nR == 1:
            return False
        # carbene not DA
        if nR == 2 and not isDA:
            return False
        # PH/SH
        if sym in ('P', 'S') and nH:
            return False
        # N/P no charge
        if sym in ('N', 'P') and val == 4 and not Z:
            return False
        # double charge
        if abs(Z) > 1:
            return False
        # DA-C
        if sym == 'C' and isDA:
            # only carbenes and carbanions are allowed
            if not ( (val == 2 and nR == 2 and Z == 0) or \
                     (val ==3 and nR == 0 and Z == -1) ):
                return False
            # not N-carbenes
            if val == 2 and nR == 2 and Z == 0:
                if 'N' not in [_.GetSymbol() for  _ in a.GetNeighbors()]:
                    return False
    
    return True



#%% Main code

# paths
path_dir = os.path.dirname(__file__)
path_sp = os.path.join(path_dir, 'csd/SP.smi')
path_mn = os.path.join(path_dir, 'csd/OH_Mn.smi')
path_ru = os.path.join(path_dir, 'csd/OH_Ru.smi')
path_out = os.path.join(path_dir, 'cleared_ligands/csd_ligands.csv')
path_png = os.path.join(path_dir, 'cleared_ligands')

# read data
SP = read_smi(path_sp, ['Pd', 'Pt'], [4])
OH = read_smi(path_mn, ['Mn'], [6]) + read_smi(path_ru, ['Ru'], [6])

# SP
spls = []
for refcode, smiles in SP:
    spls += [(refcode, ligand) for ligand in extract_ligands(smiles)]
spls = [(refcode, pretify_ligand(smiles)) for refcode, smiles in spls]
spls = [(refcode, smiles) for refcode, smiles in spls if is_good_ligand(smiles)]
spls_unique = {}
for refcode, smiles in spls:
    new_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    if new_smiles in spls_unique:
        spls_unique[new_smiles].append(refcode)
    else:
        spls_unique[new_smiles] = [refcode]

# OH
ohls = []
for refcode, smiles in OH:
    ohls += [(refcode, ligand) for ligand in extract_ligands(smiles)]
ohls = [(refcode, pretify_ligand(smiles)) for refcode, smiles in ohls]
ohls = [(refcode, smiles) for refcode, smiles in ohls if is_good_ligand(smiles)]
ohls_unique = {}
for refcode, smiles in ohls:
    new_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    if new_smiles in ohls_unique:
        ohls_unique[new_smiles].append(refcode)
    else:
        ohls_unique[new_smiles] = [refcode]

# make table with ligands info
text = ['id,geom,n,smiles,refcodes']
for i, (smiles, refcodes) in enumerate(spls_unique.items()):
    mol = Chem.MolFromSmiles(smiles)
    refcodes = ';'.join(refcodes)
    n = 0
    for a in mol.GetAtoms():
        if a.GetAtomMapNum():
            n += 1
    text.append(f'SP{i},SP,{n},{smiles},{refcodes}')
# same for OH
for i, (smiles, refcodes) in enumerate(ohls_unique.items()):
    mol = Chem.MolFromSmiles(smiles)
    refcodes = ';'.join(refcodes)
    n = 0
    for a in mol.GetAtoms():
        if a.GetAtomMapNum():
            n += 1
    text.append(f'OH{i},OH,{n},{smiles},{refcodes}')
# save file
with open(path_out, 'w') as outf:
    outf.write('\n'.join(text)+'\n')

# check ligands
from rdkit.Chem import Draw
for i, (smiles, refcodes) in enumerate(spls_unique.items()):
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, f'{path_png}/square_planar/{i}.png', (300, 300) )
for i, (smiles, refcodes) in enumerate(ohls_unique.items()):
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, f'{path_png}/octahedral/{i}.png', (300, 300) )


