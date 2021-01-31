'''
Functions for SMILES parsing
'''

#%% Imports

from rdkit import Chem

from .smiles_parsing import MolFromSmiles
from .complex_object import Complex


#%% Functions

def ComplexFromMol(mol, geom, maxResonanceStructures = 1):
    '''
    Initializes Complex from RDKit Mol. Use it if mol->smiles transform is unwanted,
    e.g. if molecule contains *=DA->CA fragment with stereospecified double bond
    '''
    X = Complex('[*:1]->[*]', geom, maxResonanceStructures)
    if not mol or not Chem.MolToSmiles(mol):
        raise ValueError('Bad molecule: None or not RDKit Mol')
    X.mol = mol
    X._CheckMol()
    X._SetComparison()
    X.mol3D = Chem.AddHs(X.mol)
    X._embedding_prepared = False
    X._ff_prepared = False
    
    return X


def ComplexFromLigands(ligands, CA, geom, maxResonanceStructures = 1):
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
    DAs = [atom for atom in mol.GetAtoms() if atom.GetAtomMapNum()]
    CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
    CTs = [Chem.rdchem.BondStereo.STEREOCIS, Chem.rdchem.BondStereo.STEREOTRANS]
    dummies = {}
    doubles = []
    for DA in DAs:
        idx_DA = DA.GetIdx()
        # is it DA->[*] fragment
        idx_dummy = None
        for i, b in enumerate(DA.GetBonds()):
            if str(b.GetBondType()) != 'DATIVE':
                continue
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
    
    return ComplexFromMol(mol, geom, maxResonanceStructures)



