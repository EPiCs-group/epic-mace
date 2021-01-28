'''
Functions reading aliases info
'''

# imports
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../mace')))
from Complex import MolFromSmiles

# path to default aliases
PATH_PACKAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../aliases')


#%% Functions

def GetAliasFiles():
    '''
    Return paths to files containing Rs & ligands aliases
    '''
    global PATH_PACKAGE
    # get path var
    PATH_MACE = os.environ.get('PATH_MACE')
    if PATH_MACE:
        dirs = [('user', _) for _ in PATH_MACE.split(os.pathsep)]
    # add wd and internal and script's paths
    dirs = [('wd', os.getcwd())] + dirs + [('mace', PATH_MACE)]
    # check files existence
    Rs, ligands = None, None
    for flag, path_dir in dirs:
        if Rs is None and os.path.isfile(os.path.join(path_dir, 'Rs.txt')):
            Rs = (flag, os.path.join(path_dir, 'Rs.txt'))
        if ligands is None and os.path.isfile(os.path.join(path_dir, 'ligands.txt')):
            ligands = (flag, os.path.join(path_dir, 'ligands.txt'))
    
    return Rs, ligands


def ReadParamFile(path):
    '''
    Extracts dictionary with parameters from "key: val"-formated file
    '''
    with open(path, 'r') as inpf:
        text = [_.strip() for _ in inpf.readlines()]
    text = [_ for _ in text if _]
    # make dict
    info = {}
    for i, line in enumerate(text):
        if ':' not in line:
            raise ValueError('Bad alias file format: bad line format:\n\nLine #{i+1}: {line}\n\nIt must be in "Var: Value" format')
        idx = line.index(':')
        key = line[:idx].strip()
        val = line[idx+1:].strip()
        if key in info:
            raise ValueError('Bad alias file format: key {key} is defined several times')
        info[key] = val
    
    return info


def GetRs(path_Rs):
    '''
    Returns dictionary of predefined substituents
    '''
    # read file
    info = ReadParamFile(path_Rs)
    # check Rs
    for name, smiles in info.items():
        R = MolFromSmiles(smiles)
        if R is None:
            raise ValueError(f'Bad Rs aliases file format: SMILES of {name} is not readable: {smiles}')
        dummies = [_ for _ in R.GetAtoms() if _.GetSymbol() == '*']
        if len(dummies) != 1:
            raise ValueError(f'Bad Rs aliases file format: SMILES of {name} must contain exactly one dummy atom: {smiles}')
        if len(dummies[0].GetNeighbors()) != 1:
            raise ValueError(f'Bad Rs aliases file format: SMILES of {name} must contain dummy atom bonded to exactly one atom by single bond: {smiles}')
        if str(dummies[0].GetBonds()[0].GetBondType()) != 'SINGLE':
            raise ValueError(f'Bad Rs aliases file format: SMILES of {name} must contain dummy atom bonded to exactly one atom by single bond: {smiles}')
    
    return info


def GetLigands(path_ligands):
    '''
    Returns dictionary of predefined ligands
    '''
    # read file
    info = ReadParamFile(path_ligands)
    # check ligands
    for name, smiles in info.items():
        mol = MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f'Bad ligand aliases file format: SMILES of {name} is not readable: {smiles}')
        DAs = [_.GetIdx() for _ in mol.GetAtoms()]
        if not DAs:
            raise ValueError(f'Bad ligands aliases file format: SMILES of {name} must countain at least one donor atom')
    
    return info


def GetAliases():
    '''
    Return dictionaries containing info on aliases
    '''
    Rs, ligands = GetAliasFiles()
    if Rs is None:
        print(f'Warning: Rs aliases file was not found anywhere')
    else:
        Rs = GetRs(Rs[1])
        print(f'Rs')
    if ligands is None:
        pass
    
    return Rs, ligands


