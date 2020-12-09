'''
Command-line interface for generation of 3D coordinates of metal complexes
'''

import sys, os
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from Complex import (
    MolFromSmiles, AddSubsToMol,
    ComplexFromMol, ComplexFromLigands,
    ComplexFromXYZFile
)


#%% Globals

path_aliases = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../aliases')
path_Rs = os.path.join(path_aliases, 'Rs.txt')
path_ligands = os.path.join(path_aliases, 'Ligands.txt')



#%% Aliases

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


def GetRs():
    '''
    Returns dictionary of predefined substituents
    '''
    global path_Rs
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


def GetLigands():
    '''
    Returns dictionary of predefined ligands
    '''
    global path_ligands
    # read file
    info = ReadParamFile(path_ligands)
    # check ligands
    for name, smiles in info.items():
        mol = MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f'Bad ligand aliases file format: SMILES of {name} is not readable: {smiles}')
        DAs = [_.GetIdx() for _ in []]
        pass
    
    return info



#%% Input

def GetInputParams(path):
    '''
    Extracts task info from input file
    '''
    global path_defaults, path_Rs, path_ligands
    # read text
    if not os.path.isfile(path):
        raise ValueError('Bad input file path: file does not exist')
    with open(path, 'r') as inpf:
        text = [line.strip() for line in inpf.readlines()]
    text = [line for line in text if line]
    text = [line for line in text if line[0] == '#']
    # find block params
    starts = []
    ends = []
    for i, line in text:
        if line[0] == '%':
            starts.append(i)
        elif line == 'end':
            ends.append(i)
    # check starts/ends
    if len(starts) != len(ends):
        raise ValueError('Bad input file format: different number of info blocks\' starts and ends')
    for s, e in zip(starts, ends):
        if s > e:
            raise ValueError('Bad input file format: the info block\'s end is placed before its start')
    for s, e in zip(starts[1:], ends[:-1]):
        if s < e:
            raise ValueError('Bad input file format: the new info block\'s starts is placed before the end of the previous block')
    # extract params
    blocks = {}
    for s, e in zip(starts, ends):
        blocks[text[s][1:].lower()] = text[s+1:e]
    
    return blocks


def ParseInput(path):
    '''
    Extracts info from input blocks
    '''
    info = GetInputParams(path)
    # check input parameters
    keys = ('numconfs', 'rmsthresh', 'lowestconf', 'stereomers', 'pathcore', 
            'pathrs', 'pathligands', 'complex', 'ligands', 'ca', 'rs')
    bad_keys = list(set(info.keys()).difference(set(keys)))
    if bad_keys:
        raise ValueError(f'Bad task block format: unknown variables: {", ".join(bad_keys)}')
    # number of conformers
    if 'numconfs' in info:
        try:
            numConfs = int(info['numconfs'])
        except ValueError:
            raise ValueError(f'Bad task block format: numConfs ({info["numconfs"]}) variable must be positive integer')
        if numConfs <= 0:
            raise ValueError(f'Bad task block format: numConfs ({info["numconfs"]}) variable must be positive integer')
        info['numconfs'] = numConfs
    # rms threshold
    if 'rmsthresh' in info:
        pass
    # return only the lowest conformer?
    if 'lowestconf' in info:
        pass
    # should we generate stereomers for the complex?
    if 'stereomers' in info:
        pass
    # constrained embedding
    if 'pathcore' in info:
        pass
    # aliases paths
    if 'pathrs' in info:
        pass
    if 'pathligands' in info:
        pass
    # structure
    if 'complex' in info:
        pass
    if 'ligands' in info:
        pass
    if 'ca' in info:
        pass
    if 'rs' in info:
        pass
    
    return info


def SanitizeInput(info):
    '''
    Checks extracted parameters and sanitize them using defaults
    '''
    # check aliases
    
    
    # check structure
    
    
    # check 3D gen params
    
    
    
    return info


def PrepareRun(path):
    '''
    Prepares data for the MACE script run
    '''
    
    
    return



#%% Main function

def main():
    '''
    Main function of the command-line tool
    '''
    # get script params
    
    
    # get input
    # info = ParseInput(path)
    # info = SanitizeInput(info)
    
    
    
    return



#%% Main code

if __name__ == '__main__':
    
    main()



