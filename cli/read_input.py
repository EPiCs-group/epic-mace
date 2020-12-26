'''
Functions for reading and interpreting input files
'''

import sys, os

sys.path.insert(0, os.path.join(os.path.abspath(os.path.dirname(__file__)), '../mace'))
from Complex import (
    MolFromSmiles, AddSubsToMol,
    ComplexFromMol, ComplexFromLigands,
    ComplexFromXYZFile
)



#%% Functions

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



