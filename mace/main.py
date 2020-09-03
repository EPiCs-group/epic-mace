'''
Command-line interface for generation of 3D coordinates of metal complexes
'''

from .complex import (
    MolFromSmiles, AddSubsToMol,
    ComplexFromMol, ComplexFromLigands,
    ComplexFromXYZFile,
    Complex
)

import os


#%% Input

def ParseInput(path):
    '''
    Extracts task info from input file
    '''
    # read text
    if not os.path.isfile(path):
        raise ValueError('Bad XYZ path: file does not exist')
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
    info = {}
    for s, e in zip(starts, ends):
        info[text[s][1:]] = text[s+1:e]
    
    return info





#%% 










