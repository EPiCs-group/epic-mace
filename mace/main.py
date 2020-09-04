'''
Command-line interface for generation of 3D coordinates of metal complexes
'''

import sys, os
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from Complex import (
    MolFromSmiles, AddSubsToMol,
    ComplexFromMol, ComplexFromLigands,
    ComplexFromXYZFile,
    Complex
)


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


def ReadInput(path):
    '''
    Extracts work parameters from the input file
    '''
    info = ParseInput(path)
    # check must-have parameters
    for block in ('structure', 'task'):
        if block not in path:
            raise ValueError(f'Bad input file format: {block} block is missing')
    # task
    task = ''
    params = {'task': task}
    
    return params



#%% Ligands and Rs

path_library = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../aliases')

def GetRs():
    '''
    Returns dictionary of predefined substituents
    '''
    global path_library
    # read file
    with open(os.path.join(path_library, 'Rs.txt'), 'r') as inpf:
        text = [_.strip() for _ in inpf.readlines()]
    text = [_ for _ in text if _]
    # make dict
    info = {}
    for line in text:
        idx = line.index(':')
        name = line[:idx].strip()
        smiles = line[idx+1:].strip()
        info[name] = smiles
    
    return info


def GetLigands():
    '''
    Returns dictionary of predefined ligands
    '''
    global path_library
    # read file
    with open(os.path.join(path_library, 'Ligands.txt'), 'r') as inpf:
        text = [_.strip() for _ in inpf.readlines()]
    text = [_ for _ in text if _]
    # make dict
    info = {}
    for line in text:
        idx = line.index(':')
        name = line[:idx].strip()
        smiles = line[idx+1:].strip()
        info[name] = smiles
    
    return info



#%% Main function

def main():
    '''
    Main function of the command-line tool
    '''
    
    return



#%% Main code

if __name__ == '__main__':
    
    main()



