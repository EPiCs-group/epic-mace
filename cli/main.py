'''
Command-line interface for generation of 3D coordinates of metal complexes
'''

import sys, os

sys.path.insert(0, os.path.join(os.path.abspath(os.path.dirname(__file__)), '../mace'))
from Complex import (
    MolFromSmiles, AddSubsToMol,
    ComplexFromMol, ComplexFromLigands,
    ComplexFromXYZFile
)

from aliases import GetAliases
from read_input import ParseInput, SanitizeInput, PrepareRun


#%% Functions

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


