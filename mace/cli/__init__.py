'''Module containing command-line scripts for massive generation of 3D coordinates
for mononuclear octahedral and square-planar metal complexes
'''

from ._quickstart import main as quickstart
from ._generator import main as generator

__all__ = ['quickstart', 'generator']


