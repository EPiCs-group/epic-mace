'''This package contains classes and functions to discover possible stereomers
and generate 3D atomic coordinates for mononuclear metal complexes.

Concept
-------

The :class:`mace.Complex` class contains all required functionality for stereomer
search and 3D embedding. The corresponding object can be initialized using the
:class:`mace.Complex` constructor, or :func:`mace.ComplexFromMol`, :func:`mace.ComplexFromLigands`,
and :func:`mace.ComplexFromXYZFile` functions.

Built-in coordination geometries currently include octahedral (``OH``),
square-planar (``SP``), tetrahedral (``TET``), square-pyramidal (``SPY``),
trigonal-bipyramidal (``TBP``), and sandwich (``SAN``).

An initialized complex object may not have defined stereochemestry of the central atom (non-empty
:attr:`mace.Complex.err_init`). In this case, the :meth:`mace.Complex.GetStereomers` method
should be used to get possible stereomers. 3D atomic coordinates can be generated using
:meth:`mace.Complex.AddConformer` and :meth:`mace.Complex.AddConformers` methods.

For other features of the MACE package see the tutorial.

Classes & Functions
-------------------
'''

# imports
from ._smiles_parsing import MolFromSmiles, MolToSmiles
from ._substituents import AddSubsToMol
from ._complex_object import Complex
from ._complex_init_mols import ComplexFromMol, ComplexFromLigands
from ._complex_init_files import ComplexFromXYZFile

# package info
__author__ = "Ivan Yu. Chernyshov"
__email__ = "ivan.chernyshoff@gmail.com"
__version__ = '0.5.0'

# module functions
__all__ = [
    'Complex',
    'ComplexFromMol', 'ComplexFromLigands', 'ComplexFromXYZFile',
    'MolFromSmiles', 'MolToSmiles', 'AddSubsToMol'
]

# disable logger
from rdkit import RDLogger as _RDLogger
_RDLogger.DisableLog('rdApp.*')

