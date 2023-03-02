Input SMILES
============

.. warning::
    This section is crucial for understanding how to get the input for any program that uses **epic-mace**,
    be it a python library, a CLI tool, or a GUI program.


SMILES specifications
---------------------

The main problem that the **epic-mace** package solves is a stereomer search and generation of 3D atomic coordinates for mononuclear octahedral and square-planar metal complexes.
Thus, any program that uses the **mace** package must receive the structure of the complex or the corresponding ligands as input. **mace** uses the `SMILES`_ notation
to describe structure of input molecules, or rather two versions: `RDKit SMILES extension`_ and `ChemAxon SMILES`_.
The reason is that these are fairly popular SMILES formats that support dative bonds used in `mace` to describe metal-ligand bonds.


SMILES requirements
-------------------

SMILES of complex or ligands must satisfy the following conditions:

    1. Donor atoms have non-zero atom map numbers describing their spatial arrangement around a central atom:
    
      .. image:: images/geoms_scheme.png
    
      If one aims to generate all stereomers and is not interested in specifying the stereochemistry of the complex, when any non-zero map number can be used.
    
    2. Non-**DA** atoms must have zero map numbers.
    
    3. Bonds between central (**CA**) and donor (**DA**) atoms are dative, directed from **DA** to **CA**.
    
    4. The complex contains one central atom only.
    
    5. Substituents can be represented by dummy atoms with isotope number corresponding to the number of substituent (or R-groups for ChemAxon SMILES). Other atoms must have zero isotope numbers.
    
    6. Substituents must be monovalent (*may change in the future*) and cannot be bound to the central atom.


ChemAxon Marvin
---------------

For those who is familiar with the practical side of chemoinformatics, getting SMILES will not be a problem.
For others the easiest way is to use the `ChemAxon Marvin`_ chemical editor:

    1. To set map numbers, right click on the atom and select the **map** option.

    2. To add dummy atoms describing substituents, use the **R-group** option.

    3. To save the structure in SMILES format, select the molecule, right click on it and select **Save as**. Use ChemAxon SMILES for metal complexes and Daylight SMILES for ligands.

    .. image:: images/marvin_get_smiles.png

    4. Now the SMILES string is copied to the clipboard and can be pasted via the **Ctrl+V** keyboard shortcut.

.. _SMILES: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
.. _RDKit SMILES extension: https://www.rdkit.org/docs/RDKit_Book.html#smiles-support-and-extensions
.. _ChemAxon SMILES: https://docs.chemaxon.com/display/docs/smiles.md
.. _ChemAxon Marvin: https://chemaxon.com/marvin


