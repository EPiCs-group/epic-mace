Changelog
=========

0.5.0
-----

The first tracked version. The main features are:
  
  - The :class:`mace.Complex` class describing mononuclear octahedral and square-planar metal complexes:
  
    - the :meth:`mace.Complex.GetStereomers` method for stereomer search;
    
    - :meth:`mace.Complex.AddConformers` and :meth:`mace.Complex.AddConstrainedConformers` methods for generation of 3D atomic coordinates;
    
    - the :meth:`mace.Complex.GetRepresentativeConfs` method for filtering representative set of conformers;
    
    - the :meth:`mace.Complex.GetBondedLigand` method extracting geometry of coordinated ligand (to use in molSimplify).
  
  - The :func:`mace.ComplexFromXYZFile` function reading complex from the **epic-mace**-generated xyz-file.
  
  - The :func:`mace.AddSubsToMol` function functionalizing ligands and complexes by adding substituents.
  
  - Command-line interface:
    
    - **epic-mace**: supports all main features of the **mace** package;
    
    - **epic-mace-quickstart**: generates templates for **epic-mace**'s input and substituents files.


