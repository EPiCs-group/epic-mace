When to use
===========

.. note::
    
    **Long story short**: use *epic-mace* if your need:
      
      1. to perform stereoisomeric and conformational analysis of the ligand;
      
      2. to generate possible metal complex geometries for the computations of a catalytic cycle;
      
      3. to find the dependence of the properties of the ligand or complex on the substituents.


Any computational research starts with an atomic model. There are 3 main approaches for their preparation:

  - **manual preparation** — this most traditional way has some disadvantages:
    
    1. it cannot be automated and thus cannot be used in high-throughput computational pipelines;
    
    2. the quality of the structure can greatly depend on the experience of the computational chemist introducing an additional bias into results. However, that problem can be overcame by using conformational sampling techniques.
  
  - **use experimental geometries** — the shortcomings of this approach overlap significantly with manual generation:
    
    1. not all molecular systems have experimental geometries (mainly XRD, NMR for biomolecules);
    
    2. using the experimental geometry without additional conformational sampling may lead to an inadequate results, especially in catalysis where kinetics, not thermodynamics, plays a decisive role.
  
  - **automated 3D generation** — the ideal option if the software does not have significant flaws.


The rapid development of chemoinformatics over the past twenty years has given us a wide range for the automatic generation of atomic models of a wide range of chemical systems, including small molecules, biomolecules, and some inorganic systems (like crystal surfaces, amorphous solids, etc.). however, support for metal complexes has appeared only in recent years:

  1. `molSimplify`_ has an excellent library of ligands and tools for generating metal complexes of various molecular geometries (in addition to a wide range of other features). However, it generates the most stable stereomers and has no possibility to generate all of them which can be crucial for the catalysis. In addition, it often fails to generate atomic coordinates for ligands with a denticity of three or more (unless their geometry in coordinated state were pre-generated).

  2. `RDKit 2022.09`_ version adds support of octahedral, square-planar, and trigonal bipyramidal geometries. However, at the moment it cannot be used either to search for stereomers or to generate the correct 3D atomic coordinates.


And here comes **epic-mace**! Its main functionality is to generate possible stereomers for mononuclear square-planar and octahedral complexes, which are the main system of a homogenious catalysis. Therefore, if your need:

  1. to perform stereoisomeric and conformational analysis of the ligand;
  
  2. to generate possible metal complex geometries for the computations of a catalytic cycle;
  
  3. to find the dependence of the properties of the ligand or complex on the substituents,

then **epic-mace** is the best tool to automatically generate initial geometries for further computations and molecular modeling.


.. _molSimplify: https://molsimplify.mit.edu/
.. _RDKit 2022.09: https://www.rdkit.org/docs/RDKit_Book.html#octahedral


