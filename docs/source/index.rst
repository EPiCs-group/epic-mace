=====================================
Welcome to epic-mace's documentation!
=====================================

.. raw:: html
   :file: 3D/epic_mace_3D.html

**epic-mace** or MACE (MetAl Complexes Embedding) is an open source python library
for the automated screening and discovery of metal complexes.
MACE is developed by `Ivan Chernyshov`_ as part of the `Evgeny Pidko Group`_ in the `Department of Chemical Engineering`_ at `TU Delft`_.
Its features are to discover all possible configurations for square-planar and octahedral metal complexes,
and generate atomic 3D coordinates suitable for quantum-chemical computations.
MACE shows high performance (> 99% success rate) for complexes of ligands of high denticity (up to 6), 
and thus is well-suited for the development of a massive computational pipelines aimed at solving problems of homogeneous catalysis.
    
.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Contents:
   
   installation
   input
   tutorial
   api
   gui


Requirements
============

* Python 3.7 or higher (Python 3.7 is recommended);
* `RDKit`_ 2020.09 or higher (RDKit 2020.09 is a **must** for a correct functioning).


Performance
===========

MACE shows high performance (> 99% success rate) for complexes of ligands, extracted from Cambridge Structural Database. For more details see `performance`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. _Ivan Chernyshov: https://github.com/IvanChernyshov
.. _Evgeny Pidko Group: https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group>
.. _Department of Chemical Engineering: https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/
.. _TU Delft: https://www.tudelft.nl/en/
.. _RDKit: https://www.rdkit.org/
.. _performance: https://github.com/EPiCs-group/epic-mace/blob/master/performance/README.ipynb
