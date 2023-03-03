=====================================
Welcome to epic-mace's documentation!
=====================================

.. raw:: html
   :file: source/3D/epic_mace_3D.html


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Before you start
   
   source/when_to_use
   source/installation
   source/input
   source/view_structures


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: How To
   
   source/cookbook
   source/cli


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Documentation
   
   source/api
   source/changelog


**epic-mace** (MACE = MetAl Complexes Embedding) is an open source python library and command-line tool for the automated screening and discovery of metal complexes.
The software is developed by `Ivan Chernyshov`_ as part of the `Evgeny Pidko Group`_ in the `Department of Chemical Engineering`_ at `TU Delft`_.
Its features are to discover all possible configurations for square-planar and octahedral mononuclear metal complexes,
and generate atomic 3D coordinates suitable for quantum-chemical computations.
**epic-mace** shows high performance (> 99% success rate) for complexes of ligands of high denticity (up to 6), and is well-suited for the development of a massive computational pipelines aimed at solving problems of homogeneous catalysis.


Main features
=============

1. Stereomer search for octahedral and square-planar complexes.
2. Generation of 3D atomic coordinates, including instruments for conformer sampling.
3. Modification of ligands with predefined substituents.
4. Generation of geometry of coordinated ligands for `molSimplify`_.
5. Two available interfaces:

    - command-line interface for routine tasks;
    - Python package for organizing complex computational pipelines.


Useful links
============

1. `GitHub`_
2. `PyPI package`_
3. `conda package`_
4. `Performance testing`_
5. `CLI examples`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. _Ivan Chernyshov: https://github.com/IvanChernyshov
.. _Evgeny Pidko Group: https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group>
.. _Department of Chemical Engineering: https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/
.. _TU Delft: https://www.tudelft.nl/en/
.. _molSimplify: https://molsimplify.mit.edu/
.. _GitHub: https://github.com/EPiCs-group/epic-mace
.. _PyPI package: https://pypi.org/project/epic-mace/
.. _conda package: https://anaconda.org/grimgenius/epic-mace
.. _Performance testing: https://github.com/EPiCs-group/epic-mace/blob/master/performance/README.ipynb
.. _CLI examples: https://github.com/EPiCs-group/epic-mace/tree/master/examples
