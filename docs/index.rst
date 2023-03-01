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
   :caption: epic-mace
   
   source/cookbook
   source/cli
   source/api
   source/changelog



**epic-mace** or MACE (MetAl Complexes Embedding) is an open source python library
for the automated screening and discovery of metal complexes.
MACE is developed by `Ivan Chernyshov`_ as part of the `Evgeny Pidko Group`_ in the `Department of Chemical Engineering`_ at `TU Delft`_.
Its features are to discover all possible configurations for square-planar and octahedral metal complexes,
and generate atomic 3D coordinates suitable for quantum-chemical computations.
MACE shows high performance (> 99% success rate) for complexes of ligands of high denticity (up to 6), 
and thus is well-suited for the development of a massive computational pipelines aimed at solving problems of homogeneous catalysis.


Main features
=============

1. stereomer search for SP/OH
2. 3D generation, including instruments for conformer sampling
3. CLI
4. `performance`_


Useful links
============

1. GitHub
2. `PyPI`_ package
3. `conda`_ package
4. `performance`_
5. CLI examples


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. _Ivan Chernyshov: https://github.com/IvanChernyshov
.. _Evgeny Pidko Group: https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group>
.. _Department of Chemical Engineering: https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/
.. _TU Delft: https://www.tudelft.nl/en/
.. _conda: https://anaconda.org/grimgenius/epic-mace
.. _PyPI: https://pypi.org/project/epic-mace/
.. _performance: https://github.com/EPiCs-group/epic-mace/blob/master/performance/README.ipynb
