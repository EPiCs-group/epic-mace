Installation
============

.. warning::
    **epic-mace** requires Python 3.7 and the 2020.09 version of `RDKit`_ for a correct functioning!
    
    Earlier versions do not support dative bonds, and in later versions there are
    significant changes in the embedding and symmetry processing algorithms
    which are not well compatible with the **epic-mace**'s underlying algorithms.
    This noticeably increases number of errors for both stereomer search and 3D embedding.


conda
-----

We highly recommend to install MACE via the **conda** package management system (`download`_, `install`_),
since it's the easiest way to get RDKit 2020.09. The following commands will create a new conda environment
with Python 3.7, RDKit 2020.09, and the latest version of **epic-mace**: ::

    > conda create -n mace -c rdkit python=3.7 rdkit=2020.09
    > conda install -n mace -c grimgenius epic-mace

Do not forget to activate the environment before using **epic-mace**: ::
    
    > conda activate mace


pip
---

**epic-mace** can be installed via **pip**: ::

    > pip install rdkit
    > pip install epic-mace

However, we strongly recommend installation via conda, since the earliest available RDKit version on PyPI is 2022.03 which does not ensure the stable operation of the **epic-mace** package.
Though it is enough for demonstrational purposes or automatic documentation generation.

In extreme cases, one can install **epic-mace** via pip to the conda environment with preinstalled RDKit 2020.09: ::

    > conda create -n mace -c rdkit python=3.7 rdkit=2020.09.1
    > conda activate mace
    > pip install epic-mace

.. warning::
    Please note, that PyPI epic-mace package does not contain rdkit in the requirements list to avoid possible conflicts between conda and pip RDKit installations.
    Therefore, RDKit must be installed beforehand.


.. _RDKit: https://www.rdkit.org/
.. _download: https://www.anaconda.com/products/distribution
.. _install: https://docs.anaconda.com/anaconda/install/


