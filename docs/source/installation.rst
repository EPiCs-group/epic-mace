Installation
============

conda
-----

We highly recommend to install MACE via the `conda`_ package management system.
The following commands will create a new conda environment with Python 3.7, RDKit 2020.09,
and the latest version of MACE: ::

    > conda create -n mace -c rdkit python=3.7 rdkit=2020.09
    > conda install -n mace -c grimgenius epic-mace

The reason for the strong preference for installation via conda is that only the RDKit 2020.09 version ensures
correct and error-free operation of the MACE package. Earlier versions do not support dative bonds,
and in later versions there are significant changes in the embedding and symmetry processing algorithms
which are not well compatible with the MACE's underlying algorithms. This noticeably increases number of errors
for both stereomer search and 3D embedding.

pip
---

MACE can be installed via `pip`_: ::

    > pip install rdkit
    > pip install epic-mace

However, we strongly recommend installation via conda, since the earliest available RDKit version on PyPI is 2022.03 which does not ensure the stable operation of the MACE package.
Though it is enough for demonstrational purposes or automatic documentation generation.

In extreme cases, one can install MACE via pip to the conda environment with preinstalled RDKit 2020.09: ::

    > conda create -n mace -c rdkit python=3.7 rdkit=2020.09.1
    > conda activate mace
    > pip install epic-mace

Please note, that PyPI epic-mace package does not contain rdkit in the requirements list to avoid possible conflicts between conda and pip RDKit installations.
Therefore, you must install RDKit manually beforehand.


.. _conda: https://anaconda.org/grimgenius/epic-mace
.. _pip: https://pypi.org/project/epic-mace/
