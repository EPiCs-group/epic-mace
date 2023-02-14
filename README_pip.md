# MACE: MetAl Complexes Embedding

MACE is an open source toolkit for the automated screening and discovery of octahedral and square-planar mononuclear complexes. MACE is developed by the [Evgeny Pidko Group](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group) in the [Department of Chemical Engineering](http://web.mit.edu/cheme/) at [TU Delft](https://www.tudelft.nl/en/). The software generates all possible configurations for square-planar and octahedral metal complexes and atomic 3D coordinates suitable for quantum-chemical computations. It supports ligands of high complexity and can be used for the development of a massive computational pipelines aimed at solving problems of homogenious catalysis.

For more details see the [GitHub page](https://github.com/EPiCs-group/mace).

## Installation

### conda

We highly recommend to install MACE via the [conda](https://conda.io/docs/) package management system. The following command will create new conda environment with Python 3.7, RDKit 2020.09, and the latest version of MACE:

```ssh
> conda create -n mace epic-mace -c grimgenius
```

The reason for the strong preference for installation via conda is that only the RDKit 2020.09 version ensures failure- and error-free operation of the MACE package. Earlier versions do not support dative bonds, and in later versions there are significant changes in the embedding and symmetry processing algorithms which are not well compatible with the MACE's underlying algorithms.

### pip

MACE can be installed via pip ([ref](https://pypi.org/project/epic-mace/)):

```bash
> pip install epic-mace
```

However, we strongly recommend installation via conda, since the earliest available RDKit version on PyPI is 2022.03 which does not ensure the stable operation of the MACE package.

In extreme cases, one can install MACE via pip to the conda environment with preinstalled RDKit 2020.09:

```bash
> conda create -n mace python=3.7 rdkit=2020.09.1 -c rdkit
> conda activate mace
> pip install epic-mace
```

Please note, that setup.py does not contain rdkit in the requirements list to avoid possible conflicts between conda and pip RDKit installations.
