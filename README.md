# MACE: MetAl Complexes Embedding

MACE is an open source toolkit for the automated screening and discovery of metal complexes. MACE is developed by [Ivan Chernyshov](https://github.com/IvanChernyshov) as part of the [Evgeny Pidko Group](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group) in the [Department of Chemical Engineering](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/) at [TU Delft](https://www.tudelft.nl/en/). The main features of the MACE package are to discover all possible configurations for square-planar and octahedral metal complexes, and generate atomic 3D coordinates suitable for quantum-chemical computations. MACE shows high performance for complexes of ligands of high denticity (up to 6), and thus is well-suited for the development of a massive computational pipelines aimed at solving problems of homogeneous catalysis.


## Installation

> **epic-mace** requires Python 3.7 and the 2020.09 version of RDKit for a correct functioning!

Earlier versions do not support dative bonds, and in later versions there are significant changes in the embedding and symmetry processing algorithms which are not well compatible with the **epic-mace**’s underlying algorithms. This noticeably increases number of errors for both stereomer search and 3D embedding.

### conda

We highly recommend to install MACE via the [conda](https://anaconda.org/grimgenius/epic-mace) package management system. The following commands will create a new conda environment with Python 3.7, RDKit 2020.09, and the latest version of **epic-mace**:

```bash
> conda create -n mace -c rdkit python=3.7 rdkit=2020.09
> conda install -n mace -c grimgenius epic-mace
```

Do not forget to activate the environment before using epic-mace:

```bash
> conda activate mace
```

### pip

**epic-mace** can be installed via ([pip](https://pypi.org/project/epic-mace/)):

```bash
> pip install rdkit
> pip install epic-mace
```

However, we strongly recommend installation via conda, since the earliest available RDKit version on PyPI is 2022.03 which does not ensure the stable operation of **epic-mace**. Though it is enough for demonstrational purposes or automatic documentation generation.

In extreme cases, one can install MACE via pip to the conda environment with preinstalled RDKit 2020.09:

```bash
> conda create -n mace python=3.7 rdkit=2020.09.1 -c rdkit
> conda activate mace
> pip install epic-mace
```

Please note, that PyPI epic-mace package does not contain rdkit in the requirements list to avoid possible conflicts between conda and pip RDKit installations. Therefore, you must install RDKit manually beforehand.


## Main features

1. Stereomer search for octahedral and square-planar complexes.

2. Generation of 3D atomic coordinates, including instruments for conformer sampling.

3. Modification of ligands with predefined substituents.

4. Generation of geometry of coordinated ligands for [molSimplify](https://molsimplify.mit.edu/).

5. Two available interfaces:
    
    - command-line interface for routine tasks;
    
    - Python package for organizing complex computational pipelines.


## Useful links

1. [Documentation](https://epic-mace.readthedocs.io/en/latest/)

2. [Performance](https://github.com/EPiCs-group/epic-mace/blob/master/performance/README.ipynb)

3. [CLI examples](https://github.com/EPiCs-group/epic-mace/tree/master/examples)

