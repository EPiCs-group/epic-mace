# MACE: MetAl Complexes Embedding

MACE is an open source toolkit for the automated screening and discovery of metal complexes. MACE is developed by the [Evgeny Pidko Group](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group) in the [Department of Chemical Engineering](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/) at [TU Delft](https://www.tudelft.nl/en/). The main features of the MACE package are to discover all possible configurations for square-planar and octahedral metal complexes, and generate atomic 3D coordinates suitable for quantum-chemical computations. MACE shows high performance (> 99% success rate) for complexes of ligands of high denticity (up to 6), and thus is well-suited for the development of a massive computational pipelines aimed at solving problems of homogeneous catalysis.

## Requirements

- Python 3.7 or higher (Python 3.7 is recommended);

- [RDKit](https://www.rdkit.org/) 2020.09 or higher (RDKit 2020.09 is a **must** for a correct functioning).


## Installation

### conda

We highly recommend to install MACE via the [conda](https://anaconda.org/grimgenius/epic-mace) package management system. The following commands will create a new conda environment with Python 3.7, RDKit 2020.09, and the latest version of MACE:

```ssh
> conda create -n mace -c rdkit python=3.7 rdkit=2020.09
> conda install -n mace -c grimgenius epic-mace
```

The reason for the strong preference for installation via conda is that only the RDKit 2020.09 version ensures correct and error-free operation of the MACE package. Earlier versions do not support dative bonds, and in later versions there are significant changes in the embedding and symmetry processing algorithms which are not well compatible with the MACE's underlying algorithms. This noticeably increases number of errors for both stereomer search and 3D embedding.

### pip

MACE can be installed via ([pip](https://pypi.org/project/epic-mace/)):

```bash
> pip install rdkit
> pip install epic-mace
```

However, we strongly recommend installation via conda, since the earliest available RDKit version on PyPI is 2022.03 which does not ensure the stable operation of the MACE package. Though it is enough for demonstrational purposes or automatic documentation generation.

In extreme cases, one can install MACE via pip to the conda environment with preinstalled RDKit 2020.09:

```bash
> conda create -n mace python=3.7 rdkit=2020.09.1 -c rdkit
> conda activate mace
> pip install epic-mace
```

Please note, that PyPI epic-mace package does not contain rdkit in the requirements list to avoid possible conflicts between conda and pip RDKit installations. Therefore, you must install RDKit manually beforehand.

## Performance

MACE shows high performance (> 99% success rate) for complexes of ligands, extracted from Cambridge Structural Database. For more details see [performance](https://github.com/EPiCs-group/epic-mace/blob/master/performance/README.ipynb).

## GUI

For convenient interactive research of metal complexes, as well as for a better understanding of MACE features, one can use [web applications](https://github.com/IvanChernyshov/mace-notebooks) built on IPython notebooks.

