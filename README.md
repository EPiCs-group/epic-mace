# MACE: MetAl Complexes Embedding

MACE is an open source toolkit for the automated screening and discovery of octahedral and square-planar mononuclear complexes. MACE is developed by the [Evgeny Pidko Group](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group) in the [Department of Chemical Engineering](http://web.mit.edu/cheme/) at [TU Delft](https://www.tudelft.nl/en/). The software generates all possible configurations for square-planar and octahedral metal complexes and atomic 3D coordinates suitable for quantum-chemical computations. It supports ligands of high complexity and can be used for the development of a massive computational pipelines aimed at solving problems of homogenious catalysis.

### Installation

MACE can be installed via pip ([ref](https://pypi.org/project/epic-mace/)):

```bash
> pip install epic-mace
```

Since this package uses RDKit as an 3D embedding engine, its performance depends on RDKit build. We found that conda's 2020.09.1 version gives the lowest error rate. Therefore, if you are unhappy with the standard pip installation, try installing RDKit before the mace package:

```bash
> conda create -n mace
> conda activate mace
> conda install -c rdkit rdkit=2020.09.1
> conda install pip
> pip install epic-mace
```

### Tutorials

- [MACE package](tutorials/mace.ipynb): describes MACE functionality for stereoisomer searching and 3D embedding.

### GUI

For convenient interactive research of metal complexes, as well as for a better understanding of MACE features, one can use [web applications](https://github.com/IvanChernyshov/mace-notebooks) built on IPython notebooks.

### Performance

MACE shows high performance (> 99% success rate) for complexes of ligands, extracted from Cambridge Structural Database. For more details see [performance](performance).
