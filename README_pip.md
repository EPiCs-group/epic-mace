# MACE: MetAl Complexes Embedding

MACE is an open source toolkit for the automated screening and discovery of octahedral and square-planar mononuclear complexes. MACE is developed by the [Evgeny Pidko Group](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/chemical-engineering/principal-scientists/evgeny-pidko/evgeny-pidko-group) in the [Department of Chemical Engineering](http://web.mit.edu/cheme/) at [TU Delft](https://www.tudelft.nl/en/). The software can generate all possible configurations for square-planar and octahedral metal complexes and atomic 3D coordinates suitable for quantum-chemical computations. It supports ligands of high complexity and can be used for the development of a massive computational pipelines aimed at solving problems of homogenious catalysis.

For more details see the [GitHub page](https://github.com/EPiCs-group/mace).

## Installation

This package uses RDKit as an 3D embedding engine, so its performance depends on RDKit build. We found that conda's 2020.09.1 version gives the lowest error rate. Therefore, if you are unhappy with the standard pip installation, try to install RDKit before the mace package:

```bash
> conda create -n mace
> conda activate mace
> conda install -c rdkit rdkit=2020.09.1
> conda install pip # if required
> pip install epic-mace
```

Enjoy!
