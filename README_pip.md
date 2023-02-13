# MACE: MetAl Complexes Embedding

Python library (*current version 0.4.0*) and command-line tool (under development) for generation of 3D coordinates for complexes of d-/f-elements.

For more details including usage examples see the [GitHub page](https://github.com/EPiCs-group/mace).

## Installation

This project uses RDKit as an 3D embedding engine, so its performance depends on RDKit build. We found that conda's 2020.09.1 version gives the lowest error rate. Therefore, if you are unhappy with the standard pip installation, try to install RDKit before the mace package:

```
> conda create -n mace
> conda activate mace
> conda install -c rdkit rdkit=2020.09.1
> conda install pip # if required
> pip install epic-mace
```

Enjoy!


