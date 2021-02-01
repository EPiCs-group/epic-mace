# MACE: MetAl Complexes Embedding

Python library (*current version 0.4.0*) and command-line tool (*under development*) for generation of 3D coordinates for complexes of d-/f-elements.

## Installation

MACE can be easily installed with conda package manager ([ref](https://anaconda.org/epics.group/epic-mace)):

```
> conda create -n mace
> conda activate mace
> conda install -c epics.group epic-mace
```

You can also use pip ([ref](https://pypi.org/project/epic-mace/)), but in this case [RDKit](https://www.rdkit.org/) 2020.03 or later must be already installed, e.g. using conda:

```
> conda create -n mace
> conda activate mace
> conda install -c rdkit rdkit=2020.09
> conda install pip
> pip install epic-mace
```

Such scenario may be useful if there are some problems with compatibility of libraries during direct conda installation.

## CookBook

This section briefly describes the main possibilities of the package. For more details and usage cases see the [manual](manual/manual.md) (under development).

### Generate 3D coordinates for the complex

First, you need to initialize the complex. The easiest way to do it is to draw the complex in [Marvin Sketch](https://chemaxon.com/products/marvin) and copy its [ChemAxon SMILES](https://docs.chemaxon.com/display/docs/smiles.md):

<img src="manual/pics/README/marvin_copy_smiles.png" width="50%" />

Bonds between donor atoms and the central ion must be encoded as dative bonds. Atomic map numbers (blue numbers close to donor atoms) are used to describe the spatial arrangement of ligands.

After copying SMILES, initialize the Complex object:

```python
import mace

# copied SMILES
smiles = '[Ru++]12([H-:3])([C-:5]#[O+])([C-:6]#[O+])[NH:1](CC[P:4]1(C)C)CC[P:2]2(C)C |C:9.9,14.15,6.5,1.0,2.1,4.3|'
geom = 'OH' # OH - octahedral, SP - square planar
X = mace.Complex(smiles, geom)
X.mol
```

<img src="manual/pics/README/X1.png" width="40%" />

To generate atomic coordinates, use `AddConformer` or `AddConformers` methods:

```python
conf_idxs = X.AddConformers(numConfs = 5)
X.ToXYZ('test.xyz', confId = 'all')
```

Here's the result:

<img src="manual/pics/README/X1_3D.png" width="50%"/>

### Search of stereomers

In homogeneous catalysis you do not know structure of the complex *a priori*. Thus, you need to analyze all possible stereomers. For the purpose, first initialize complex, and this time let's do it from ligands:

```python
import mace

# SMILES of ligands copied from Marvin Sketch
# donor atoms must have any non-zero map number
ligands = ['C[P:1](C)CC[NH:1]CC[P:1](C)C', '[H-:1]', '[C-:1]#[O+]', '[C-:1]#[O+]']
CA = '[Ru+2]' # SMILES of central atom
geom = 'OH'
X = mace.ComplexFromLigands(ligands, CA, geom)
```

You can not generate atomic coordinates for this complex as it contains incorrect information on the spatial arrangement of ligands (all donor atoms are in the first position):

```python
X.AddConformer()
##> Bad SMILES: isotopic labels are not unique
##> 
##> The initial SMILES contains insufficient or erroneous info
##> on the positions of the ligands around the central atom
##> encoded with isotopic labels.
##> To use 3D generation and other features, generate
##> possible stereomers using GetStereomers method.
print(X.GetNumConformers())
##> 0
```

Next, let's find all possible stereomers:

```python
Xs = X.GetStereomers(regime = 'all', dropEnantiomers = False)
print(len(Xs))
##> 9
Xs = X.GetStereomers(regime = 'all', dropEnantiomers = True)
print(len(Xs))
##> 7
```

So, this complex has 7 stereomers, and 2 of them are enantiomeric. Let's generate atomic coordinates for them and save to XYZ-files:

```python
for i, X in enumerate(Xs):
    X.AddConformers(numConfs = 5)
    X.ToXYZ(f'X1_{i}.xyz', confId = 'min')
```

Here's the result (aliphatic hydrogens removed for clarity):

<img src="manual/pics/README/X1s_3D.png" width="100%"/>

### Introduction of substituents

Imagine, that you need to generate a lot of complexes with the same core and different substituents. Generation of all structures will result in serious time loss due to the QM geometry optimization. There are more tricky way. First, you need to draw ligands containing dummy atoms (or hydrogens) with isotopic labels equal to the number of substituent (can be encoded as *R1*, *R2*, etc. in Marvin Sketch):

<img src="manual/pics/README/subs.png" width="50%"/>

After that, you need to generate atomic coordinates for the core structure:

```python
import mace

# core
ligands = ['C[P:3](C)CC1=CC([*])=C([*])C(C[P:1](C)C)=[N:2]1 |$;;;;;;;_R1;;_R2;;;;;;$,c:14,t:4,7|', '[Cl-:4]']
CA = '[Rh+]'
geom = 'SP'
core = mace.ComplexFromLigands(ligands, CA, geom)

# 3D embedding
core.AddConformer()
core.ToXYZ('core.xyz', confId = 'min')
```

Next, we need to substitute dummies using the `mace.AddSubsToMol` function:

```python
# substituents
subs = {'R1': mace.MolFromSmiles('[*]OC'),
        'R2': mace.MolFromSmiles('[*]C#N')}
# add subs
X = mace.ComplexFromMol(mace.AddSubsToMol(core.mol, subs), core._geom)
```

And generate atomic coordinates using constrained embedding using `AddConstrainedConformer` method:

```python
X.AddConstrainedConformer(core, confId = core.GetMinEnergyConfId(0))
X.ToXYZ('X.xyz', confId = 'min')
```

`core` must be a substructure of `X`, and the number of donor atoms in `core` and `X` must be equal. In that case we'll obtain atomic coordinates, and the positions of matched atoms will be almost the same as those of a template molecule:

<img src="manual/pics/README/X_and_core.png" width="80%"/>

This method can be easily used to generate a big set of related complexes, which may be useful for QSAR research.

