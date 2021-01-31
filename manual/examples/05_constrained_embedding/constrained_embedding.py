'''
Direct 3D generation of all complexes
'''

import mace

#%% Prepare cores

# main ligand
ligand_smiles = '[*][P:1]([*])[C@@H]1CCCC[C@@H]1C[N-:1]CCN1C=C[N+]([*])=[C-:1]1 |r,$_R1;;_R1;;;;;;;;;;;;;;;_R2;$,c:15,18|'
ligand_smiles = '[1H][P:1]([1H])[C@@H]1CCCC[C@@H]1C[N-:1]CCN1C=C[N+]([2H])=[C-:1]1'

# core complex
ligands = [ligand_smiles, '[C-:1]#[O+]', '[C-:1]#[O+]', '[C-:1]#[O+]']
CA = '[Mn+]'
geom = 'OH'
X = mace.ComplexFromLigands(ligands, CA, geom)

# stereomers
Xs = X.GetStereomers()

# generate conformers
for i, X in enumerate(Xs):
    X.AddConformers(numConfs = 20, rmsThresh = 1.0)
    X.ToXYZ(f'xyz/X{i}.xyz', 'all')


#%% Constrained embedding

# subs
Rs = {'R1': mace.MolFromSmiles('[*]CC'),
      'R2': mace.MolFromSmiles('[*]c1c(C)cc(C)cc(C)1')}
Rs = {'R1': mace.MolFromSmiles('[*]C'),
      'R2': mace.MolFromSmiles('[*]C')}

# complexes
X0 = Xs[0]
X1 = mace.ComplexFromMol(mace.AddSubsToMol(X0.mol, Rs), 'OH')
X1.AddConformer()

# embedding
confId = 2 # check XYZ file for the best one
X1.AddConstrainedConformer(X0, confId)


