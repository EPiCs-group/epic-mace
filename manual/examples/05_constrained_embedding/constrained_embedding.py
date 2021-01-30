'''
Direct 3D generation of all complexes
'''

import mace

#%% Prepare cores

# main ligand
ligand_smiles = '[*][P:1]([*])[C@@H]1CCCC[C@@H]1C[N-:1]CCN1C=C[N+]([*])=[C-:1]1 |r,$_R1;;_R1;;;;;;;;;;;;;;;_R2;$,c:15,18|'

# subs
Rs = {'R1': mace.MolFromSmiles('[*]CC'),
      'R2': mace.MolFromSmiles('[*]c1c(C)cc(C)cc(C)1')}

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
    X.ToXYZ(f'X{i}.xyz', 'all')


#%% Constrained embedding




