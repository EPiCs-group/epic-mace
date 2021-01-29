'''
Generates
'''

#%% Imports

import os
from itertools import product

from rdkit import Chem # remove after adding MolToSmiles to mace
import mace


#%% Define complexes

# ligands
ligands_NH = {'PNP': '[*][P:1]([*])CC[NH:1]CC[P:1]([*])[*] |$_R1;;_R1;;;;;;;_R1;_R1$|',
              'CNPhP': 'CC1=CC(C)=C(N2[C:1]N(CC[NH:1]CC3=C(C=CC=C3)[P:1]([*])[*])C=C2)C(C)=C1 |$;;;;;;;;;;;;;;;;;;;;_R1;_R1;;;;;$,c:15,17,23,27,t:1,4,13,^3:7|',
              'CNzCyP': 'CC1=CC(C)=C(N2[C:1]N(CC[NH:1]C[C@H]3CCCC[C@H]3[P:1]([*])[*])C=C2)C(C)=C1 |r,$;;;;;;;;;;;;;;;;;;;;_R1;_R1;;;;;$,c:23,27,t:1,4,^3:7|',
              'CNeCyP': 'CC1=CC(C)=C(N2[C:1]N(CC[NH:1]C[C@@H]3CCCC[C@H]3[P:1]([*])[*])C=C2)C(C)=C1 |r,$;;;;;;;;;;;;;;;;;;;;_R1;_R1;;;;;$,c:23,27,t:1,4,^3:7|',
              'PNpyP': '[*][P:1]([*])CC1=CC=CC(C[P:1]([*])[*])=[N:1]1 |$_R1;;_R1;;;;;;;;;_R1;_R1;$,c:6,12,t:4|'}
ligands_N = {'PNP': '[*][P:1]([*])CC[N-:1]CC[P:1]([*])[*] |$_R1;;_R1;;;;;;;_R1;_R1$|',
             'CNPhP': 'CC1=CC(C)=C(N2[C:1]N(CC[N-:1]CC3=C(C=CC=C3)[P:1]([*])[*])C=C2)C(C)=C1 |$;;;;;;;;;;;;;;;;;;;;_R1;_R1;;;;;$,c:15,17,23,27,t:1,4,13,^3:7|',
             'CNzCyP': 'CC1=CC(C)=C(N2[C:1]N(CC[N-:1]C[C@H]3CCCC[C@H]3[P:1]([*])[*])C=C2)C(C)=C1 |r,$;;;;;;;;;;;;;;;;;;;;_R1;_R1;;;;;$,c:23,27,t:1,4,^3:7|',
             'CNeCyP': 'CC1=CC(C)=C(N2[C:1]N(CC[N-:1]C[C@@H]3CCCC[C@H]3[P:1]([*])[*])C=C2)C(C)=C1 |r,$;;;;;;;;;;;;;;;;;;;;_R1;_R1;;;;;$,c:23,27,t:1,4,^3:7|',
             'PNpyP': '[*][P:1]([*])CC1=CC=C\C([N-:1]1)=C\[P:1]([*])[*] |$_R1;;_R1;;;;;;;;;;_R1;_R1$,c:6,t:4|'}

# subs
Rs = {'Et' : {'R1': '[*]CC'},
      'iPr': {'R1': '[*]C(C)C'},
      'Cy' : {'R1': '[*]C1CCCCC1'},
      'Ph' : {'R1': '[*]c1ccccc1'}}

# CA and aux ligands
CA = '[Mn+]'
auxs_NH = {'MnH(CO)2': ['[H-:1]', '[C-:1]#[O+]', '[C-:1]#[O+]']}
auxs_N = {'Mn(CO)2': ['[C-:1]#[O+]', '[C-:1]#[O+]'],
          'Mn(CO)3': ['[C-:1]#[O+]', '[C-:1]#[O+]', '[C-:1]#[O+]']}

# other params
geom = 'OH'
merRule = True
numConfs = 3
confId = -1 # save all confs

# output path
path_out = 'xyz'


#%% Generate NHs

for (lig_name, lig_smiles), (subs_name, subs) in product(ligands_NH.items(), Rs.items()):
    # modify ligand
    mol = mace.MolFromSmiles(lig_smiles)
    mol = mace.AddSubsToMol(mol, subs)
    ligands = [Chem.MolToSmiles(mol)]
    # add auxilaries
    for aux_name, aux in auxs_NH.items():
        # output
        basename = f'{lig_name}_{subs_name}_{aux_name}'
        if not os.path.exists(f'{path_out}/{basename}'):
            os.mkdir(f'{path_out}/{basename}')
        # get stereomers
        X = mace.ComplexFromLigands(ligands + aux, CA, geom)
        Xs = X.GetStereomers(merRule = merRule)
        # make 3D
        absent = []
        for i, Xi in enumerate(Xs):
            cis = Xi.AddConformers(numConfs = numConfs)
            if not cis:
                absent.append(Chem.MolToSmiles(Xi.mol))
                continue
            Xi.ToXYZ(f'{path_out}/{basename}/{basename}_{i}.xyz', confId = confId)
        if absent:
            with open(f'{path_out}/{basename}/absent.txt', 'w') as outf:
                outf.write('\n'.join(absent)+'\n')



#%% Generate N-

for (lig_name, lig_smiles), (subs_name, subs) in product(ligands_N.items(), Rs.items()):
    # modify ligand
    mol = mace.MolFromSmiles(lig_smiles)
    mol = mace.AddSubsToMol(mol, subs)
    ligands = [Chem.MolToSmiles(mol)]
    # add auxilaries
    for aux_name, aux in auxs_N.items():
        # output
        basename = f'{lig_name}_{subs_name}_anion_{aux_name}'
        if not os.path.exists(f'{path_out}/{basename}'):
            os.mkdir(f'{path_out}/{basename}')
        # get stereomers
        X = mace.ComplexFromLigands(ligands + aux, CA, geom)
        Xs = X.GetStereomers(merRule = merRule)
        # make 3D
        absent = []
        for i, Xi in enumerate(Xs):
            cis = Xi.AddConformers(numConfs = numConfs)
            if not cis:
                absent.append(Chem.MolToSmiles(Xi.mol))
                continue
            Xi.ToXYZ(f'{path_out}/{basename}/{basename}_{i}.xyz', confId = confId)
        if absent:
            with open(f'{path_out}/{basename}/absent.txt', 'w') as outf:
                outf.write('\n'.join(absent)+'\n')


