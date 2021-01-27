'''
MACE demonstration - not needed, rework to README
'''

#%% Imports

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace


#%% Paths

# root
if not os.path.isdir('outputs/demo'):
    os.mkdir('outputs/demo')
path_dir = 'outputs/demo'


#%% Stereo

smiles = '[Co+3](N)(N)(N)([F-])([F-])[F-] |C:4.3,5.4,1.0,6.5,2.1,3.2|'
geom = 'OH'
X = mace.Complex(smiles, geom)
Xs = X.GetStereomers()
for i, x in enumerate(Xs):
    x.AddConformer()
    x.ToXYZ(f'{path_dir}/demo_stereo_{i}.xyz')


#%% 3D from ChemAxon Marvin SMILES

smiles = '[Co+3]([NH3:4])([NH3:5])([NH3:3])([F-:1])([F-:6])[F-:2] |C:4.3,5.4,1.0,6.5,2.1,3.2|'
geom = 'OH'
X = mace.Complex(smiles, geom)
X.AddConformer()
X.ToXYZ(f'{path_dir}/demo_chem_axon.xyz')


#%% 3D from ligands

ligands = ['CN1C=C[N+](CC[NH:3]CC[P:2](C)C)=[C-:4]1 |c:2,12|',
           '[C-:6]#[O+]', '[H-:1]', '[H-:5]']
CA = '[Ru+2]'
geom = 'OH'
smiles = '[Ru++]12([C-:4]3=[N+](CC[NH:3]1CC[P:2]2(C)C)C=CN3C)([C-:6]#[O+])([H-:1])[H-:5] |c:1,13,C:8.9,5.5,1.0,17.19,15.17,18.20|'
#X = mace.ComplexFromLigands(ligands, CA, geom)
X = mace.Complex(smiles, geom)
X.AddConformers(numConfs = 5, rmsThresh = 2.0)
X.ToXYZ(f'{path_dir}/demo_from_ligands.xyz')


#%% Substituents support

smiles = '[*]C1=C([*])N([C-:4]2=[N+]1CC[NH:3]1[Ru++]2([C-:6]#[O+])([P:2](CC1)(C)C)([H-:1])[H-:5])C |$_R1;;;_R2;;;;;;;;;;;;;;;;;$,c:1,5,C:13.14,9.10,5.11,18.20,11.12,19.21|'
mol = mace.MolFromSmiles(smiles)
subs = {'R1': '[*]C', 'R2': '[*]N(=O)=O'}
mol = mace.AddSubsToMol(mol, subs)
geom = 'OH'
X = mace.ComplexFromMol(mol, geom)
X.AddConformers(numConfs = 5, rmsThresh = 2.0)
X.ToXYZ(f'{path_dir}/demo_subs.xyz')


#%% Constrained embedding

smiles = '[Co+3]([NH3:4])([NH3:5])([NH3:3])([F-:1])([F-:6])[F-:2] |C:4.3,5.4,1.0,6.5,2.1,3.2|'
geom = 'OH'
X = mace.Complex(smiles, geom)
X.AddConformer()
X.ToXYZ(f'{path_dir}/demo_constrained_dummy.xyz')

smiles = '[N:5]([Co+3]([NH3:4])([NH3:3])([F-:1])([F-:6])[F-:2])(C)(C)C |C:2.1,0.0,3.2,4.3,5.4,6.5|'
geom = 'OH'
X = mace.Complex(smiles, geom)
X.AddConstrainedConformerFromXYZ(pathCore = f'{path_dir}/x.xyz', ignoreHs = False)
X.ToXYZ(f'{path_dir}/demo_constrained_final.xyz')


