'''
Tests MACE code
'''

import re, sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import mace

from rdkit import Chem


#%% Support functions

def MolToMolFile(mol, path = 'C:/Users/ivan-/Desktop/x.mol', confId = 0):
    '''
    Substitutes 'R' to 'X' for ChemCraft visualization
    '''
    text = Chem.MolToMolBlock(mol, confId = confId)
    text = re.sub(' R ', ' X ', text)
    with open(path, 'w') as outf:
        outf.write(text)



#%% Geometry tests

# 1


# 2


# 3


# 4





#%% Raw

# full 3D

# Zhenya's ligand
#ligands = ['c1ccccc1[P:2](c1ccccc1)CC[NH:3]CCN1C=CN(c2c(C)cc(C)cc(C)2)[C:4]1', '[O+]#[C-:1]', '[H-:5]', '[H-:6]']
ligands = ['c1ccccc1[P:2](c1ccccc1)CC[NH:3]CC[N+]=1C=CN(c2c(C)cc(C)cc(C)2)[C-:4]1', '[O+]#[C-:1]', '[H-:5]', '[H-:6]']
ligands = ['C[P:2](C)CC[NH:3]CC[N+]=1C=CN(C)[C-:4]1', '[O+]#[C-:1]', '[O+]#[C-:5]', '[O+]#[C-:6]']
#ligands = ['c1ccccc1[P:2](c1ccccc1)CC[NH:3]CC[N+]=1C=CN(c2c(C)cc(C)cc(C)2)[C-:4]1']
#ligands = ['c1ccccc1[P:2](c1ccccc1)CC[NH:3]CCc1cccc(c2c(C)cc(C)cc(C)2)[n:4]1', '[O+]#[C-:1]', '[H-:5]', '[H-:6]']
#ligands = ['c1ccccc1[P:2](c1ccccc1)CC[NH:3]CCc1cccc(c2c(C)cc(C)cc(C)2)[c-:4]1', '[O+]#[C-:1]', '[H-:5]', '[H-:6]']
#ligands = ['c1ccccc1[P:2](c1ccccc1)CC[NH:3]CCC1=CC=C(c2c(C)cc(C)cc(C)2)[N-:4]1', '[O+]#[C-:1]', '[H-:5]', '[H-:6]']
#ligands = ['c1ccccc1[P:2](c1ccccc1)CC[NH:3]CCC1=CN(C)C(c2c(C)cc(C)cc(C)2)=[N:4]1', '[O+]#[C-:1]', '[H-:5]', '[H-:6]']
ligands = ['c1ccccc1[C-:4]=CC=[C-:5]c1ccccc1', '[O+]#[C-:1]', '[O+]#[C-:2]', '[O+]#[C-:3]', '[O+]#[C-:6]']
#ligands = ['C[P:2](C)CCN1C=C[N+](CC[P:4](C)C)=[C-:3]1', '[O+]#[C-:1]', '[O+]#[C-:5]', '[O+]#[C-:6]']
CA = '[*+2]'
geom = 'OH'
regime = 'all'
drop_enantiomers = True
X = mace.ComplexFromLigands(ligands, CA, geom)
Xs = X.GetStereomers(regime = regime, drop_enantiomers = drop_enantiomers)
print(len(Xs))

# # 3D
# for i, x in enumerate(Xs):
#     x.AddConformers(n = 3, rms_thresh = 2.0)
#     x.ToXYZ(f'C:/Users/ivan-/Desktop/x_{i}.xyz', -2)

# # mol3Dx coords
# idx = 2
# x0 = Xs[idx]
# #x0.AddConformers(n = 3, rms_thresh = 2.0)
# for _ in range(3):
#     x0.AddConformer(clear_confs = False)
# for i in range(x0.mol3Dx.GetNumConformers()):
#     MolToMolFile(x0.mol3Dx, f'C:/Users/ivan-/Desktop/xx_{i}_init.mol', confId = i)
# x0.ToXYZ(f'C:/Users/ivan-/Desktop/x_{idx}.xyz')
# x1 = ComplexFromXYZ(f'C:/Users/ivan-/Desktop/x_{idx}.xyz')
# for i in range(x1.mol3Dx.GetNumConformers()):
#     MolToMolFile(x1.mol3Dx, f'C:/Users/ivan-/Desktop/xx_{i}_file.mol', confId = i)

# # constrained embedding

# # core partial
# ligands = ['c1ccccc1[P:2](c1ccccc1)CC[N@H:3](->[*])CC[N+]=1C=CN(c2c(C)cc(C)cc(C)2)[C-:4]1']
# #ligands = ['CC#[N:1]','CC#[N:2]','CC#[N:3]']
# CA = '[Mn]'
# geom = 'OH'
# core = ComplexFromLigands(ligands, CA, geom)
# core.AddConformer()
# core.ToXYZ('C:/Users/ivan-/Desktop/core.xyz')
# path_core = 'C:/Users/ivan-/Desktop/core.xyz'
# # modified
# ligands = ['c1ccccc1[P:2](c1ccccc1)CC[N@H:3](->[*])CC[N+]=1C(C#N)=C(C#CC(C)(C)C)N(c2c(C)cc(C)cc(C)2)[C-:4]1', '[O+]#[C-:1]',  '[O+]#[C-:5]', '[O+]#[C-:6]']
# #ligands = ['CC#[N:1]','CC#[N:2]','CC#[N:3]', '[O+]#[C-:4]',  '[O+]#[C-:5]', '[O+]#[C-:6]']
# CA = '[Mn]'
# geom = 'OH'
# X = ComplexFromLigands(ligands, CA, geom)
# # # constrained embedding
# ignoreHs = True
# clear_confs = True
# #flag = X.AddConstrainedConformer(core, confId = 0, ignoreHs = ignoreHs, clear_confs = clear_confs)
# flag = X.AddConstrainedConformerFromXYZ(path_core, ignoreHs = ignoreHs, clear_confs = clear_confs)
# print(flag)
# if flag != -1:
#     X.ToXYZ('C:/Users/ivan-/Desktop/full.xyz')



# smiles = 'C1=CC=CC(=C1)N1C=C[N+]2=[C-:2]1[Mn+]1([C-:5]#[O+])([C-:6]#[O+])([P:4](CC[N@@:3]1(CC2)[H])(C)C)[H-:1] |r,c:0,2,4,8,10,C:19.21,16.17,10.12,25.28,12.13,14.15|'
# smiles = 'C1=CC=CC(=C1)N1C=C[N+]2=[C-:2]1[Mn+]1([C-:5]#[O+])([C-:6]#[O+])([P:4](CC[N@:3]1(CC2)[H])(C)C)[H-:1] |r,c:0,2,4,8,10,C:19.21,16.17,10.12,25.28,12.13,14.15|'
# #smiles = '[Mn+]12([H-:1])([C-:5]#[O+])([C-:6]#[O+])[C-:2]3=[N+](C[C@H](C)[N@:3]1([H])CC[P:4]2(C)CC)C=CN3C1=CC=CC=C1 |r,c:6,21,27,29,t:25,C:11.11,15.16,6.5,1.0,2.1,4.3|'
# #smiles = 'N1(C2=CC=CC=C2)C=C[N+]2=[C-:2]1[Mn+]1([N@:3]([C@H](C2)C)(CC[P:4]1(C)C)[H])([C-:5]#[O+])([C-:6]#[O+])[H-:1] |r,c:3,5,8,10,t:1,C:12.13,18.21,10.12,26.29,22.25,24.27|'
# X = Complex(smiles, 'OH')
# X.mol
# ps = Chem.SmilesParserParams()
# ps.sanitize = False
# ps.removeHs = False
# Chem.MolFromSmiles(smiles.split()[0], ps) # inversion happens

# X.AddConformer()
# X.ToXYZ('C:/Users/ivan-/Desktop/x.xyz')



# # 1
# ligands = ['[H-:1]->*','*<-[C-:5]#[O+]','*<-[C-:6]#[O+]','*<-[C-:2]1=[N+](C[C@H](C)[N@:3](->*)([H])CC[P:4](->*)(C)C)C=CN1C1=CC=CC=C1']
# CA = '[Mn+]'
# X = ComplexFromLigands(ligands, CA, 'OH')
# X.mol

# # 2
# ligands = ['[H-:1]->*','*<-[C-:5]#[O+]','*<-[C-:6]#[O+]','*<-[C-:2]1=[N+](C[C@H](C)[N@H:3](->*)CC[P:4](->*)(C)C)C=CN1C1=CC=CC=C1']
# Chem.MolFromSmiles(ligands[-1])
# CA = '[Mn+]'
# X = ComplexFromLigands(ligands, CA, 'OH')
# X.mol

# # 3 E/Z
# ligands = ['C/C=[CH-:1]\\[*] |C:2.2|']
# CA = '[Mn+]'
# geom = 'OH'
# X = ComplexFromLigands(ligands, CA, geom)
# X.AddConformer()

# # 4 - Bug: fixed
# ligands = ['C1=CC=CC=C1C(=[O:1])[O-:2]','C1=CC=CC=C1C(=[O:3])[O-:4]','C1=CC=CC=C1C(=[O:5])[O-:6]']
# CA = '[Mn+2]'
# geom = 'OH'
# X = ComplexFromLigands(ligands, CA, geom)
# X.AddConformer()
# X.ToXYZ('C:/Users/ivan-/Desktop/x.xyz')
# Xs = X.GetStereomers(drop_enantiomers = False)

