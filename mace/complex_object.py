'''
Chemoinformatic support for mononuclear complexes. Generates possible
stereoisomers and 3D coordinates. Supports square planar and octahedral
geometries
'''

#%% Imports

import json
from copy import deepcopy
from itertools import product, combinations

from rdkit import Chem
from rdkit.Chem import AllChem
#from rdkit.Geometry.rdGeometry import Point3D
import rdkit.Chem.rdDistGeom as rdDG

from .smiles_parsing import MolFromSmiles
from .parameters import params
from .supporting_functions import _CalcTHVolume, _RemoveRs


#%% Complex object

class Complex():
    '''
    Wrapper around RDKit Mol object for mononuclear complexes
    '''
    
    # object internal parameters
    _FFParams = params.FFParams
    _Rcov = params.Rcov
    _Syms = params.Syms
    _Geoms = params.Geoms
    _Bounds = params.Bounds
    _PosVs = params.PosVs
    _MinVs = params.MinVs
    _EqOrs = params.EqOrs
    _Nears = params.Nears
    _Angles = params.Angles
    
    
#%% Initialization
    
    def _CheckMol(self):
        # find and check dative bonds
        info = [(b.GetIdx(), b.GetBeginAtom(), b.GetEndAtom()) for b in self.mol.GetBonds() \
                if str(b.GetBondType()) == 'DATIVE']
        info.sort(key = lambda x: x[0])
        if len(info) == 0:
            raise ValueError('Bad SMILES: no dative bonds')
        # check number of central atoms
        CAs = set([_[2].GetIdx() for _ in info])
        if len(CAs) > 1:
            raise ValueError('Bad SMILES: there are several acceptors of dative bonds (central atoms)')
        # check CA's bonds
        self._idx_CA = list(CAs)[0]
        CA = self.mol.GetAtomWithIdx(self._idx_CA)
        if len(CA.GetBonds()) > len(info):
            raise ValueError('Bad SMILES: some bonds with central atom are not dative')
        if CA.GetNumImplicitHs() > 0:
            raise ValueError('Bad SMILES: all hydrogens (hydrides) bonded to central atom must be encoded explicitly with isotopic label')
        # check donor atoms labelling
        self._DAs = {_[1].GetIdx(): _[1].GetAtomMapNum() for _ in info}
        labs = list(self._DAs.values())
        if len(labs) > len([_ for _ in self._Geoms[self._geom] if str(_).isdigit()]):
            raise ValueError('Bad SMILES: number of donor atoms exceeds maximal possible for given geometry')
        # check donor atoms labelling
        if 0 in labs:
            self.err_init = 'Bad SMILES: some donor atoms don\'t have an isotopic label'
            return
        elif len(set(labs)) != len(labs):
            self.err_init = 'Bad SMILES: isotopic labels are not unique'
            return
        elif max(labs) > max([_ for _ in self._Geoms[self._geom] if str(_).isdigit()]):
            self.err_init = 'Bad SMILES: maximal isotopic label exceeds number of ligands'
            return
    
    
    def _SetComparison(self):
        '''
        Prepares set of all possible smiles to compare with other complexes
        '''
        mol_norm = deepcopy(self.mol)
        # fix resonance issues without ResonanceMolSupplier
        Chem.Kekulize(mol_norm, clearAromaticFlags = True)
        # modify X<-[C-]=[N+] fragments to X<-[C]-[N] (carbenes)
        for i, j, k in mol_norm.GetSubstructMatches(Chem.MolFromSmarts('[*]<-[C-]=[N+]')):
            mol_norm.GetAtomWithIdx(j).SetFormalCharge(0)
            mol_norm.GetAtomWithIdx(k).SetFormalCharge(0)
            mol_norm.GetBondBetweenAtoms(j, k).SetBondType(Chem.rdchem.BondType.SINGLE)
        # C([O-]->[*])=O->[*] to [C+]([O-]->[*])-[O-]->[*]
        for i, j, k in mol_norm.GetSubstructMatches(Chem.MolFromSmarts('[O-]C=[O]')):
            # TODO: check bonding better
            if mol_norm.GetAtomWithIdx(i).GetAtomMapNum() and \
               mol_norm.GetAtomWithIdx(k).GetAtomMapNum():
                mol_norm.GetAtomWithIdx(j).SetFormalCharge(1)
                mol_norm.GetAtomWithIdx(k).SetFormalCharge(-1)
                mol_norm.GetBondBetweenAtoms(j, k).SetBondType(Chem.rdchem.BondType.SINGLE)
        mol_norm = Chem.MolFromSmiles(Chem.MolToSmiles(mol_norm, canonical = False))
        # generate all ligand orientations
        mol = deepcopy(mol_norm)
        Chem.SetBondStereoFromDirections(mol)
        # remove dative bonds # HINT: remove after RDKit fix
        drop = []
        CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
        for b in mol.GetBonds():
            if str(b.GetBondType()) == 'DATIVE':
                if b.GetBeginAtom().GetChiralTag() in CHIs:
                    continue
                drop.append( (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) )
        # create bond with dummy
        ed = Chem.EditableMol(mol)
        for i, j in drop:
            ed.RemoveBond(i, j)
        mol = ed.GetMol()
        Chem.SanitizeMol(mol)
        # invert mol
        mol_inv = deepcopy(mol)
        CHIs = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]
        for atom in mol_inv.GetAtoms():
            # revert all chiral centers
            tag = atom.GetChiralTag()
            if tag in CHIs:
                atom.SetChiralTag(CHIs[not CHIs.index(tag)])
        # generate structure descriptor
        EqOrs = self._EqOrs[self._geom]
        if 'enant' + self._geom in self._EqOrs:
            EqOrsInv = self._EqOrs['enant' + self._geom]
        else:
            EqOrsInv = self._EqOrs[self._geom]
        _ID = []
        _eID = []
        _DAs = {_.GetIdx(): _.GetAtomMapNum() for _ in mol.GetAtoms() if _.GetAtomMapNum()}
        _DAs_inv = {_.GetIdx(): _.GetAtomMapNum() for _ in mol_inv.GetAtoms() if _.GetAtomMapNum()}
        for i in range(len(EqOrs[1])):
            for idx, num in _DAs.items():
                mol.GetAtomWithIdx(idx).SetAtomMapNum(EqOrs[num][i])
                mol.GetAtomWithIdx(idx).SetIsotope(EqOrs[num][i])
            for idx, num in _DAs_inv.items():
                mol_inv.GetAtomWithIdx(idx).SetAtomMapNum(EqOrsInv[num][i])
                mol_inv.GetAtomWithIdx(idx).SetIsotope(EqOrsInv[num][i])
            # add basic mol
            _ID.append(Chem.CanonSmiles(Chem.MolToSmiles(mol)))
            _eID.append(Chem.CanonSmiles(Chem.MolToSmiles(mol_inv)))
            # add resonance structures to mols
            if self.maxResonanceStructures > 1:
                idx = 0
                for m in Chem.ResonanceMolSupplier(mol):
                    if idx >= self.maxResonanceStructures:
                        break
                    _ID.append( Chem.CanonSmiles(Chem.MolToSmiles(m)) )
                    idx += 1
                # same for inv mol
                idx = 0
                for m in Chem.ResonanceMolSupplier(mol_inv):
                    if idx >= self.maxResonanceStructures:
                        break
                    _eID.append( Chem.CanonSmiles(Chem.MolToSmiles(m)) )
                    idx += 1
        self._ID = set(_ID)
        self._eID = set(_eID)
    
    
    def __init__(self, smiles, geom, maxResonanceStructures = 1):
        '''
        Generates Complex object from SMILES with stereo info encoded as
        atomic map numbers of donor atoms
        '''
        self._smiles_init = smiles
        self._geom = geom
        try:
            self.maxResonanceStructures = int(maxResonanceStructures)
        except TypeError:
            raise TypeError('Bad maximal number of resonance structures: must be an integer')
        if self.maxResonanceStructures < 0:
            raise ValueError('Bad maximal number of resonance structures: must be zero or positive')
        self.err_init = None
        # check geom
        if self._geom not in self._Geoms:
            raise ValueError(f'Unknown geometry type: {self._geom}')
        # check smiles
        self.mol = MolFromSmiles(self._smiles_init)
        if not self.mol:
            raise ValueError('Bad SMILES: not readable')
        # check mol
        self._CheckMol()
        # other stuff
        self._SetComparison()
        self.mol3D = Chem.AddHs(self.mol)
        self._embedding_prepared = False
        self._ff_prepared = False
    
    
    def _RaiseErrorInit(self):
        '''
        Raises warning if complex does not have enough stereo info
        '''
        if self.err_init:
            message = [self.err_init,
                       '',
                       'The initial SMILES contains insufficient or erroneous info',
                       'on the positions of the ligands around the central atom',
                       'encoded with isotopic labels.', '',
                       'To use 3D generation and other features, first generate',
                       'possible stereomers', '']
            raise ValueError('\n'.join(message))
        
        return False
    
    
    def IsEqual(self, X):
        '''
        Compares two Complex objects for equivalence
        '''
        self._RaiseErrorInit()
        
        return bool(self._ID.intersection(X._ID))
    
    
    def IsEnantiomeric(self):
        '''
        Checks is the complex enantiomeric
        '''
        self._RaiseErrorInit()
        
        return not bool(self._ID.intersection(self._eID))
    
    
    def IsEnantiomer(self, X):
        '''
        Checks if two complexes are enantiomers
        '''
        self._RaiseErrorInit()
        
        return bool(self._ID.intersection(X._eID))
    
    
#%% Stereomers search
    
    def _FindNeighboringDAs(self, minTransCycle = None):
        '''
        Finds restrictions of multidentate ligands.
        Returns pairs of DAs which must be near each other
        '''
        # get DAs
        DAs = list(self._DAs.keys())
        # get rings
        rings = [list(r) for r in Chem.GetSymmSSSR(self.mol)]
        rings = [r for r in rings if self._idx_CA in r]
        # find restrictions
        restrictions = []
        for r in rings:
            if minTransCycle and len(r) >= minTransCycle:
                continue
            r = [idx for idx in r if idx in DAs]
            if len(r) == 2:
                restrictions.append(r)
        
        return restrictions
    
    
    def _FindMerOnly(self):
        '''
        Finds restrictions for rigid X-Y-Z fragments (fac- geometry is impossible)
        Returns pairs of DAs which must not be near each other
        '''
        # get DAs
        DAs = list(self._DAs.keys())
        # get rings
        rings = [list(r) for r in Chem.GetSymmSSSR(self.mol)]
        rings = [r for r in rings if self._idx_CA in r]
        # get neighboring DAs and the corresponding paths
        neighbors = {}
        paths = {}
        for r in rings:
            i, j = [idx for idx in r if idx in DAs]
            # set neighbors
            if i not in neighbors:
                neighbors[i] = [j]
            else:
                neighbors[i] += [j]
            if j not in neighbors:
                neighbors[j] = [i]
            else:
                neighbors[j] += [i]
            # path
            idx = r.index(self._idx_CA)
            path = r[idx+1:] + r[:idx]
            # TODO: check path
            if path[0] == i:
                paths[(i,j)] = path
                paths[(j,i)] = path[::-1]
            elif path[0] == j:
                paths[(i,j)] = path[::-1]
                paths[(j,i)] = path
        # drop DAs with less than 2 neighbors
        drop = [idx for idx, ns in neighbors.items() if len(ns) < 2]
        for idx in drop:
            neighbors.pop(idx)
        # TODO: improve rotability check
        restricted = lambda i, j: i + j < 4 or i < 2 and j == 3 or j < 2 and i == 3
        # check rigidity of central DA
        restrictions = []
        for idx, ns in neighbors.items():
            a = self.mol.GetAtomWithIdx(idx)
            flag = False
            # sp2 or carbenes
            if a.GetSymbol() in ('C', 'N') and str(a.GetHybridization()) == 'SP2' or \
               a.GetSymbol() == 'C' and a.GetNumRadicalElectrons() == 2:
                flag = True
            # check conjugated pyrrol-like anion
            if a.GetSymbol() == 'N' and a.GetFormalCharge() == -1:
                hybr = [(_.GetSymbol(), str(_.GetHybridization())) for _ in a.GetNeighbors() if _.GetIdx() != self._idx_CA]
                if ('N', 'SP2') in hybr or ('C', 'SP2') in hybr:
                    flag = True
            # drop flexible
            if not flag:
                continue
            # check number of rotable bonds between "idx" and "ns"
            rot_bonds = {}
            for n in ns:
                path = paths[(idx, n)][1:-1]
                if len(path) < 2:
                    rot_bonds[n] = 0
                counter = 0
                for i in range(len(path)-1):
                    b = self.mol.GetBondBetweenAtoms(path[i], path[i+1])
                    if str(b.GetBondType()) == 'SINGLE':
                        counter += 1
                rot_bonds[n] = counter
            # add final restrictions
            restrictions += [(i, j) for i, j in combinations(ns, r = 2) if restricted(rot_bonds[i], rot_bonds[j])]
        
        return restrictions
    
    
    def GetStereomers(self, regime = 'all', dropEnantiomers = True,
                      minTransCycle = None, merRule = True):
        '''
        Generates all possible stereomers of a complex.
        Saves stereochemistry of existing centers.
        Three regimes are available:
            - "CA": only stereochemistry of central atom;
            - "ligands": only stereochemistry of ligands;
            - "all": both "CA" and "ligands" regimes.
        '''
        if regime not in ('CA', 'ligands', 'all'):
            raise ValueError('Regime variable bad value: must be one of "CA", "ligands", "all"')
        if type(merRule) is not bool:
            raise ValueError('Bad meridial-rule: must be True or False')
        # set Mol object
        mol = deepcopy(self.mol)
        if regime == 'ligands':
            # check numbering
            if self.err_init:
                raise ValueError('Stereo info for the central atom is not specified correctly. Use "CA" or "all" regimes to fix that')
        else:
            # randomly set isotopic numbers
            DAs = list(self._DAs.keys())
            for num, idx in enumerate(DAs):
                mol.GetAtomWithIdx(idx).SetAtomMapNum(num + 1)
                mol.GetAtomWithIdx(idx).SetIsotope(num + 1)
        # generate needed stereomers
        if regime == 'CA':
            mols = [mol]
        else:
            idxs = [idx for idx, chi in Chem.FindMolChiralCenters(mol, includeUnassigned = True) if chi == '?']
            idxs = [idx for idx in idxs if idx != self._idx_CA]
            # drop P/As with nH > 2 # HINT: remove after fixing RDKit #3773
            drop = []
            for idx in idxs:
                a = mol.GetAtomWithIdx(idx)
                if a.GetSymbol() not in ('P', 'As'):
                    continue
                nHs = a.GetNumExplicitHs() + len([_ for _ in a.GetNeighbors() if _.GetSymbol() == 'H'])
                nXs = len([_ for _ in a.GetNeighbors() if _.GetSymbol() == '*'])
                if nHs + nXs >= 2:
                    drop.append(idx)
            idxs = [_ for _ in idxs if _ not in drop]
            # generate all possible combinations of stereocentres
            mols = []
            for chis in product([Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW], repeat = len(idxs)):
                m = deepcopy(mol)
                for idx, chi in zip(idxs, chis):
                    m.GetAtomWithIdx(idx).SetChiralTag(chi)
                mols.append(m)
            # drop similar ones
            smiles = [Chem.MolToSmiles(m) for m in mols]
            drop = []
            for i in range(len(smiles)-1):
                if i in drop:
                    continue
                for j in range(i+1, len(smiles)):
                    if smiles[i] == smiles[j]:
                        drop.append(j)
            mols = [m for i, m in enumerate(mols) if i not in drop]
        # generate all CA isomers
        if regime == 'ligands':
            # transform mols to Complex objects and return them
            stereomers = []
            for m in mols:
                stereomers.append( Complex(Chem.MolToSmiles(m), self._geom, self.maxResonanceStructures) )
            return stereomers
        # find restrictions on DA positions
        pairs = self._FindNeighboringDAs(minTransCycle)
        mers = self._FindMerOnly() if merRule else []
        # generate all possible CA orientations
        stereomers = []
        for m in mols:
            addend = []
            for idx_sym in sorted(list(self._Syms[self._geom].keys())):
                sym = self._Syms[self._geom][idx_sym]
                m1 = deepcopy(m)
                # set new isotopes
                info = {}
                for idx, a_idx in enumerate(DAs):
                    num = sym.index(idx + 1) + 1
                    m1.GetAtomWithIdx(a_idx).SetAtomMapNum(num)
                    m1.GetAtomWithIdx(a_idx).SetIsotope(num)
                    info[a_idx] = num
                # check neighboring DAs restriction
                drop = False
                for idx_a, idx_b in pairs:
                    if info[idx_b] not in self._Nears[self._geom][info[idx_a]]:
                        drop = True
                if drop:
                    continue
                # check mer DAs restriction
                drop = False
                for idx_a, idx_b in mers:
                    if info[idx_b] in self._Nears[self._geom][info[idx_a]]:
                        drop = True
                if drop:
                    continue
                addend.append(m1)
            addend = [Complex(Chem.MolToSmiles(m), self._geom, self.maxResonanceStructures) for m in addend]
            # filter uniques
            drop = []
            for i in range(len(addend)-1):
                if i in drop:
                    continue
                for j in range(i+1, len(addend)):
                    if addend[i].IsEqual(addend[j]):
                        drop.append(j)
            addend = [compl for i, compl in enumerate(addend) if i not in drop]
            # add to main
            stereomers += addend
        # final filtering
        drop = []
        for i in range(len(stereomers)-1):
            if i in drop:
                continue
            for j in range(i+1, len(stereomers)):
                if stereomers[i].IsEqual(stereomers[j]):
                    drop.append(j)
                elif dropEnantiomers and stereomers[i].IsEnantiomer(stereomers[j]):
                    drop.append(j)
        stereomers = [compl for i, compl in enumerate(stereomers) if i not in drop]
        
        return stereomers    
    
    
#%% 3D Generation
    
    def _SetEmbedding(self):
        '''
        Sets cordinates' constraints
        '''
        # find dummies-helpers
        add = [idx for idx in self._Geoms[self._geom] if 'X' in str(idx)]
        must = set([idx for idx in self._Geoms[self._geom] if str(idx).isdigit()])
        have = set([num for num in self._DAs.values()])
        add += list(must.difference(have))
        # prepare coordMap
        self._coordMap = {self._idx_CA: self._Geoms[self._geom]['CA']}
        for idx, num in self._DAs.items():
            self._coordMap[idx] = self._Geoms[self._geom][num]
        # add dummies-helpers to mol3Dx and coordMap
        self._dummies = {}
        ed = Chem.EditableMol(self.mol3D)
        for num in add:
            idx = ed.AddAtom(Chem.Atom(0))
            ed.AddBond(idx, self._idx_CA, Chem.BondType.DATIVE)
            self._dummies[idx] = num
            self._coordMap[idx] = self._Geoms[self._geom][num]
        self.mol3Dx = ed.GetMol()
        Chem.SanitizeMol(self.mol3Dx)
        # prepare bounds matrix
        X = rdDG.GetMoleculeBoundsMatrix(self.mol3Dx)
        CS = [(self._idx_CA, 'CA')] + list(self._DAs.items()) + list(self._dummies.items())
        for (i, num1), (j, num2) in combinations(CS, r = 2):
            dmax = self._Bounds[self._geom][num1][num2]
            dmin = self._Bounds[self._geom][num2][num1]
            X[min(i,j)][max(i,j)] = max(dmin, dmax)
            X[max(i,j)][min(i,j)] = min(dmin, dmax)
        self._boundsMatrix = X
        # final
        self._embedding_prepared = True
    
    
    def _SetCentralAtomAngles(self):
        '''
        Sets restrictions to L->X<-L angles
        '''
        DAs = list(self._DAs.items()) + list(self._dummies.items())
        for i, j in combinations(range(len(DAs)), r = 2):
            a1_idx, a1_num = DAs[i]
            a2_idx, a2_num = DAs[j]
            angle = self._Angles[self._geom][a1_num][a2_num]
            if a1_num in self._Nears[self._geom][a2_num]:
                k = self._FFParams['kZ-LXL']
            else:
                k = self._FFParams['kE-LXL']
            constraint = [a1_idx, self._idx_CA, a2_idx, False, angle, angle, k]
            self._angle_params.append(constraint)
            self._ff.UFFAddAngleConstraint(*constraint)
    
    
    def _SetCentralAtomBonds(self):
        '''
        Sets restrictions to X<-L bonds
        '''
        for idx in self._DAs:
            dist = self._Rcov[self.mol3Dx.GetAtomWithIdx(self._idx_CA).GetAtomicNum()] + \
                   self._Rcov[self.mol3Dx.GetAtomWithIdx(idx).GetAtomicNum()]
            constraint = [idx, self._idx_CA, False, dist, dist, self._FFParams['kXL']]
            self._bond_params.append(constraint)
            self._ff.UFFAddDistanceConstraint(*constraint)
        # dummies-helpers
        for idx, num in self._dummies.items():
            if 'X' in str(num):
                constraint = [idx, self._idx_CA, False, self._FFParams['X*'], self._FFParams['X*'], self._FFParams['kX*']]
                self._bond_params.append(constraint)
                self._ff.UFFAddDistanceConstraint(*constraint)
            else:
                dist = self._Rcov[self.mol3Dx.GetAtomWithIdx(self._idx_CA).GetAtomicNum()] + \
                       self._Rcov[self.mol3Dx.GetAtomWithIdx(idx).GetAtomicNum()]
                constraint = [idx, self._idx_CA, False, dist, dist, self._FFParams['kXL']]
                self._bond_params.append(constraint)
                self._ff.UFFAddDistanceConstraint(*constraint)
    
    
    def _SetDonorAtomsAngles(self):
        '''
        Sets restrictions to X<-L-A angles
        '''
        for DA in self._DAs:
            # get neighbors
            ns = self.mol3Dx.GetAtomWithIdx(DA).GetNeighbors()
            ns = [n.GetIdx() for n in ns if n.GetIdx() != self._idx_CA]
            if not ns or len(ns) > 3:
                continue
            # set angles
            if len(ns) == 1 and str(self.mol3Dx.GetAtomWithIdx(ns[0]).GetHybridization()) == 'SP2':
                angle = self._FFParams['XLO']
                k = self._FFParams['kXLO']
            else:
                angle = {1: self._FFParams['XLA'],
                         2: self._FFParams['XLA2'],
                         3: self._FFParams['XLA3']}[len(ns)]
                k = self._FFParams['kXLA']
            for n in ns:
                constraint = [self._idx_CA, DA, n, False, angle, angle, k]
                self._angle_params.append(constraint)
                self._ff.UFFAddAngleConstraint(*constraint)
    
    
    def _SetDonorAtomsParams(self):
        '''
        Sets DA parameters without CA
        '''
        for DA in self._DAs:
            # get neighbors
            ns = [_.GetIdx() for _ in self.mol3Dx.GetAtomWithIdx(DA).GetNeighbors()]
            ns = [_ for _ in ns if _ != self._idx_CA]
            if not ns:
                continue
            N = len(ns)
            # bonds
            for n in ns:
                d = self._Rcov[self.mol3Dx.GetAtomWithIdx(n).GetAtomicNum()] + \
                    self._Rcov[self.mol3Dx.GetAtomWithIdx(DA).GetAtomicNum()]
                constraint = [DA, n, False, d, d, self._FFParams['kLA']]
                self._bond_params.append(constraint)
                self._ff.UFFAddDistanceConstraint(*constraint)
            # X<-L-A angles
            if N < 2 or N > 3:
                continue
            for n1, n2 in combinations(ns, r = 2):
                a = self._FFParams['XLA2'] if N == 2 else self._FFParams['XLA3']
                constraint = [n1, DA, n2, False, a, a, self._FFParams['kALA']]
                self._angle_params.append(constraint)
                self._ff.UFFAddAngleConstraint(*constraint)
            # L-A-B angles
            for n in ns:
                atom = self.mol3Dx.GetAtomWithIdx(n)
                n2s = [_.GetIdx() for _ in atom.GetNeighbors()]
                N = len(n2s)
                if N < 3 or N > 4:
                    continue
                for n21, n22 in combinations(n2s, r = 2):
                    if str(atom.GetHybridization()) == 'SP2':
                        a = self._FFParams['XLA2']
                        constraint = [n21, n, n22, False, a, a, self._FFParams['kALA']]
                        self._angle_params.append(constraint)
                        self._ff.UFFAddAngleConstraint(*constraint)
                    elif str(atom.GetHybridization()) == 'SP3':
                        a = self._FFParams['XLA3']
                        constraint = [n21, n, n22, False, a, a, self._FFParams['kALA']]
                        self._angle_params.append(constraint)
                        self._ff.UFFAddAngleConstraint(*constraint)
    
    
    def _SetDummiesBonds(self):
        '''
        Sets DA parameters without CA
        '''
        # get list of dummies
        for a in self.mol3Dx.GetAtoms():
            if a.GetSymbol() != '*' or a.GetAtomMapNum():
                continue
            # dummies bonded to DA was already treated
            flags = [n.GetIdx() in self._DAs.keys() for n in a.GetNeighbors()]
            if True in flags:
                continue
            for n in a.GetNeighbors():
                if n.GetIdx() in self._DAs.keys():
                    continue
                # set k
                d = self._Rcov[n.GetAtomicNum()] + self._Rcov[a.GetAtomicNum()]
                constraint = [a.GetIdx(), n.GetIdx(), False, d, d, self._FFParams['kA*']]
                self._bond_params.append(constraint)
                self._ff.UFFAddDistanceConstraint(*constraint)
    
    
    def _SetForceField(self, confId):
        '''
        Sets force field for the geometry optimization
        '''
        self._ff = AllChem.UFFGetMoleculeForceField(self.mol3Dx, confId = confId)
        if self._ff_prepared:
            # restore them from saved params
            for constraint in self._angle_params:
                self._ff.UFFAddAngleConstraint(*constraint)
            for constraint in self._bond_params:
                self._ff.UFFAddDistanceConstraint(*constraint)
            return
        # set ff parameters
        self._angle_params = []
        self._bond_params = []
        self._SetCentralAtomAngles()
        self._SetCentralAtomBonds()
        self._SetDonorAtomsAngles()
        self._SetDonorAtomsParams()
        self._SetDummiesBonds()
        self._ff_prepared = True
    
    
    def _CheckStereoCA(self, confId):
        '''
        Checks that generated coordinates corresponds to initial chirality
        '''
        conf = self.mol3Dx.GetConformer(confId)
        DAs = {val: key for key, val in self._DAs.items()}
        DAs['CA'] = self._idx_CA
        for key, val in self._dummies.items():
            DAs[val] = key
        Vs = []
        for idxs in self._PosVs[self._geom]:
            Vs.append(_CalcTHVolume(conf, [DAs[idx] for idx in idxs]))
        
        return sum(Vs) > self._MinVs[self._geom] # sum(Vs) > 0
    
    
    def Optimize(self, confId = 0, maxIts = 1000):
        '''
        Optimizes Complex
        '''
        self._RaiseErrorInit()
        # optimization
        self._SetForceField(confId)
        self._ff.Initialize()
        flag = self._ff.Minimize(maxIts = maxIts)
        # energy
        E = self._ff.CalcEnergy()
        conf3Dx = self.mol3Dx.GetConformer(confId)
        conf3Dx.SetDoubleProp('E', E)
        # synchronize with mol3D
        conf3D = self.mol3D.GetConformer(confId)
        conf3D.SetDoubleProp('E', E)
        for atom in self.mol3D.GetAtoms():
            idx = atom.GetIdx()
            conf3D.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
        # synchronize with mol
        conf = self.mol.GetConformer(confId)
        conf.SetDoubleProp('E', E)
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            conf.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
        
        return flag
    
    
    def AddConformer(self, clearConfs = True, useRandomCoords = True, maxAttempts = 10):
        '''
        Generates complex conformer using constrained embedding
        '''
        self._RaiseErrorInit()
        # check embedding prerequisites
        if not self._embedding_prepared:
            self._SetEmbedding()
        # set embedding parameters
        params = rdDG.EmbedParameters()
        params.clearConfs = clearConfs
        params.enforceChirality = True
        params.useRandomCoords = useRandomCoords
        #params.embedFragmentsSeparately = False
        params.SetBoundsMat(self._boundsMatrix)
        # embedding
        flag = -1
        attempt = maxAttempts
        while flag == -1 and attempt > 0:
            attempt -= 1
            # embedding
            flag = AllChem.EmbedMolecule(self.mol3Dx, params)
            if flag == -1:
                continue
            # optimization # HINT: do not use self.Optimize as we need to apply self._CheckStereoCA after
            self._SetForceField(flag)
            self._ff.Initialize()
            self._ff.Minimize(maxIts = 1000)
            # check chiral centers from 3D
            if not self._CheckStereoCA(flag):
                self.mol3Dx.RemoveConformer(flag)
                flag = -1
                continue
            # refresh other confs
            if clearConfs:
                self.mol.RemoveAllConformers()
                self.mol3D.RemoveAllConformers()
            # move CA to (0,0,0)
            conf3Dx = self.mol3Dx.GetConformer(flag)
            r = deepcopy(conf3Dx.GetAtomPosition(self._idx_CA))
            for i in range(conf3Dx.GetNumAtoms()):
                conf3Dx.SetAtomPosition(i, conf3Dx.GetAtomPosition(i) - r)
            # energy
            E = self._ff.CalcEnergy()
            conf3Dx.SetDoubleProp('E', E)
            conf3Dx.SetDoubleProp('EmbedRMS', -1)
            # synchronize with mol3D
            conf3D = AllChem.Conformer()
            conf3D.SetDoubleProp('E', E)
            conf3D.SetDoubleProp('EmbedRMS', -1)
            for atom in self.mol3D.GetAtoms():
                idx = atom.GetIdx()
                conf3D.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            conf3D.SetId(flag)
            self.mol3D.AddConformer(conf3D, assignId = True)
            # synchronize with mol
            conf = AllChem.Conformer()
            conf.SetDoubleProp('E', E)
            conf.SetDoubleProp('EmbedRMS', -1)
            for atom in self.mol.GetAtoms():
                idx = atom.GetIdx()
                conf.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            conf.SetId(flag)
            self.mol.AddConformer(conf, assignId = True)
        
        return flag
    
    
    def AddConstrainedConformer(self, core, confId = 0, clearConfs = True,
                                useRandomCoords = True, maxAttempts = 10,
                                engine = 'coordMap', deltaR = 0.01):
        '''
        Constrained embedding using other complex geometry
        Core complex must be a substructure of complex and
        must contain CA and all DAs (excluding supportive dummies)
        '''
        if len(Chem.GetMolFrags(core.mol)) != 1:
            raise ValueError('Bad core: core must contain exactly one fragment')
        if engine not in ('coordMap', 'boundsMatrix'):
            raise ValueError('Unknown engine: must be one of "coordMap" or "boundsMatrix"')
        # make mol3Dx and mol3D
        if not self._embedding_prepared:
            self._SetEmbedding()
        # prepare molecule for substructure search
        core_mol = _RemoveRs(deepcopy(core.mol)) # HINT: cannot convert to SMARTS: RDKIT #3774
        # substructure check
        match = self.mol3Dx.GetSubstructMatch(core_mol, useChirality = True)
        if not match:
            raise ValueError('Bad core: core is not a substructure of the complex')
        if self._idx_CA not in match:
            raise ValueError('Bad core: core must contain the central atom')
        # prepare embed params
        coordMap = {}
        coreConf = core_mol.GetConformer(confId)
        for i, idxI in enumerate(match):
            coordMap[idxI] = coreConf.GetAtomPosition(i)
        # check if some DAs missed in new coordMap
        add = [idx for idx in self._coordMap if idx not in match]
        if add:
            # make dummy mol
            dummy = Chem.MolFromSmiles('.'.join(['[*]']*len(self._coordMap))) # CA
            dummyMap = {}
            conf = Chem.Conformer()
            for i, (idx, point) in enumerate(self._coordMap.items()):
                dummyMap[i] = idx
                conf.SetAtomPosition(i, point)
            dummy.AddConformer(conf)
            # orient dummy mol over core
            algMap = [(key, match.index(val)) for key, val in dummyMap.items() if val in match]
            AllChem.AlignMol(dummy, core_mol, atomMap = algMap, maxIters = 200)
            # renew coordmap
            dummyMap = {val: key for key, val in dummyMap.items()}
            for idx in add:
                coordMap[idx] = dummy.GetConformer().GetAtomPosition(dummyMap[idx])
        # set bounds matrix
        BM = rdDG.GetMoleculeBoundsMatrix(self.mol3Dx)
        for (i, ri), (j, rj) in combinations(coordMap.items(), r = 2):
            d = sum([_**2 for _ in list(ri-rj)])**0.5
            BM[min(i,j)][max(i,j)] = d + deltaR
            BM[max(i,j)][min(i,j)] = d - deltaR
        # embedding parameters
        params = rdDG.EmbedParameters()
        params.clearConfs = clearConfs
        params.enforceChirality = True
        params.useRandomCoords = useRandomCoords
        #params.embedFragmentsSeparately = False
        params.SetBoundsMat(BM)
        # embedding
        flag = -1
        attempt = maxAttempts
        while flag == -1 and attempt > 0:
            attempt -= 1
            if engine == 'coordMap':
                flag = AllChem.EmbedMolecule(self.mol3Dx, coordMap = coordMap,
                                             clearConfs = clearConfs,
                                             useRandomCoords = useRandomCoords,
                                             enforceChirality = True)
            elif engine == 'boundsMatrix':
                flag = AllChem.EmbedMolecule(self.mol3Dx, params)
            if flag == -1:
                continue
            # set ff
            self._SetForceField(flag)
            # reorient core
            algMap = [(j, i) for i, j in enumerate(match)]
            AllChem.AlignMol(self.mol3Dx, core_mol, atomMap = algMap, maxIters = 200)
            # add tethers
            conf = core_mol.GetConformer(confId)
            for i in range(core_mol.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                pIdx = self._ff.AddExtraPoint(p.x, p.y, p.z, fixed = True) - 1
                self._ff.UFFAddDistanceConstraint(pIdx, match[i], False, 0, 0, 100.)
            # optimize
            self._ff.Initialize()
            self._ff.Minimize(maxIts = 1000)
            rms = AllChem.AlignMol(self.mol3Dx, core_mol, atomMap = algMap, maxIters = 200)
            # check chiral centers from 3D
            if not self._CheckStereoCA(flag):
                self.mol3Dx.RemoveConformer(flag)
                flag = -1
                continue
            # refresh other confs
            if clearConfs:
                self.mol.RemoveAllConformers()
                self.mol3D.RemoveAllConformers()
            # energy
            E = self._ff.CalcEnergy()
            conf3Dx = self.mol3Dx.GetConformer(flag)
            conf3Dx.SetDoubleProp('E', E)
            conf3Dx.SetDoubleProp('EmbedRMS', rms)
            # synchronize with mol3D
            conf3D = AllChem.Conformer()
            conf3D.SetDoubleProp('E', E)
            conf3D.SetDoubleProp('EmbedRMS', rms)
            for atom in self.mol3D.GetAtoms():
                idx = atom.GetIdx()
                conf3D.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            self.mol3D.AddConformer(conf3D, assignId = True)
            # synchronize with mol
            conf = AllChem.Conformer()
            conf.SetDoubleProp('E', E)
            conf.SetDoubleProp('EmbedRMS', rms)
            for atom in self.mol.GetAtoms():
                idx = atom.GetIdx()
                conf.SetAtomPosition(idx, conf3Dx.GetAtomPosition(idx))
            self.mol.AddConformer(conf, assignId = True)        
        
        return flag
    
    
#%% 3D support methods
    
    def GetNumConformers(self):
        '''
        Returns number of conformers
        '''
        
        return self.mol.GetNumConformers()
    
    
    def RemoveConformer(self, confId):
        '''
        Removes conformer with given idx
        '''
        self.mol.RemoveConformer(confId)
        self.mol3D.RemoveConformer(confId)
        self.mol3Dx.RemoveConformer(confId)
    
    
    def RemoveAllConformers(self):
        '''
        Removes all conformers
        '''
        self.mol.RemoveAllConformers()
        self.mol3D.RemoveAllConformers()
        self.mol3Dx.RemoveAllConformers()
    
    
    def GetMinEnergyConfId(self, i):
        '''
        Returns idx of Conformer with the i-th lowest energy.
        Returns None if i >= NumConfs
        '''
        N = self.GetNumConformers()
        if i >= N:
            raise ValueError('Bad conformer ID')
        # order Es
        Es = [(self.mol.GetConformer(idx).GetDoubleProp('E'), idx) for idx in range(N)]
        Es.sort()
        
        return Es[i][1]
    
    
    def AddConformers(self, numConfs = 10, clearConfs = True, useRandomCoords = True,
                      maxAttempts = 10, rmsThresh = -1):
        '''
        Generates several conformers
        '''
        self._RaiseErrorInit()
        # generate 3D
        flags = []
        for i in range(numConfs):
            clearConfsIter = False if flags else clearConfs
            flag = self.AddConformer(clearConfs = clearConfsIter,
                                     useRandomCoords = useRandomCoords,
                                     maxAttempts = maxAttempts)
            # check flag and rms
            if flag == -1:
                continue
            if rmsThresh == -1:
                flags.append(flag)
                continue
            # check rms with previous conformers
            remove_conf = False
            for cid in flags:
                rms = AllChem.GetConformerRMS(self.mol3D, cid, flag)
                if rms < rmsThresh:
                    remove_conf = True
                    break
            if remove_conf:
                self.RemoveConformer(flag)
            else:
                flags.append(flag)
        
        return flags
    
    
    def AddConstrainedConformers(self, core, confId = 0, numConfs = 10,
                                 clearConfs = True, useRandomCoords = True,
                                 maxAttempts = 10, engine = 'coordMap',
                                 deltaR = 0.01, rmsThresh = -1):
        '''
        Generates several conformers
        '''
        self._RaiseErrorInit()
        # generate 3D
        flags = []
        for i in range(numConfs):
            clearConfsIter = False if flags else clearConfs
            flag = self.AddConstrainedConformer(core, confId = confId,
                                                clearConfs = clearConfsIter,
                                                useRandomCoords = useRandomCoords,
                                                maxAttempts = maxAttempts,
                                                engine = engine, deltaR = deltaR)
            # check flag and rms
            if flag == -1:
                continue
            if rmsThresh == -1:
                flags.append(flag)
                continue
            # check rms with previous conformers
            remove_conf = False
            for cid in flags:
                rms = AllChem.GetConformerRMS(self.mol3D, cid, flag)
                if rms < rmsThresh:
                    remove_conf = True
                    break
            if remove_conf:
                self.RemoveConformer(flag)
            else:
                flags.append(flag)
        
        return flags
    
    
#%% MolSimplify helper
    
    def GetBondedLigand(self, num):
        '''
        Extracts ligand with available 3D from the embedded complex
        num is an isotopic number of DA owned by target ligand
        '''
        self._RaiseErrorInit()
        N = self.mol3D.GetNumConformers()
        if not N:
            raise ValueError('Complex has no conformers')
        # is num in DAs
        idx_DA = None
        for DA, isotope in self._DAs.items():
            if isotope == num:
                idx_DA = DA
        if idx_DA is None:
            raise ValueError(f'Complex has donor atom with {num} order number')
        # remove dative bonds
        ed = Chem.EditableMol(self.mol3D)
        for DA in self._DAs:
            ed.RemoveBond(DA, self._idx_CA)
        mol = ed.GetMol()
        # remove unneeded fragments
        frags = Chem.GetMolFrags(mol)
        drop = []
        for frag in frags:
            if idx_DA not in frag:
                drop += list(frag)
        drop = sorted(drop, reverse = True)
        ed = Chem.EditableMol(mol)
        for idx in drop:
            ed.RemoveAtom(idx)
        mol = ed.GetMol()
        Chem.SanitizeMol(mol)
        
        return mol
    
    
#%% Output
    
    def _ConfToXYZ(self, confId):
        '''
        Generates text of XYZ file of conformer
        '''
        # coordinates
        xyz = []
        conf = self.mol3Dx.GetConformer(confId) # not mol3D as it uses in AlignMol
        for atom in self.mol3D.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol == '*':
                symbol = 'X'
            pos = conf.GetAtomPosition(atom.GetIdx())
            line = f'{symbol:2} {pos.x:>-10.4f} {pos.y:>-10.4f} {pos.z:>-10.4f}'
            xyz.append(line)
        # conf params
        E = conf.GetDoubleProp('E')
        rms = conf.GetDoubleProp('EmbedRMS')
        # mol smiles
        mol = deepcopy(self.mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        smiles = Chem.MolToSmiles(mol, canonical = False)
        # mol3D smiles
        mol3D = deepcopy(self.mol3D)
        for atom in mol3D.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        smiles3D = Chem.MolToSmiles(mol3D, canonical = False)
        # mol3Dx smiles
        mol3Dx = deepcopy(self.mol3Dx)
        for atom in mol3Dx.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        smiles3Dx = Chem.MolToSmiles(mol3Dx, canonical = False)
        # dummies' coords
        dummies = []
        conf3Dx = mol3Dx.GetConformer(confId)
        for i in range(mol3D.GetNumAtoms(), mol3Dx.GetNumAtoms()):
            if i < mol3D.GetNumAtoms():
                continue
            p = conf3Dx.GetAtomPosition(i)
            dummies += [p.x, p.y, p.z]
        # make text
        info = {'conf': confId, 'E': float(f'{E:.2f}'),
                'rms': float(f'{rms:.4f}'), 'geom': self._geom,
                'total_charge': sum([a.GetFormalCharge() for a in self.mol.GetAtoms()]),
                'CA_charge': self.mol.GetAtomWithIdx(self._idx_CA).GetFormalCharge(),
                'smiles': smiles, 'smiles3D': smiles3D,
                'smiles3Dx': smiles3Dx, 'dummies': dummies}
        text = [str(len(xyz)), json.dumps(info)] + xyz
        
        return '\n'.join(text)+'\n'
    
    
    def ToXYZBlock(self, confId = 'all'):
        '''
        Returns XYZ as text block
        If confId is 'min' or -1, conformer with the lowest energy will be saved
        If confId is in 'min(i)', i > 0 format, conformer with the i-th lowest
        energy will be saved
        If confId is 'all' or -2, all conformers will be saved
        '''
        self._RaiseErrorInit()
        N = self.mol3D.GetNumConformers()
        if not N:
            raise ValueError('Bad conformer ID: complex has no conformers')
        # prepare conf idxs
        if confId in (-2, 'all'):
            conf_idxs = list(range(N))
        elif confId in (-1, 'min'):
            conf_idxs = [self.GetMinEnergyConfId(0)]
        elif str(confId)[:4] == 'min(' and str(confId)[-1] == ')':
            idx = str(confId)[4:-1]
            if not idx.isdigit():
                raise ValueError(f'Bad confId values: {confId}')
            conf_idxs = [self.GetMinEnergyConfId(int(idx))]
        else:
            conf_idxs = [confId]
        text = ''
        for conf_id in conf_idxs:
            text += self._ConfToXYZ(conf_id)
        
        return text
    
    
    def ToXYZ(self, path, confId = -2):
        '''
        Saves found conformers in XYZ format
        If confId is 'min' or -1, conformer with the lowest energy will be saved
        If confId is in 'min(i)', i > 0 format, conformer with the i-th lowest
        energy will be saved
        If confId is 'all' or -2, all conformers will be saved
        '''
        text = self.ToXYZBlock(confId)
        with open(path, 'w') as outf:
            outf.write(text)


