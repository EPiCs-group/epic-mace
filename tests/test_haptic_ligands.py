"""Tests for haptic (eta-type) and sigma ligand support.

Four complex types are covered:

1.  W(CO)5(eta2-H2)  [Kubas-type sigma-H2 complex]
    An octahedral W(0) complex with a side-on coordinated H2 molecule.
    Simplest possible haptic ligand (eta2, two haptic atoms), used to verify
    the centroid detection and force-field constraint logic end-to-end.
    CO-only ancillary ligands are used to avoid H-H vdW conflicts that arise
    when PH3/NH3 ligands fill all six octahedral positions simultaneously.

2.  (eta5-Cp)Mn(CO)3  [cymantrene]
    A half-sandwich Mn(I) complex.  The Cp centroid [*:1] is bonded to one
    ring carbon in the molecular graph, but MACE now expands that anchor to
    the full five-membered Cp ring internally. Tests verify that the centroid
    is detected as an eta5 haptic donor, that the non-haptic CO donors are
    unaffected, and that 3D embedding succeeds.

3.  cis-[PtCl2(NH3)2]  [cisplatin]
    A classical square-planar Pt(II) complex with no haptic ligands.
    Regression test: sigma-donor atoms (N, Cl) must NOT be classified as
    haptic centroids, two stereomers (cis + trans) must be found, and the
    embedded geometry must match expected bond lengths and angles.

4.  [PtCl2(NH3)(eta2-C2H4)]  [Zeise's-salt analogue]
    A square-planar Pt(II) complex with one pi-bound ethylene ligand, one
    NH3, and two Cl- ligands.  Tests verify that the ethylene centroid is
    detected as an eta2 haptic donor (two C atoms), that sigma donors are
    not misclassified, that 3D embedding succeeds, and that exactly two SP
    stereomers are found.

SMILES conventions
------------------
CO ligand   : [C-:N]#[O+]   — C-donor, dative bond; formal charges maintain
                               the correct Lewis structure (N = position label)
sigma-H2    : [*:N]([H])[H] — centroid bonded to two H atoms
Cp ring     : [C-]1C=CC=C1  — non-aromatic cyclopentadienyl anion;
                               the [C-] carbon (no H) bonds to the centroid
NH3 ligand  : [NH3:N]
Cl ligand   : [Cl-:N]

NOTE: eta5 ring autodetection
-----------------------------
Connecting a centroid to all five ring atoms requires multi-centre bond
notation that is not directly available in SMILES. The supported workflow is
to bond the centroid to one ring carbon ([C-] or aromatic c); MACE then
expands that anchor internally to the full Cp/Cp* ring and treats it as eta5.

Running
-------
    pytest tests/test_haptic_ligands.py -v
"""

import numpy as np
import pytest

from mace import Complex, ComplexFromLigands


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _count_xyz_element(xyz_block, symbol):
    """Count element lines in an XYZ block."""
    return sum(
        1 for line in xyz_block.splitlines()[2:]
        if line.split() and line.split()[0] == symbol
    )

def _dist(conf, idx_a, idx_b):
    """Euclidean distance in Angstrom between two atoms in a conformer."""
    a = np.array(conf.GetAtomPosition(idx_a))
    b = np.array(conf.GetAtomPosition(idx_b))
    return float(np.linalg.norm(a - b))


def _angle_deg(conf, idx_a, idx_ca, idx_b):
    """L-M-L angle in degrees (M = central atom at idx_ca)."""
    ca = np.array(conf.GetAtomPosition(idx_ca))
    va = np.array(conf.GetAtomPosition(idx_a)) - ca
    vb = np.array(conf.GetAtomPosition(idx_b)) - ca
    cos = np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


# ---------------------------------------------------------------------------
# W(CO)3(PH3)2(eta2-H2)  —  Kubas-type sigma-H2 complex
# ---------------------------------------------------------------------------

class TestSigmaH2:
    """W(CO)5(eta2-H2): simplest haptic ligand, eta2 (two H atoms).

    Electron count: W(0) 6d + eta2-H2 2e + 5 CO 10e = 18e.
    The sigma-H2 centroid occupies position 1; CO at 2-6.
    All six OH positions are filled — no dummy helpers needed.

    Note: CO-only ancillary ligands avoid the H-H vdW conflict that prevents
    the RDKit distance-geometry embedder from converging when H-bearing ligands
    (PH3, NH3) simultaneously occupy all six coordination sites.
    """

    SMILES = (
        "[W]"
        "(<-[C-:2]#[O+])"
        "(<-[C-:3]#[O+])"
        "(<-[C-:4]#[O+])"
        "(<-[C-:5]#[O+])"
        "(<-[C-:6]#[O+])"
        "<-[*:1]([H])[H]"
    )
    GEOM = "OH"

    @pytest.fixture(autouse=True)
    def complex_obj(self):
        self.c = Complex(self.SMILES, self.GEOM)

    # --- initialisation ---

    def test_initialisation_succeeds(self):
        assert self.c.err_init is None

    def test_six_donor_atoms(self):
        """All six OH positions must be occupied."""
        assert len(self.c._DAs) == 6

    def test_donor_map_numbers(self):
        assert set(self.c._DAs.values()) == {1, 2, 3, 4, 5, 6}

    # --- haptic-DA detection ---

    def test_one_haptic_da_detected(self):
        """The H2 centroid must be the only haptic DA."""
        assert len(self.c._haptic_DAs) == 1

    def test_hapticity_is_2(self):
        """Two H atoms must be listed as haptic neighbours of the centroid."""
        haptic_atoms = list(self.c._haptic_DAs.values())[0]
        assert len(haptic_atoms) == 2

    def test_haptic_atoms_are_hydrogen(self):
        """The haptic neighbours of the centroid must be H (atomic number 1)."""
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        for idx in haptic_idxs:
            atom = self.c.mol.GetAtomWithIdx(idx)
            assert atom.GetAtomicNum() == 1, (
                f"Haptic atom at idx {idx} is not H (anum={atom.GetAtomicNum()})"
            )

    def test_centroid_is_dummy(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomicNum() == 0

    def test_centroid_map_number(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomMapNum() == 1

    # --- 3D embedding ---

    def test_conformer_generation_succeeds(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0, "Conformer generation failed (returned -1)"

    def test_conformer_energy_is_finite(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        assert np.isfinite(self.c.GetConfEnergy(flag))

    def test_w_centroid_distance(self):
        """W–centroid distance should be near HapticDist[2] = 1.85 A."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        # CA is moved to origin after embedding
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        dist = _dist(conf, self.c._idx_CA, centroid_idx)
        expected = Complex._HapticDist[2]
        assert abs(dist - expected) < 0.5, (
            f"W-centroid distance {dist:.2f} A differs from HapticDist[2]={expected} A by more than 0.5 A"
        )

    def test_centroid_h_distance(self):
        """Centroid–H distance should be near HapticCentR[2] = 0.55 A."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        expected = Complex._HapticCentR[2]
        for h_idx in haptic_idxs:
            dist = _dist(conf, centroid_idx, h_idx)
            assert abs(dist - expected) < 0.4, (
                f"Centroid–H distance {dist:.2f} A differs from HapticCentR[2]={expected} A"
            )


# ---------------------------------------------------------------------------
# (eta5-Cp)Mn(CO)3  —  cymantrene
# ---------------------------------------------------------------------------

class TestCymantrene:
    """(eta5-Cp)Mn(CO)3: cymantrene, half-sandwich OH complex.

    Electron count: Mn(I) 6e + eta5-Cp- 6e + 3 CO 6e = 18e.
    The Cp centroid [*:1] occupies position 1; CO ligands at 2, 3, 4.
    Positions 5 and 6 are filled with dummy atoms during 3D embedding.

    The Cp ring is encoded as a non-aromatic cyclopentadienyl anion:
        [C-]1C=CC=C1
    The [C-] carbon carries the -1 charge (no H) and bonds to the centroid.
    MACE detects that anchor as part of a five-membered pi ring and expands
    the centroid donor to the full eta5 Cp haptic set internally.
    """

    SMILES = (
        "[Mn+]"
        "(<-[C-:2]#[O+])"
        "(<-[C-:3]#[O+])"
        "(<-[C-:4]#[O+])"
        "<-[*:1][C-]1C=CC=C1"
    )
    GEOM = "OH"

    @pytest.fixture(autouse=True)
    def complex_obj(self):
        self.c = Complex(self.SMILES, self.GEOM)

    # --- initialisation ---

    def test_initialisation_succeeds(self):
        assert self.c.err_init is None

    def test_geometry_is_oh(self):
        assert self.c.geom == "OH"

    def test_four_donor_atoms(self):
        """One Cp centroid + three CO donors = four DAs total."""
        assert len(self.c._DAs) == 4

    def test_donor_map_numbers(self):
        assert set(self.c._DAs.values()) == {1, 2, 3, 4}

    # --- haptic-DA detection ---

    def test_one_haptic_da_detected(self):
        """Exactly the Cp centroid dummy must be classified as haptic."""
        assert len(self.c._haptic_DAs) == 1

    def test_centroid_is_dummy(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomicNum() == 0

    def test_centroid_map_number_is_1(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomMapNum() == 1

    def test_hapticity_is_5(self):
        """The Cp centroid must expand to all five ring carbons."""
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        assert len(haptic_idxs) == 5

    def test_haptic_atoms_are_ring_carbons(self):
        """All haptic atoms must be carbon atoms from the Cp ring."""
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        for idx in haptic_idxs:
            atom = self.c.mol.GetAtomWithIdx(idx)
            assert atom.GetAtomicNum() == 6, (
                f"Haptic atom expected C (anum=6), got anum={atom.GetAtomicNum()}"
            )

    def test_cp_ring_keeps_five_hydrogens_in_xyz_output(self):
        """Eta5 Cp should export as C5H5 in XYZ output."""
        flag = self.c.AddConformer(maxAttempts = 20)
        assert flag >= 0
        xyz = self.c.ToXYZBlock(flag)
        assert _count_xyz_element(xyz, 'H') == 5

    def test_co_donors_are_not_haptic(self):
        """CO carbon donors (positions 2, 3, 4) must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num in (2, 3, 4):
                assert idx not in haptic_set, (
                    f"CO donor at position {num} was incorrectly classified as haptic"
                )

    def test_co_donors_are_carbon(self):
        for idx, num in self.c._DAs.items():
            if num in (2, 3, 4):
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 6

    # --- 3D embedding ---

    def test_conformer_generation_succeeds(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0, "Conformer generation failed (returned -1)"

    def test_conformer_energy_is_finite(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        assert np.isfinite(self.c.GetConfEnergy(flag))

    def test_mn_centroid_distance_reasonable(self):
        """Mn-centroid distance should be close to the eta5 Cp default."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        dist = _dist(conf, self.c._idx_CA, centroid_idx)
        expected = Complex._HapticDist[5]
        assert abs(dist - expected) < 0.5, (
            f"Mn-centroid distance {dist:.2f} A differs from HapticDist[5]={expected} A by more than 0.5 A"
        )

    # --- stereomer search ---

    def test_stereomers_can_be_generated(self):
        """GetStereomers must complete without error for a haptic complex."""
        stereomers = self.c.GetStereomers(regime="CA", dropEnantiomers=True)
        assert isinstance(stereomers, list)
        assert len(stereomers) >= 1


# ---------------------------------------------------------------------------
# cis-[PtCl2(NH3)2]  —  cisplatin  (regression: classical SP complex)
# ---------------------------------------------------------------------------

class TestCisplatin:
    """cis-[PtCl2(NH3)2]: square-planar Pt(II), no haptic ligands.

    Electron count: Pt(II) d8 + 2 NH3 4e + 2 Cl- 4e = 16e (SP).
    Positions: NH3 at 1, 2 (adjacent = cis); Cl- at 3, 4 (adjacent = cis).
    In the SP geometry positions 1 and 3 are trans (180 deg), as are 2 and 4.

    Regression checks
    -----------------
    * _haptic_DAs must be empty (sigma donors must not be misclassified)
    * GetStereomers must return exactly 2 isomers: cis (cisplatin) + trans (transplatin)
    * Embedded geometry: Pt-Cl ~2.2-2.6 A, Pt-N ~1.9-2.4 A
    * Cl-Pt-Cl and N-Pt-N angles ~90 deg (cis); N(1)-Pt-Cl(3) ~180 deg (trans)
    """

    GEOM    = "SP"
    LIGANDS = ["[NH3:1]", "[NH3:2]", "[Cl-:3]", "[Cl-:4]"]
    CA      = "[Pt+2]"

    @pytest.fixture(autouse=True)
    def complex_obj(self):
        self.c = ComplexFromLigands(self.LIGANDS, self.CA, self.GEOM)

    # --- initialisation ---

    def test_initialisation_succeeds(self):
        assert self.c.err_init is None

    def test_geometry_is_sp(self):
        assert self.c.geom == "SP"

    def test_four_donor_atoms(self):
        assert len(self.c._DAs) == 4

    def test_donor_map_numbers(self):
        assert set(self.c._DAs.values()) == {1, 2, 3, 4}

    # --- regression: no haptic DAs ---

    def test_no_haptic_das(self):
        """Sigma-donor ligands must NOT be classified as haptic centroids."""
        assert len(self.c._haptic_DAs) == 0

    def test_n_donors_are_nitrogen(self):
        """Positions 1 and 2 (NH3) must be N (atomic number 7)."""
        for idx, num in self.c._DAs.items():
            if num in (1, 2):
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 7

    def test_cl_donors_are_chlorine(self):
        """Positions 3 and 4 (Cl-) must be Cl (atomic number 17)."""
        for idx, num in self.c._DAs.items():
            if num in (3, 4):
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 17

    # --- stereomer search ---

    def test_two_stereomers_found(self):
        """[PtCl2(NH3)2] has exactly two SP isomers: cis and trans."""
        stereomers = self.c.GetStereomers(regime="CA", dropEnantiomers=True)
        assert len(stereomers) == 2, (
            f"Expected 2 stereomers (cis + trans), got {len(stereomers)}"
        )

    # --- 3D embedding ---

    def test_conformer_generation_succeeds(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0, "Conformer generation failed (returned -1)"

    def test_conformer_energy_is_finite(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        assert np.isfinite(self.c.GetConfEnergy(flag))

    def test_pt_cl_distance(self):
        """Pt-Cl bond length: 2.0-3.0 A (crystal value ~2.33 A)."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        cl_idxs = [idx for idx, num in self.c._DAs.items() if num in (3, 4)]
        assert len(cl_idxs) == 2
        for idx in cl_idxs:
            dist = _dist(conf, self.c._idx_CA, idx)
            assert 2.0 < dist < 3.0, (
                f"Pt-Cl distance {dist:.2f} A is outside 2.0-3.0 A"
            )

    def test_pt_n_distance(self):
        """Pt-N bond length: 1.7-2.7 A (crystal value ~2.05 A)."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        n_idxs = [idx for idx, num in self.c._DAs.items() if num in (1, 2)]
        assert len(n_idxs) == 2
        for idx in n_idxs:
            dist = _dist(conf, self.c._idx_CA, idx)
            assert 1.7 < dist < 2.7, (
                f"Pt-N distance {dist:.2f} A is outside 1.7-2.7 A"
            )

    def test_cl_pt_cl_angle_cis(self):
        """Cl-Pt-Cl angle must be ~90 deg (cis isomer, positions 3 and 4 are adjacent)."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf  = self.c.mol3D.GetConformer(flag)
        cl_idxs = [idx for idx, num in self.c._DAs.items() if num in (3, 4)]
        assert len(cl_idxs) == 2
        angle = _angle_deg(conf, cl_idxs[0], self.c._idx_CA, cl_idxs[1])
        assert 70.0 < angle < 110.0, (
            f"Cl-Pt-Cl angle {angle:.1f} deg; expected ~90 deg for cis isomer"
        )

    def test_n_pt_n_angle_cis(self):
        """N-Pt-N angle must be ~90 deg (cis isomer, positions 1 and 2 are adjacent)."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf   = self.c.mol3D.GetConformer(flag)
        n_idxs = [idx for idx, num in self.c._DAs.items() if num in (1, 2)]
        assert len(n_idxs) == 2
        angle = _angle_deg(conf, n_idxs[0], self.c._idx_CA, n_idxs[1])
        assert 70.0 < angle < 110.0, (
            f"N-Pt-N angle {angle:.1f} deg; expected ~90 deg for cis isomer"
        )

    def test_n1_pt_cl3_angle_trans(self):
        """N(pos 1)-Pt-Cl(pos 3) angle must be ~180 deg (trans arrangement in SP)."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf   = self.c.mol3D.GetConformer(flag)
        idx_n  = next(idx for idx, num in self.c._DAs.items() if num == 1)
        idx_cl = next(idx for idx, num in self.c._DAs.items() if num == 3)
        angle  = _angle_deg(conf, idx_n, self.c._idx_CA, idx_cl)
        assert 150.0 < angle < 180.0, (
            f"N(1)-Pt-Cl(3) angle {angle:.1f} deg; expected ~180 deg (trans)"
        )


# ---------------------------------------------------------------------------
# [PtCl2(NH3)(eta2-C2H4)]  —  Zeise's-salt analogue  (pi-complex)
# ---------------------------------------------------------------------------

class TestPtEthylene:
    """[PtCl2(NH3)(eta2-C2H4)]: square-planar Pt(II) with one pi-bound ethylene.

    Electron count: Pt(II) d8 + eta2-C2H4 2e + NH3 2e + 2 Cl- 4e = 16e (SP).
    Positions: Cl- at 1, Cl- at 2, NH3 at 3, C2H4-centroid at 4.

    The ethylene ligand is encoded via the centroid dummy atom [*:4] bonded to
    the two alkene carbons:  [*:4]([CH2])[CH2].  This gives hapticity 2 and
    uses HapticDist[2] = 1.85 A for the Pt-centroid distance.

    Two SP stereomers exist:
    * Isomer A: Cl-trans-Cl (positions 1 and 2), NH3-trans-centroid
    * Isomer B: Cl-trans-NH3 (positions 1 and 3), Cl-trans-centroid
    """

    LIGANDS = ["[Cl-:1]", "[Cl-:2]", "[NH3:3]", "[*:4]([CH2])[CH2]"]
    CA      = "[Pt+2]"
    GEOM    = "SP"

    @pytest.fixture(autouse=True)
    def complex_obj(self):
        self.c = ComplexFromLigands(self.LIGANDS, self.CA, self.GEOM)

    # --- initialisation ---

    def test_initialisation_succeeds(self):
        assert self.c.err_init is None

    def test_geometry_is_sp(self):
        assert self.c.geom == "SP"

    def test_four_donor_atoms(self):
        """One ethylene centroid + two Cl- + one NH3 = four DAs total."""
        assert len(self.c._DAs) == 4

    def test_donor_map_numbers(self):
        assert set(self.c._DAs.values()) == {1, 2, 3, 4}

    # --- haptic-DA detection ---

    def test_one_haptic_da_detected(self):
        """Exactly the ethylene centroid must be classified as haptic."""
        assert len(self.c._haptic_DAs) == 1

    def test_hapticity_is_2(self):
        """Two C atoms must be listed as haptic neighbours of the centroid."""
        haptic_atoms = list(self.c._haptic_DAs.values())[0]
        assert len(haptic_atoms) == 2

    def test_haptic_atoms_are_carbon(self):
        """The haptic neighbours of the ethylene centroid must be C (atomic number 6)."""
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        for idx in haptic_idxs:
            atom = self.c.mol.GetAtomWithIdx(idx)
            assert atom.GetAtomicNum() == 6, (
                f"Haptic atom at idx {idx} is not C (anum={atom.GetAtomicNum()})"
            )

    def test_centroid_is_dummy(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomicNum() == 0

    def test_centroid_map_number(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomMapNum() == 4

    # --- regression: no spurious haptic classification ---

    def test_cl_donors_are_not_haptic(self):
        """Cl- ligands at positions 1 and 2 must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num in (1, 2):
                assert idx not in haptic_set, (
                    f"Cl donor at position {num} was incorrectly classified as haptic"
                )

    def test_n_donor_is_not_haptic(self):
        """NH3 at position 3 must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num == 3:
                assert idx not in haptic_set, (
                    f"NH3 donor at position {num} was incorrectly classified as haptic"
                )

    def test_cl_donors_are_chlorine(self):
        """Positions 1 and 2 (Cl-) must be Cl (atomic number 17)."""
        for idx, num in self.c._DAs.items():
            if num in (1, 2):
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 17

    def test_n_donor_is_nitrogen(self):
        """Position 3 (NH3) must be N (atomic number 7)."""
        for idx, num in self.c._DAs.items():
            if num == 3:
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 7

    # --- stereomer search ---

    def test_two_stereomers_found(self):
        """[PtCl2(NH3)(C2H4)] has exactly two SP stereomers."""
        stereomers = self.c.GetStereomers(regime="CA", dropEnantiomers=True)
        assert len(stereomers) == 2, (
            f"Expected 2 stereomers (Cl-trans-Cl and Cl-trans-NH3), got {len(stereomers)}"
        )

    # --- 3D embedding ---

    def test_conformer_generation_succeeds(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0, "Conformer generation failed (returned -1)"

    def test_conformer_energy_is_finite(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        assert np.isfinite(self.c.GetConfEnergy(flag))

    def test_pt_centroid_distance(self):
        """Pt-centroid distance should be near HapticDist[2] = 1.85 A."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        dist = _dist(conf, self.c._idx_CA, centroid_idx)
        expected = Complex._HapticDist[2]
        assert abs(dist - expected) < 0.5, (
            f"Pt-centroid distance {dist:.2f} A differs from HapticDist[2]={expected} A by more than 0.5 A"
        )

    def test_centroid_c_distances(self):
        """Centroid-C distances should be near 0.67 A (half C=C bond)."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        for c_idx in haptic_idxs:
            dist = _dist(conf, centroid_idx, c_idx)
            assert 0.4 < dist < 1.2, (
                f"Centroid-C distance {dist:.2f} A is outside 0.4-1.2 A"
            )

    def test_cc_distance(self):
        """C-C distance should be in a chemically reasonable range (1.0-2.0 A)."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        haptic_idxs = list(self.c._haptic_DAs.values())[0]
        assert len(haptic_idxs) == 2
        dist = _dist(conf, haptic_idxs[0], haptic_idxs[1])
        assert 1.0 < dist < 2.0, (
            f"C-C distance {dist:.2f} A is outside 1.0-2.0 A"
        )

    def test_ethylene_perpendicular_to_coord_plane(self):
        """C=C bond must be roughly perpendicular to the coordination plane.

        The angle between the C1-C2 vector and the Pt-centroid axis should be
        ~90 deg (C atoms on opposite sides of the centroid, perpendicular to
        the M-centroid line).  This ensures the ethylene pi-system faces the
        metal, not the C-H bonds.
        """
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        h1, h2 = self.c._haptic_DAs[centroid_idx]
        pt  = np.array(conf.GetAtomPosition(self.c._idx_CA))
        xc  = np.array(conf.GetAtomPosition(centroid_idx))
        c1  = np.array(conf.GetAtomPosition(h1))
        c2  = np.array(conf.GetAtomPosition(h2))
        v_cc  = c2 - c1
        v_ptx = xc - pt
        cos_a = np.dot(v_cc, v_ptx) / (np.linalg.norm(v_cc) * np.linalg.norm(v_ptx))
        angle = np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0)))
        assert 70.0 < angle < 110.0, (
            f"Angle(C=C, Pt-centroid) = {angle:.1f} deg; expected ~90 deg"
        )


# ---------------------------------------------------------------------------
# Ni(eta3-C3H5)(H)(PMe3)2(CO)2  —  Ni(II) allyl complex
# ---------------------------------------------------------------------------

class TestNiAllylComplex:
    """Ni(eta3-C3H5)(H)(PMe3)2(CO)2: octahedral Ni(II) with eta3-allyl.

    Electron count: Ni(II) d8 + eta3-allyl- 4e + H- 2e + 2 PMe3 4e + 2 CO 4e = 18e.
    Positions: allyl centroid at 1; hydride at 2 (trans to allyl); PMe3 at 3, 4;
    CO at 5, 6.

    The allyl ligand is encoded with all three C atoms bonded directly to the
    centroid dummy atom: [*:1]([CH2])([CH])[CH2].  This gives hapticity 3 and
    uses HapticDist[3] = 2.00 A for the Ni-centroid distance and
    HapticCentR[3] = 1.20 A for the centroid-C distances.

    Note on SMILES encoding
    -----------------------
    The documented eta3-allyl example "[*:N]([CH2])[CH][CH2]" bonds the centroid
    to only the first CH2 and the central CH (2 direct neighbours => eta2 in the
    molecular graph).  To obtain true eta3 detection, all three allyl carbons
    must be in direct branches of the centroid:
        [*:N]([CH2])([CH])[CH2]
    No explicit C-C bonds appear in the molecular graph; the allyl geometry is
    maintained entirely through the centroid-C distance constraints.
    """

    LIGANDS = [
        "[*:1]([CH2])([CH])[CH2]",  # eta3-allyl centroid at position 1
        "[H-:2]",                   # hydride at position 2
        "[P:3](C)(C)C",             # PMe3 at position 3
        "[P:4](C)(C)C",             # PMe3 at position 4
        "[C-:5]#[O+]",              # CO at position 5
        "[C-:6]#[O+]",              # CO at position 6
    ]
    CA   = "[Ni+2]"
    GEOM = "OH"

    @pytest.fixture(autouse=True)
    def complex_obj(self):
        self.c = ComplexFromLigands(self.LIGANDS, self.CA, self.GEOM)

    # --- initialisation ---

    def test_initialisation_succeeds(self):
        assert self.c.err_init is None

    def test_geometry_is_oh(self):
        assert self.c.geom == "OH"

    def test_six_donor_atoms(self):
        """One allyl centroid + hydride + 2 PMe3 + 2 CO = six DAs total."""
        assert len(self.c._DAs) == 6

    def test_donor_map_numbers(self):
        assert set(self.c._DAs.values()) == {1, 2, 3, 4, 5, 6}

    # --- haptic-DA detection ---

    def test_one_haptic_da_detected(self):
        """Exactly the allyl centroid must be classified as haptic."""
        assert len(self.c._haptic_DAs) == 1

    def test_hapticity_is_3(self):
        """Three C atoms must be listed as haptic neighbours of the centroid."""
        haptic_atoms = list(self.c._haptic_DAs.values())[0]
        assert len(haptic_atoms) == 3

    def test_haptic_atoms_are_carbon(self):
        """All three haptic neighbours of the allyl centroid must be C (anum 6)."""
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        for idx in haptic_idxs:
            atom = self.c.mol.GetAtomWithIdx(idx)
            assert atom.GetAtomicNum() == 6, (
                f"Haptic atom at idx {idx} is not C (anum={atom.GetAtomicNum()})"
            )

    def test_centroid_is_dummy(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomicNum() == 0

    def test_centroid_map_number_is_1(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomMapNum() == 1

    # --- regression: no spurious haptic classification ---

    def test_hydride_is_not_haptic(self):
        """H- at position 2 must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num == 2:
                assert idx not in haptic_set, (
                    f"Hydride donor at position {num} was incorrectly classified as haptic"
                )

    def test_p_donors_are_not_haptic(self):
        """PMe3 donors at positions 3 and 4 must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num in (3, 4):
                assert idx not in haptic_set, (
                    f"P donor at position {num} was incorrectly classified as haptic"
                )

    def test_co_donors_are_not_haptic(self):
        """CO donors at positions 5 and 6 must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num in (5, 6):
                assert idx not in haptic_set, (
                    f"CO donor at position {num} was incorrectly classified as haptic"
                )

    def test_hydride_is_hydrogen(self):
        """Position 2 (H-) must be H (atomic number 1)."""
        for idx, num in self.c._DAs.items():
            if num == 2:
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 1

    def test_p_donors_are_phosphorus(self):
        """Positions 3 and 4 (PMe3) must be P (atomic number 15)."""
        for idx, num in self.c._DAs.items():
            if num in (3, 4):
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 15

    def test_co_donors_are_carbon(self):
        """Positions 5 and 6 (CO) must be C (atomic number 6)."""
        for idx, num in self.c._DAs.items():
            if num in (5, 6):
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 6

    # --- stereomer search ---

    def test_stereomers_can_be_generated(self):
        """GetStereomers must complete without error for an eta3-allyl complex."""
        stereomers = self.c.GetStereomers(regime="CA", dropEnantiomers=True)
        assert isinstance(stereomers, list)
        assert len(stereomers) >= 1

    # --- 3D embedding ---

    def test_conformer_generation_succeeds(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0, "Conformer generation failed (returned -1)"

    def test_conformer_energy_is_finite(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        assert np.isfinite(self.c.GetConfEnergy(flag))

    def test_ni_centroid_distance(self):
        """Ni-centroid distance should be near HapticDist[3] = 2.00 A."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        dist = _dist(conf, self.c._idx_CA, centroid_idx)
        expected = Complex._HapticDist[3]
        assert abs(dist - expected) < 0.5, (
            f"Ni-centroid distance {dist:.2f} A differs from HapticDist[3]={expected} A by >0.5 A"
        )

    def test_centroid_c_distances(self):
        """Centroid-C distances should be near HapticCentR[3] = 1.20 A."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        expected = Complex._HapticCentR[3]
        for c_idx in haptic_idxs:
            dist = _dist(conf, centroid_idx, c_idx)
            assert abs(dist - expected) < 0.5, (
                f"Centroid-C distance {dist:.2f} A differs from HapticCentR[3]={expected} A by >0.5 A"
            )

    def test_ni_h_distance(self):
        """Ni-H distance should be in the expected range for a terminal hydride."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        h_idx = next(idx for idx, num in self.c._DAs.items() if num == 2)
        dist = _dist(conf, self.c._idx_CA, h_idx)
        assert 1.2 < dist < 2.5, (
            f"Ni-H distance {dist:.2f} A is outside the expected range 1.2-2.5 A"
        )


# ---------------------------------------------------------------------------
# [Cp*(Mn)(NH3)(CO)2]  —  Mn(I) piano-stool complex
# ---------------------------------------------------------------------------

class TestCpStarMnComplex:
    """CpStar(Mn)(NH3)(CO)2: Mn(I) half-sandwich (piano-stool) complex.

    Electron count: Mn(0/I) d7/d6 + eta5-Cp* 6e + NH3 2e + 2 CO 4e = 18e.
    The Cp* centroid [*:1] occupies position 1 in the OH framework; NH3 at 2;
    CO at 3 and 4.  Positions 5 and 6 are filled with dummy atoms during 3D
    embedding (as for cymantrene).

    Cp* SMILES encoding
    -------------------
    Cp* (pentamethylcyclopentadienyl) is encoded as a 5-methyl aromatic ring:
        [*:1]c1(C)c(C)c(C)c(C)c1C
    All five ring carbons carry a methyl group.  The centroid bonds to one ring
    carbon; MACE detects the full aromatic five-membered ring and treats the
    centroid as an eta5 Cp* donor internally.
    """

    LIGANDS = [
        "[*:1]c1(C)c(C)c(C)c(C)c1C",  # eta5-Cp* centroid at position 1
        "[NH3:2]",                       # NH3 at position 2
        "[C-:3]#[O+]",                   # CO at position 3
        "[C-:4]#[O+]",                   # CO at position 4
    ]
    CA   = "[Mn]"
    GEOM = "OH"

    @pytest.fixture(autouse=True)
    def complex_obj(self):
        self.c = ComplexFromLigands(self.LIGANDS, self.CA, self.GEOM)

    # --- initialisation ---

    def test_initialisation_succeeds(self):
        assert self.c.err_init is None

    def test_geometry_is_oh(self):
        assert self.c.geom == "OH"

    def test_four_donor_atoms(self):
        """Cp* centroid + NH3 + 2 CO = four DAs total."""
        assert len(self.c._DAs) == 4

    def test_donor_map_numbers(self):
        assert set(self.c._DAs.values()) == {1, 2, 3, 4}

    # --- haptic-DA detection ---

    def test_one_haptic_da_detected(self):
        """Exactly the Cp* centroid must be classified as haptic."""
        assert len(self.c._haptic_DAs) == 1

    def test_centroid_is_dummy(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomicNum() == 0

    def test_centroid_map_number_is_1(self):
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        assert self.c.mol.GetAtomWithIdx(centroid_idx).GetAtomMapNum() == 1

    def test_hapticity_is_5(self):
        """The Cp* centroid must expand to all five ring carbons."""
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        assert len(haptic_idxs) == 5

    def test_haptic_atoms_are_ring_carbons(self):
        """All haptic atoms in the expanded Cp* set must be carbon."""
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        haptic_idxs  = self.c._haptic_DAs[centroid_idx]
        for idx in haptic_idxs:
            atom = self.c.mol.GetAtomWithIdx(idx)
            assert atom.GetAtomicNum() == 6, (
                f"Haptic atom expected C (anum=6), got anum={atom.GetAtomicNum()}"
            )

    def test_cpstar_retains_all_methyl_hydrogens_in_xyz_output(self):
        """Cp*Mn(NH3)(CO)2 should export all 18 hydrogens in XYZ output."""
        flag = self.c.AddConformer(maxAttempts = 20)
        assert flag >= 0
        xyz = self.c.ToXYZBlock(flag)
        assert _count_xyz_element(xyz, 'H') == 18

    # --- regression: no spurious haptic classification ---

    def test_nh3_is_not_haptic(self):
        """NH3 at position 2 must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num == 2:
                assert idx not in haptic_set, (
                    f"NH3 donor at position {num} incorrectly classified as haptic"
                )

    def test_co_donors_are_not_haptic(self):
        """CO donors at positions 3 and 4 must NOT be in _haptic_DAs."""
        haptic_set = set(self.c._haptic_DAs.keys())
        for idx, num in self.c._DAs.items():
            if num in (3, 4):
                assert idx not in haptic_set, (
                    f"CO donor at position {num} incorrectly classified as haptic"
                )

    def test_nh3_donor_is_nitrogen(self):
        """Position 2 (NH3) must be N (atomic number 7)."""
        for idx, num in self.c._DAs.items():
            if num == 2:
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 7

    def test_co_donors_are_carbon(self):
        """Positions 3 and 4 (CO) must be C (atomic number 6)."""
        for idx, num in self.c._DAs.items():
            if num in (3, 4):
                assert self.c.mol.GetAtomWithIdx(idx).GetAtomicNum() == 6

    # --- stereomer search ---

    def test_stereomers_can_be_generated(self):
        """GetStereomers must complete without error for a piano-stool complex."""
        stereomers = self.c.GetStereomers(regime="CA", dropEnantiomers=True)
        assert isinstance(stereomers, list)
        assert len(stereomers) >= 1

    # --- 3D embedding ---

    def test_conformer_generation_succeeds(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0, "Conformer generation failed (returned -1)"

    def test_conformer_energy_is_finite(self):
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        assert np.isfinite(self.c.GetConfEnergy(flag))

    def test_mn_centroid_distance_reasonable(self):
        """Mn-centroid distance should be close to the eta5 Cp default."""
        flag = self.c.AddConformer(maxAttempts=20)
        assert flag >= 0
        conf = self.c.mol3D.GetConformer(flag)
        centroid_idx = list(self.c._haptic_DAs.keys())[0]
        dist = _dist(conf, self.c._idx_CA, centroid_idx)
        expected = Complex._HapticDist[5]
        assert abs(dist - expected) < 0.5, (
            f"Mn-centroid distance {dist:.2f} A differs from HapticDist[5]={expected} A by more than 0.5 A"
        )
