"""Tests for sigma (η²-H₂) and pi (η⁵-Cp) Ru complexes with SNS pincer variants.

Three target complexes
----------------------

1.  Ru(SNS-deprot)(η²-H₂)(H⁻)(PPh₃)
    ─ SNS pincer, deprotonated N (amide N⁻), tridentate S/N/S
    ─ one side-on σ-H₂ ligand (Kubas type, η²)
    ─ one classic hydride (H⁻)
    ─ one PPh₃
    ─ Geometry: octahedral (OH), Ru(II), neutral complex

2.  Ru(NS-deprot)(η²-H₂)₂(H⁻)(PPh₃)
    ─ SNS pincer used as bidentate (one S + deprotonated N); second S arm is pendant
    ─ two side-on σ-H₂ ligands (η²)
    ─ one classic hydride (H⁻)
    ─ one PPh₃
    ─ Geometry: OH, Ru(II), neutral complex

3.  Ru(NS-prot)(η⁵-Cp)(H⁻)
    ─ SNS pincer used as bidentate (one S + protonated N-H); second S arm is pendant
    ─ one η⁵-cyclopentadienyl ring (Cp⁻, pi-ligand)
    ─ one classic hydride (H⁻)
    ─ Geometry: OH (half-sandwich; two phantom positions filled automatically)
    ─ Ru(II), neutral complex

SMILES conventions
------------------
    deprotonated N  :  [N-:n]      (amide, no H)
    protonated N    :  [NH:n]      (amine, explicit H kept through dative bond formation)
    free S arm      :  -S-         (no map number  →  not tagged as donor)
    sigma-H2        :  [*:n]([H])[H]    (centroid bonded to two H atoms)
    eta5-Cp         :  [*:n][C-]1C=CC=C1   (centroid → one Cp carbon; ring auto-expanded to η⁵)
    hydride         :  [H-:n]
    PPh3            :  [P:n](c1ccccc1)(c1ccccc1)c1ccccc1

Running
-------
    pytest tests/test_sigma_pi_complexes.py -v
    pytest tests/test_sigma_pi_complexes.py -v --tb=short   # for concise output

Performance note
----------------
    Each test class uses ``setup_class`` so that the expensive conformer-generation
    step is performed only **once** per class (not before every individual test).
"""

import sys
from pathlib import Path

# Ensure the local development version of mace is used, not the pip-installed one.
# The haptic-ligand features tested here are present in the source repo but are
# not yet part of the published 0.5.0 package.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import numpy as np
import pytest

from mace import Complex, ComplexFromLigands


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _dist(conf, idx_a, idx_b):
    a = np.array(conf.GetAtomPosition(idx_a))
    b = np.array(conf.GetAtomPosition(idx_b))
    return float(np.linalg.norm(a - b))


def _lowest_energy_isomer(stereomers):
    """Return the stereomer with the lowest energy conformer."""
    best = None
    best_e = float('inf')
    for X in stereomers:
        if X.GetNumConformers() > 0:
            e = X.GetConfEnergy(0)
            if e < best_e:
                best_e = e
                best = X
    return best, best_e


# ---------------------------------------------------------------------------
# Complex 1 — Ru(SNS-deprot)(η²-H₂)(H⁻)(PPh₃)
# ---------------------------------------------------------------------------

class TestRuSNS_Deprot_SigmaH2_H_PPh3:
    """Ru(II) octahedral: SNS tridentate (N deprotonated) + σ-H₂ + H⁻ + PPh₃.

    Donor atoms  : S(1) + N⁻(2) + S(3) [SNS tridentate]  +  centroid(4) [σ-H₂]
                   + H⁻(5)  +  P(6) [PPh₃]  =  6 total (OH)
    Charge       : Ru²⁺ + N⁻(−1) + H⁻(−1) = neutral
    Haptic ligand: η²-H₂  →  centroid bonded to 2 H atoms
    """

    LIGANDS = [
        'CC[S:1]CC[N-:2]CC[S:3]CC',              # SNS tridentate, N deprotonated
        '[*:4]([H])[H]',                           # sigma-H2 (eta2 centroid)
        '[H-:5]',                                  # hydride
        '[P:6](c1ccccc1)(c1ccccc1)c1ccccc1',      # PPh3
    ]
    CA   = '[Ru+2]'
    GEOM = 'OH'

    @classmethod
    def setup_class(cls):
        """Build the complex once for the whole class (avoids repeated conformer generation)."""
        cls.X0 = ComplexFromLigands(cls.LIGANDS, cls.CA, cls.GEOM)
        cls.stereomers = cls.X0.GetStereomers(
            regime='all', dropEnantiomers=True, merRule=True
        )
        for X in cls.stereomers:
            X.AddConformers(numConfs=10, rmsThresh=0.5)
            X.OrderConfsByEnergy()

    # --- initialisation ---

    def test_init_succeeds(self):
        # err_init is None means "no error" – complex assembled successfully
        assert self.X0.err_init is None
        # each resolved stereomer must also have no error
        for X in self.stereomers:
            assert X.err_init is None, f"Stereomer has err_init: {X.err_init}"

    def test_six_donor_atoms(self):
        assert len(self.X0._DAs) == 6

    def test_geometry_is_oh(self):
        assert self.X0.geom == 'OH'

    # --- haptic detection ---

    def test_one_haptic_da(self):
        """Exactly one haptic DA: the σ-H₂ centroid."""
        assert len(self.X0._haptic_DAs) == 1

    def test_hapticity_is_2(self):
        haptic_atoms = list(self.X0._haptic_DAs.values())[0]
        assert len(haptic_atoms) == 2

    def test_haptic_atoms_are_hydrogen(self):
        centroid_idx  = list(self.X0._haptic_DAs.keys())[0]
        haptic_idxs   = self.X0._haptic_DAs[centroid_idx]
        for idx in haptic_idxs:
            assert self.X0.mol.GetAtomWithIdx(idx).GetAtomicNum() == 1

    # --- deprotonated N ---

    def test_n_has_no_hydrogen(self):
        """Amide N⁻ must carry no hydrogen after complex assembly."""
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'N' and a.GetAtomMapNum() > 0:
                assert a.GetTotalNumHs() == 0, (
                    f"N donor has {a.GetTotalNumHs()} H(s); expected 0 for deprotonated N"
                )

    def test_n_is_anionic(self):
        """Amide N must carry formal charge −1."""
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'N' and a.GetAtomMapNum() > 0:
                assert a.GetFormalCharge() == -1

    # --- stereomers ---

    def test_multiple_stereomers_found(self):
        assert len(self.stereomers) >= 3, (
            f"Expected ≥3 stereomers, found {len(self.stereomers)}"
        )

    # --- 3D embedding ---

    def test_conformers_generated(self):
        n_with_confs = sum(1 for X in self.stereomers if X.GetNumConformers() > 0)
        assert n_with_confs >= 3, (
            f"Only {n_with_confs}/{len(self.stereomers)} isomers have conformers"
        )

    def test_best_energy_is_finite(self):
        _, e = _lowest_energy_isomer(self.stereomers)
        assert np.isfinite(e)

    def test_xyz_block_contains_hydrogen(self):
        """XYZ output must include the two σ-H₂ hydrogen atoms."""
        best, _ = _lowest_energy_isomer(self.stereomers)
        assert best is not None
        xyz = best.ToXYZBlock(confId=0)
        h_lines = [l for l in xyz.splitlines()[2:] if l.split() and l.split()[0] == 'H']
        # at minimum: 2 from sigma-H2 + 1 hydride + many from SNS/PPh3
        assert len(h_lines) >= 3

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Known UFF limitation: after force-field relaxation the centroid "
            "dummy atom is pulled to a bond-length distance (~0.73 Å) instead "
            "of staying at the target HapticDist[2] = 1.85 Å.  The graph "
            "topology is correct but the 3D geometry is not optimal for σ-H₂."
        ),
    )
    def test_ru_centroid_distance(self):
        """Ru–σ-H₂ centroid distance should be near HapticDist[2] = 1.85 Å.

        NOTE: This test is expected to FAIL with the current UFF parameterisation.
        The centroid is placed at ~0.73 Å from Ru rather than the target 1.85 Å.
        Mark strict=True so that an unexpected pass will be surfaced as a test error,
        signalling that the geometry has been improved.
        """
        best, _ = _lowest_energy_isomer(self.stereomers)
        assert best is not None
        conf = best.mol3D.GetConformer(0)
        centroid_idx = list(best._haptic_DAs.keys())[0]
        d = _dist(conf, best._idx_CA, centroid_idx)
        expected = Complex._HapticDist[2]
        assert abs(d - expected) < 0.6, (
            f"Ru–centroid distance {d:.2f} Å deviates more than 0.6 Å from {expected} Å"
        )


# ---------------------------------------------------------------------------
# Complex 2 — Ru(NS-deprot)(η²-H₂)₂(H⁻)(PPh₃)
# ---------------------------------------------------------------------------

class TestRuNS_Deprot_TwoSigmaH2_H_PPh3:
    """Ru(II) octahedral: SNS as bidentate NS (deprotonated N) + 2×σ-H₂ + H⁻ + PPh₃.

    The SNS pincer coordinates through one S and the deprotonated N only;
    the second S arm is pendant (no map number → not a donor).

    Donor atoms  : S(1) + N⁻(2) [NS bidentate]  +  centroid(3) [σ-H₂ #1]
                   + centroid(4) [σ-H₂ #2]  +  H⁻(5)  +  P(6) [PPh₃]  =  6 total (OH)
    Charge       : Ru²⁺ + N⁻(−1) + H⁻(−1) = neutral
    Haptic ligands: 2 × η²-H₂
    """

    LIGANDS = [
        'CC[S:1]CC[N-:2]CCSCC',                   # NS bidentate (second S arm free)
        '[*:3]([H])[H]',                            # sigma-H2 #1
        '[*:4]([H])[H]',                            # sigma-H2 #2
        '[H-:5]',                                   # hydride
        '[P:6](c1ccccc1)(c1ccccc1)c1ccccc1',       # PPh3
    ]
    CA   = '[Ru+2]'
    GEOM = 'OH'

    @classmethod
    def setup_class(cls):
        """Build the complex once for the whole class."""
        cls.X0 = ComplexFromLigands(cls.LIGANDS, cls.CA, cls.GEOM)
        cls.stereomers = cls.X0.GetStereomers(
            regime='all', dropEnantiomers=True, merRule=True
        )
        for X in cls.stereomers:
            X.AddConformers(numConfs=10, rmsThresh=0.5)
            X.OrderConfsByEnergy()

    # --- initialisation ---

    def test_init_ok(self):
        assert self.X0.err_init is None
        for X in self.stereomers:
            assert X.err_init is None

    def test_six_donor_atoms(self):
        assert len(self.X0._DAs) == 6

    # --- haptic detection ---

    def test_two_haptic_das(self):
        """Both σ-H₂ centroids must be detected as haptic DAs."""
        assert len(self.X0._haptic_DAs) == 2

    def test_both_hapticity_2(self):
        for haptic_idxs in self.X0._haptic_DAs.values():
            assert len(haptic_idxs) == 2

    def test_all_haptic_atoms_are_hydrogen(self):
        for centroid_idx, haptic_idxs in self.X0._haptic_DAs.items():
            for idx in haptic_idxs:
                assert self.X0.mol.GetAtomWithIdx(idx).GetAtomicNum() == 1

    # --- pendant S arm (must NOT be a donor) ---

    def test_second_s_is_not_donor(self):
        """The second (untagged) S must not appear in _DAs."""
        donor_idxs = set(self.X0._DAs.keys())
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'S' and a.GetAtomMapNum() == 0:
                assert a.GetIdx() not in donor_idxs, (
                    "Untagged S atom incorrectly registered as donor"
                )

    # --- deprotonated N ---

    def test_n_has_no_hydrogen(self):
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'N' and a.GetAtomMapNum() > 0:
                assert a.GetTotalNumHs() == 0

    def test_n_is_anionic(self):
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'N' and a.GetAtomMapNum() > 0:
                assert a.GetFormalCharge() == -1

    # --- stereomers ---

    def test_multiple_stereomers(self):
        assert len(self.stereomers) >= 3

    # --- 3D embedding ---

    def test_conformers_generated(self):
        n_ok = sum(1 for X in self.stereomers if X.GetNumConformers() > 0)
        assert n_ok >= 3

    def test_best_energy_finite(self):
        _, e = _lowest_energy_isomer(self.stereomers)
        assert np.isfinite(e)

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Known UFF limitation: after force-field relaxation each centroid "
            "dummy atom collapses to ~0.69 Å from Ru rather than the target "
            "HapticDist[2] = 1.85 Å.  Graph topology is correct."
        ),
    )
    def test_both_centroids_near_ru(self):
        """Both σ-H₂ centroids should be within 1.85 ± 0.6 Å of Ru.

        NOTE: Expected to FAIL with current UFF; marked xfail(strict=True).
        """
        best, _ = _lowest_energy_isomer(self.stereomers)
        assert best is not None
        conf = best.mol3D.GetConformer(0)
        expected = Complex._HapticDist[2]
        for centroid_idx in best._haptic_DAs:
            d = _dist(conf, best._idx_CA, centroid_idx)
            assert abs(d - expected) < 0.6, (
                f"Ru–centroid distance {d:.2f} Å deviates more than 0.6 Å from {expected} Å"
            )


# ---------------------------------------------------------------------------
# Complex 3 — Ru(NS-prot)(η⁵-Cp)(H⁻)
# ---------------------------------------------------------------------------

class TestRuNS_Prot_Cp_H:
    """Ru(II) half-sandwich: SNS as bidentate NS (protonated N-H) + η⁵-Cp + H⁻.

    The SNS pincer coordinates through one S and the protonated N-H only;
    the second S arm is pendant.  The η⁵-Cp ring occupies one OH position;
    mace fills the remaining two phantom positions automatically during embedding.

    Donor atoms  : S(1) + NH(2) [NS bidentate]  +  centroid(3) [η⁵-Cp]
                   + H⁻(4)  =  4 total (OH, with 2 phantom positions)
    Charge       : Ru²⁺ + Cp⁻(−1) + H⁻(−1) = neutral
    Haptic ligand: η⁵-Cp → centroid auto-expanded to 5-membered ring
    """

    LIGANDS = [
        'CC[S:1]CC[NH:2]CCSCC',         # NS bidentate (N-H protonated, second S free)
        '[*:3][C-]1C=CC=C1',             # eta5-Cp (centroid -> ring anchor; auto-expanded)
        '[H-:4]',                         # hydride
    ]
    CA   = '[Ru+2]'
    GEOM = 'OH'

    @classmethod
    def setup_class(cls):
        """Build the complex once for the whole class."""
        cls.X0 = ComplexFromLigands(cls.LIGANDS, cls.CA, cls.GEOM)
        cls.stereomers = cls.X0.GetStereomers(
            regime='all', dropEnantiomers=True, merRule=True
        )
        for X in cls.stereomers:
            X.AddConformers(numConfs=10, rmsThresh=0.5)
            X.OrderConfsByEnergy()

    # --- initialisation ---

    def test_init_ok(self):
        assert self.X0.err_init is None
        for X in self.stereomers:
            assert X.err_init is None

    def test_four_donor_atoms(self):
        """S + NH + Cp-centroid + H⁻ = 4 explicit DAs (2 phantom positions added by mace)."""
        assert len(self.X0._DAs) == 4

    # --- haptic detection ---

    def test_one_haptic_da(self):
        """Exactly one haptic DA: the η⁵-Cp centroid."""
        assert len(self.X0._haptic_DAs) == 1

    def test_hapticity_is_5(self):
        """Cp ring auto-expansion must yield 5 haptic atoms."""
        haptic_atoms = list(self.X0._haptic_DAs.values())[0]
        assert len(haptic_atoms) == 5

    def test_haptic_atoms_are_carbon(self):
        """All five haptic atoms of the Cp ring must be carbon."""
        centroid_idx = list(self.X0._haptic_DAs.keys())[0]
        haptic_idxs  = self.X0._haptic_DAs[centroid_idx]
        for idx in haptic_idxs:
            assert self.X0.mol.GetAtomWithIdx(idx).GetAtomicNum() == 6

    # --- protonated N ---

    def test_n_has_one_hydrogen(self):
        """Protonated amine N must carry exactly 1 hydrogen."""
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'N' and a.GetAtomMapNum() > 0:
                assert a.GetTotalNumHs() == 1, (
                    f"N donor has {a.GetTotalNumHs()} H(s); expected 1 for protonated N-H"
                )

    def test_n_is_neutral(self):
        """Protonated N must be neutral (not anionic)."""
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'N' and a.GetAtomMapNum() > 0:
                assert a.GetFormalCharge() == 0

    # --- pendant S arm ---

    def test_second_s_is_not_donor(self):
        donor_idxs = set(self.X0._DAs.keys())
        for a in self.X0.mol.GetAtoms():
            if a.GetSymbol() == 'S' and a.GetAtomMapNum() == 0:
                assert a.GetIdx() not in donor_idxs

    # --- stereomers ---

    def test_multiple_stereomers(self):
        assert len(self.stereomers) >= 3

    # --- 3D embedding ---

    def test_conformers_generated(self):
        n_ok = sum(1 for X in self.stereomers if X.GetNumConformers() > 0)
        assert n_ok >= 3

    def test_best_energy_finite(self):
        _, e = _lowest_energy_isomer(self.stereomers)
        assert np.isfinite(e)

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Known UFF limitation: after force-field relaxation the η⁵-Cp "
            "centroid collapses to ~0.78 Å from Ru rather than the target "
            "HapticDist[5] = 1.90 Å.  Graph topology (5 haptic C atoms) is "
            "correct; only the optimised bond geometry is affected."
        ),
    )
    def test_ru_cp_centroid_distance(self):
        """Ru–Cp centroid distance should be near HapticDist[5] = 1.90 Å.

        NOTE: Expected to FAIL with current UFF; marked xfail(strict=True).
        """
        best, _ = _lowest_energy_isomer(self.stereomers)
        assert best is not None
        conf = best.mol3D.GetConformer(0)
        centroid_idx = list(best._haptic_DAs.keys())[0]
        d = _dist(conf, best._idx_CA, centroid_idx)
        expected = Complex._HapticDist[5]
        assert abs(d - expected) < 0.6, (
            f"Ru–Cp centroid distance {d:.2f} Å deviates more than 0.6 Å from {expected} Å"
        )

    def test_xyz_contains_cp_carbons(self):
        """XYZ output must include at least 5 carbon atoms (from the Cp ring)."""
        best, _ = _lowest_energy_isomer(self.stereomers)
        assert best is not None
        xyz = best.ToXYZBlock(confId=0)
        c_lines = [l for l in xyz.splitlines()[2:] if l.split() and l.split()[0] == 'C']
        assert len(c_lines) >= 5
