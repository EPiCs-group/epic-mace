"""Regression tests for additional coordination modes."""

import numpy as np

from mace import Complex, ComplexFromLigands


def _count_xyz_element(xyz_block, symbol):
    """Count element lines in an XYZ block."""
    return sum(
        1 for line in xyz_block.splitlines()[2:]
        if line.split() and line.split()[0] == symbol
    )


def _angle_deg(conf, idx_a, idx_ca, idx_b):
    """L-M-L angle in degrees (M = central atom at idx_ca)."""
    ca = np.array(conf.GetAtomPosition(idx_ca))
    va = np.array(conf.GetAtomPosition(idx_a)) - ca
    vb = np.array(conf.GetAtomPosition(idx_b)) - ca
    cos = np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def test_new_geometry_tables_are_registered():
    """The new geometries must expose the expected symmetry group sizes."""
    assert len(Complex._Syms['TET']) == 2
    assert len(Complex._EqOrs['TET'][1]) == 12
    assert 'enantTET' in Complex._EqOrs

    assert len(Complex._Syms['SPY']) == 30
    assert len(Complex._EqOrs['SPY'][1]) == 4
    assert 'enantSPY' in Complex._EqOrs

    assert len(Complex._Syms['TBP']) == 20
    assert len(Complex._EqOrs['TBP'][1]) == 6
    assert 'enantTBP' in Complex._EqOrs

    assert len(Complex._Syms['SAN']) == 1
    assert len(Complex._EqOrs['SAN'][1]) == 2


def test_tetrahedral_enantiomer_pair_is_detected():
    """Four distinct ligands in TET must produce one enantiomeric pair."""
    ligands = ['[F-:1]', '[Cl-:2]', '[Br-:3]', '[I-:4]']
    c = ComplexFromLigands(ligands, '[Zn+2]', 'TET')

    assert c.err_init is None
    assert len(c.GetStereomers(regime='CA', dropEnantiomers=False)) == 2
    assert len(c.GetStereomers(regime='CA', dropEnantiomers=True)) == 1

    flag = c.AddConformer(maxAttempts=20)
    assert flag >= 0

    conf = c.mol3D.GetConformer(flag)
    donor_idxs = [idx for idx, _ in sorted(c._DAs.items(), key=lambda item: item[1])]
    angles = [
        _angle_deg(conf, donor_idxs[i], c._idx_CA, donor_idxs[j])
        for i in range(len(donor_idxs) - 1)
        for j in range(i + 1, len(donor_idxs))
    ]
    for angle in angles:
        assert 95.0 < angle < 125.0, (
            f"Tetrahedral angle {angle:.1f} deg is outside the expected range"
        )


def test_square_pyramidal_stereomer_count():
    """Five distinct ligands in SPY must give 30 arrangements, 15 after dropping enantiomers."""
    ligands = ['[F-:1]', '[Cl-:2]', '[Br-:3]', '[I-:4]', '[NH3:5]']
    c = ComplexFromLigands(ligands, '[Zn+2]', 'SPY')

    assert c.err_init is None
    assert len(c.GetStereomers(regime='CA', dropEnantiomers=False)) == 30
    assert len(c.GetStereomers(regime='CA', dropEnantiomers=True)) == 15
    assert c.AddConformer(maxAttempts=20) >= 0


def test_trigonal_bipyramidal_stereomer_count():
    """Five distinct ligands in TBP must give 20 arrangements, 10 after dropping enantiomers."""
    ligands = ['[F-:1]', '[Cl-:2]', '[Br-:3]', '[I-:4]', '[NH3:5]']
    c = ComplexFromLigands(ligands, '[Zn+2]', 'TBP')

    assert c.err_init is None
    assert len(c.GetStereomers(regime='CA', dropEnantiomers=False)) == 20
    assert len(c.GetStereomers(regime='CA', dropEnantiomers=True)) == 10
    assert c.AddConformer(maxAttempts=20) >= 0


def test_cpstar_mn_can_use_tetrahedral_starting_geometry():
    """Cp*(Mn)(NH3)(CO)2 should embed successfully with the new TET geometry."""
    ligands = [
        "[*:1]c1(C)c(C)c(C)c(C)c1C",
        "[NH3:2]",
        "[C-:3]#[O+]",
        "[C-:4]#[O+]",
    ]
    c = ComplexFromLigands(ligands, '[Mn]', 'TET')

    assert c.err_init is None
    assert len(c._DAs) == 4
    assert len(c._haptic_DAs) == 1

    flag = c.AddConformer(maxAttempts=20)
    assert flag >= 0

    conf = c.mol3D.GetConformer(flag)
    donor_idxs = [idx for idx, _ in sorted(c._DAs.items(), key=lambda item: item[1])]
    angles = [
        _angle_deg(conf, donor_idxs[i], c._idx_CA, donor_idxs[j])
        for i in range(len(donor_idxs) - 1)
        for j in range(i + 1, len(donor_idxs))
    ]
    for angle in angles:
        assert 90.0 < angle < 125.0, (
            f"Tetrahedral-start angle {angle:.1f} deg is outside the expected range"
        )


def test_ferrocene_sandwich_support():
    """Ferrocene should embed as a two-centroid sandwich with eta5 Cp rings."""
    ligands = [
        "[*:1][C-]1C=CC=C1",
        "[*:2][C-]1C=CC=C1",
    ]
    c = ComplexFromLigands(ligands, '[Fe+2]', 'SAN')

    assert c.err_init is None
    assert len(c._DAs) == 2
    assert len(c._haptic_DAs) == 2
    assert all(len(haptic_idxs) == 5 for haptic_idxs in c._haptic_DAs.values())
    assert len(c.GetStereomers(regime='CA', dropEnantiomers=True)) == 1

    flag = c.AddConformer(maxAttempts=20)
    assert flag >= 0
    xyz = c.ToXYZBlock(flag)
    assert _count_xyz_element(xyz, 'H') == 10

    conf = c.mol3D.GetConformer(flag)
    centroid_idxs = [idx for idx, _ in sorted(c._DAs.items(), key=lambda item: item[1])]
    angle = _angle_deg(conf, centroid_idxs[0], c._idx_CA, centroid_idxs[1])
    assert 170.0 < angle <= 180.0, (
        f"Centroid-Fe-centroid angle {angle:.1f} deg is outside the expected sandwich range"
    )
    expected = Complex._HapticDist[5]
    for idx in centroid_idxs:
        dist = np.linalg.norm(
            np.array(conf.GetAtomPosition(idx)) -
            np.array(conf.GetAtomPosition(c._idx_CA))
        )
        assert abs(dist - expected) < 0.5, (
            f"Fe-centroid distance {dist:.2f} A differs from HapticDist[5]={expected} A by more than 0.5 A"
        )
