# [PtCl₂(NH₃)(η²-C₂H₄)] — Zeise's-Salt Analogue (pi-complex test case)

## Description

This test case generates the square-planar Pt(II) complex [PtCl₂(NH₃)(η²-C₂H₄)],
a Zeise's-salt analogue in which one chlorido ligand of Zeise's anion [PtCl₃(C₂H₄)]⁻
is replaced by NH₃. It validates that epic-mace can correctly handle η²-alkene (pi)
ligands encoded via the centroid dummy-atom approach.

## Input

- **Geometry**: Square-planar (SP)
- **Central atom**: Pt²⁺
- **Ligands** (all donors carry `:1` — no coordination positions pre-assigned):
  - 2 × Cl⁻
  - 1 × NH₃
  - 1 × η²-C₂H₄ centroid `[*:1]([CH2])[CH2]`
- **Stereomer search**: `regime: CA` (central-atom permutations, enantiomers dropped)

All ligands share map number `:1` so that MACE does not pre-define the coordination
geometry and instead enumerates all distinct square-planar arrangements.

## Results

MACE found **2 stereoisomers**, as expected for a SP complex of type MA₂BC:

| File | UFF Energy (kcal/mol) | Geometry |
|------|----------------------|----------|
| `Pt_Cl2_NH3_C2H4_iso0.xyz` | 345.52 | **Isomer A**: C₂H₄ trans to Cl⁻, NH₃ trans to Cl⁻ (two Cl cis) |
| `Pt_Cl2_NH3_C2H4_iso1.xyz` | 349.70 | **Isomer B**: Cl⁻ trans to Cl⁻, C₂H₄ trans to NH₃ |

Each xyz file contains one representative conformer (lowest UFF energy from 10 attempts,
RMSD threshold 0.5 Å). The dummy atom `X` in the xyz output is the ethylene centroid.

## Geometry validation

The η²-alkene constraint logic produces proper pi-complex geometry:

- **C=C perpendicular to coordination plane**: angle between C=C vector and plane
  normal is ~8° (ideal = 0°)
- **Pt–centroid**: ~1.82 Å (target `HapticDist[2]` = 1.85 Å)
- **C–C distance**: ~1.58 Å (longer than free C=C due to centroid bonding; DFT
  refinement recommended)
- **Pt–C**: ~1.94 Å (both C atoms equidistant from Pt)

## Notes

- The η²-ethylene centroid occupies one SP coordination position and is detected as
  a haptic donor atom (hapticity 2, two C neighbours) by the haptic constraint logic.
- The force-field distinguishes σ-H₂ (both haptic atoms H) from η²-alkene (heavy atoms)
  using element-specific constraints: C=C distance [1.20–1.50 Å], M–C [2.0–2.5 Å],
  centroid–C = 0.67 Å, plus an equidistance constraint from a cis donor that pins the
  C=C perpendicular to the coordination plane.
- Real Pt–alkene centroid distances are ~2.0 Å — always follow with DFT optimisation.
- This test case is covered by `TestPtEthylene` in `tests/test_haptic_ligands.py`.
