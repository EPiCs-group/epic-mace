# W(PPh₃)₂(CO)₃(η²-H₂) Test Case

## Description
This test case generates the octahedral tungsten complex W(PPh₃)₂(CO)₃(η²-H₂), featuring a sigma-H₂ ligand coordinated in η² fashion.

## Input
- **Geometry**: Octahedral (OH)
- **Ligands**:
  - 2 × PPh₃ (triphenylphosphine)
  - 3 × CO (carbonyl)
  - 1 × η²-H₂ (sigma dihydrogen)
- **Central Atom**: W (tungsten)

## Expected Outcome
- **Stereoisomers**: 1 (no stereoisomer search due to regime: none)
- **Geometry**: Proper octahedral coordination with W-H₂ maintaining η² bonding
- **H₂ Coordination**: Side-on coordination with H-H distance ~0.85 Å and appropriate W-H distances
- **Files Generated**: `W_PPh3_2_CO_3_H2_iso0.xyz`

## Notes
This case validates that sigma-H₂ complexes maintain proper η² coordination geometry when appropriate constraints are applied to prevent UFF optimization from collapsing to hydride-like structures.