# Mn(PNP)(CO)₂(η²-H₂) Test Case

## Description
This test case generates the octahedral manganese complex Mn(PNP)(CO)₂(η²-H₂), where PNP is a tridentate phosphine-amine-phosphine ligand with iPr substituents.

## Input
- **Geometry**: Octahedral (OH)
- **Ligands**:
  - 1 × PNP (CC(C)[P:1](C(C)C)CC[N:1]CC[P:1](C(C)C)C(C)C - tridentate with P,N,P donors)
  - 2 × CO (carbonyl)
  - 1 × η²-H₂ (sigma dihydrogen)
- **Central Atom**: Mn (manganese)

## Expected Outcome
- **Stereoisomers**: 5 (found via stereoisomer search)
- **Geometry**: Proper octahedral coordination with meridional PNP arrangement
- **H₂ Coordination**: Side-on coordination with H-H distance ~0.85 Å and appropriate Mn-H distances
- **Files Generated**: `Mn_PNP_CO_2_H2_iso0.xyz` through `iso4.xyz`

## Notes
This case demonstrates stereoisomer generation for complexes with multi-dentate ligands and sigma-H₂. The PNP ligand occupies three coordination sites in meridional fashion, leaving three sites for CO and η²-H₂. The added constraints ensure proper η²-H₂ geometry is maintained.