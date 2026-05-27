"""Standalone runner: generate and save XYZ files for three Ru sigma/pi complexes.

Usage
-----
    python tests/run_sigma_pi_complexes.py [output_dir]

    output_dir defaults to  ./sigma_pi_output

Output layout
-------------
    sigma_pi_output/
    ├── RuSNS_deprot_sigmaH2_H_PPh3/
    │   ├── README.md
    │   └── RuSNS_deprot_sigmaH2_H_PPh3_iso0.xyz  (lowest-energy isomers)
    │   └── ...
    ├── RuNS_deprot_2sigmaH2_H_PPh3/
    │   └── ...
    └── RuNS_prot_Cp_H/
        └── ...

Each XYZ file holds the lowest-energy conformer of that isomer.
High-energy isomers (> E_REL_MAX kJ/mol above the minimum) are excluded.
"""

import sys
import datetime
from pathlib import Path

# Ensure the local development version of mace is used, not the pip-installed one.
# The haptic-ligand features tested here are present in the source repo but are
# not yet part of the published 0.5.0 package.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import mace


# ── settings ──────────────────────────────────────────────────────────────
NUM_CONFS    = 10      # conformers attempted per isomer
RMS_THRESH   = 0.5    # Å – RMSD filter during generation
E_REL_MAX    = 150.0  # kJ/mol – drop isomers with energy this far above minimum
                       # (UFF energies are large; use relative differences here)

COMPLEXES = [
    {
        # ── Complex 1 ───────────────────────────────────────────────────────
        # Ru(II) octahedral
        # SNS pincer, BOTH S coordinating, N deprotonated (amide, N⁻)
        # + one side-on σ-H₂ (Kubas, η²)
        # + one hydride (H⁻)
        # + PPh₃
        # Charge balance: Ru²⁺  +  N⁻ (−1)  +  H⁻ (−1) = neutral
        'name': 'RuSNS_deprot_sigmaH2_H_PPh3',
        'description': (
            'Ru(II) octahedral, SNS tridentate (N deprotonated), '
            'eta2-sigma-H2 (Kubas), hydride, PPh3'
        ),
        'ligands': [
            'CC[S:1]CC[N-:2]CC[S:3]CC',             # SNS tridentate, N⁻
            '[*:4]([H])[H]',                          # eta2-sigma-H2
            '[H-:5]',                                 # hydride
            '[P:6](c1ccccc1)(c1ccccc1)c1ccccc1',     # PPh3
        ],
        'CA'  : '[Ru+2]',
        'geom': 'OH',
    },
    {
        # ── Complex 2 ───────────────────────────────────────────────────────
        # Ru(II) octahedral
        # SNS pincer used as BIDENTATE (one S + N⁻); second S arm pendant
        # + TWO side-on σ-H₂ (η²)
        # + one hydride (H⁻)
        # + PPh₃
        # Charge balance: Ru²⁺  +  N⁻ (−1)  +  H⁻ (−1) = neutral
        'name': 'RuNS_deprot_2sigmaH2_H_PPh3',
        'description': (
            'Ru(II) octahedral, SNS as bidentate NS (N deprotonated, one S free), '
            '2x eta2-sigma-H2, hydride, PPh3'
        ),
        'ligands': [
            'CC[S:1]CC[N-:2]CCSCC',                  # NS bidentate, N⁻, second S free
            '[*:3]([H])[H]',                           # eta2-sigma-H2 #1
            '[*:4]([H])[H]',                           # eta2-sigma-H2 #2
            '[H-:5]',                                  # hydride
            '[P:6](c1ccccc1)(c1ccccc1)c1ccccc1',      # PPh3
        ],
        'CA'  : '[Ru+2]',
        'geom': 'OH',
    },
    {
        # ── Complex 3 ───────────────────────────────────────────────────────
        # Ru(II) half-sandwich (piano-stool)
        # SNS pincer used as BIDENTATE (one S + N-H); second S arm pendant
        # + η⁵-cyclopentadienyl (Cp⁻, pi-ligand)
        # + one hydride (H⁻)
        # Two phantom coordination positions are filled automatically by mace.
        # Charge balance: Ru²⁺  +  Cp⁻ (−1)  +  H⁻ (−1) = neutral
        'name': 'RuNS_prot_Cp_H',
        'description': (
            'Ru(II) half-sandwich, SNS as bidentate NS (N protonated, one S free), '
            'eta5-Cp, hydride'
        ),
        'ligands': [
            'CC[S:1]CC[NH:2]CCSCC',         # NS bidentate, N-H, second S free
            '[*:3][C-]1C=CC=C1',             # eta5-Cp (auto-expanded to full ring)
            '[H-:4]',                         # hydride
        ],
        'CA'  : '[Ru+2]',
        'geom': 'OH',
    },
]


# ── helpers ────────────────────────────────────────────────────────────────

def _run_complex(spec, out_root: Path) -> dict:
    name        = spec['name']
    description = spec['description']
    ligands     = spec['ligands']
    ca          = spec['CA']
    geom        = spec['geom']

    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"  {description}")
    print(f"{'='*60}")

    # ── assemble & find stereomers ─────────────────────────────────────────
    X0 = mace.ComplexFromLigands(ligands, ca, geom)

    print(f"  Donor atoms : {len(X0._DAs)}")
    print(f"  Haptic DAs  : {dict(X0._haptic_DAs)}")

    stereomers = X0.GetStereomers(regime='all', dropEnantiomers=True, merRule=True)
    print(f"  Stereomers  : {len(stereomers)}")

    # ── generate conformers ────────────────────────────────────────────────
    for X in stereomers:
        X.AddConformers(numConfs=NUM_CONFS, rmsThresh=RMS_THRESH)
        X.OrderConfsByEnergy()

    # ── energy filter ──────────────────────────────────────────────────────
    energies = [X.GetConfEnergy(0) if X.GetNumConformers() else float('inf')
                for X in stereomers]
    e_min = min(e for e in energies if e < float('inf'))
    kept = [X for X, e in zip(stereomers, energies) if e - e_min <= E_REL_MAX]
    dropped = len(stereomers) - len(kept)

    print(f"  Kept after ΔE ≤ {E_REL_MAX} kJ/mol filter : {len(kept)} "
          f"(dropped {dropped} high-energy isomers)")

    # ── save XYZ ───────────────────────────────────────────────────────────
    out_dir = out_root / name
    out_dir.mkdir(parents=True, exist_ok=True)

    saved_files = []
    for i, X in enumerate(kept):
        if X.GetNumConformers() == 0:
            print(f"    Isomer {i}: no conformers – skipped")
            continue
        fpath = out_dir / f'{name}_iso{i}.xyz'
        X.ToMultipleXYZ(str(fpath), confIds=[0])   # save only lowest-energy conf
        e = X.GetConfEnergy(0)
        de = e - e_min
        saved_files.append((str(fpath), e, de))
        print(f"    Isomer {i}: E = {e:.1f} kJ/mol  (ΔE = {de:.1f})  → {fpath.name}")

    # ── write README ───────────────────────────────────────────────────────
    _write_readme(out_dir, name, description, spec, stereomers, kept,
                  e_min, saved_files)

    return {
        'name': name,
        'n_stereomers': len(stereomers),
        'n_saved': len(saved_files),
        'e_min': e_min,
    }


def _write_readme(out_dir, name, description, spec, all_stereomers, kept,
                  e_min, saved_files):
    ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    lines = [
        f'# epic-mace: {name}',
        f'',
        f'**Generated:** {ts}  ',
        f'**epic-mace version:** {mace.__version__}',
        f'',
        f'## Description',
        f'',
        f'{description}',
        f'',
        f'## Input',
        f'',
        f'| Parameter | Value |',
        f'|-----------|-------|',
        f'| Geometry  | {spec["geom"]} (Octahedral) |',
        f'| Central atom | `{spec["CA"]}` |',
        f'',
        f'**Ligands:**',
        f'',
    ]
    for lig in spec['ligands']:
        lines.append(f'- `{lig}`')

    lines += [
        f'',
        f'## Stereomer generation',
        f'',
        f'| Parameter | Value |',
        f'|-----------|-------|',
        f'| Regime | `all` |',
        f'| Drop enantiomers | `True` |',
        f'| Mer rule | `True` |',
        f'| Conformers attempted | `{NUM_CONFS}` |',
        f'| RMSD threshold | `{RMS_THRESH}` Å |',
        f'| Energy window | `{E_REL_MAX}` kJ/mol |',
        f'',
        f'## Results',
        f'',
        f'Total stereomers found: **{len(all_stereomers)}**  ',
        f'Isomers saved (ΔE ≤ {E_REL_MAX} kJ/mol): **{len(saved_files)}**',
        f'',
        f'| File | E (kJ/mol) | ΔE (kJ/mol) |',
        f'|------|:----------:|:-----------:|',
    ]
    for fpath, e, de in saved_files:
        fname = Path(fpath).name
        lines.append(f'| `{fname}` | {e:.1f} | {de:.1f} |')

    p = out_dir / 'README.md'
    p.write_text('\n'.join(lines))
    print(f"    README → {p}")


# ── main ───────────────────────────────────────────────────────────────────

def main():
    out_root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('./sigma_pi_output')
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"\nepic-mace v{mace.__version__}")
    print(f"Output directory: {out_root.resolve()}")

    summary = []
    for spec in COMPLEXES:
        result = _run_complex(spec, out_root)
        summary.append(result)

    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}")
    for r in summary:
        print(f"  {r['name']}")
        print(f"    stereomers found : {r['n_stereomers']}")
        print(f"    isomers saved    : {r['n_saved']}")
        print(f"    lowest E         : {r['e_min']:.1f} kJ/mol")
    print(f"\nAll output written to: {out_root.resolve()}")


if __name__ == '__main__':
    main()
