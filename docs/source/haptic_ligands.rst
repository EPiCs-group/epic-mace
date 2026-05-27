.. _haptic link:

Haptic and sigma complexes
==========================

.. note::

    This feature is available from version 0.6.0 onwards.

**epic-mace** supports η²–η⁶ haptic ligands (pi-complexes) and σ-bound
two-electron donors such as σ-H₂ (Kubas-type) complexes. They are encoded
using a **centroid dummy atom** in SMILES, and they work within the existing
octahedral (``OH``), square-planar (``SP``), and sandwich (``SAN``) geometry
frameworks without any changes to stereomer generation or 3D embedding workflows.


Background
----------

Classical ligands in **epic-mace** donate through a single atom (N, O, P, C,
etc.). Haptic ligands such as η⁵-cyclopentadienyl or η²-ethylene coordinate
through a *face* or *bond* of the ligand, not a single atom. **epic-mace**
handles this via a phantom **centroid atom**: a dummy atom (``*``) with a map
number placed at the geometric centre of the coordinating fragment. The dummy
occupies one coordination position just like any other donor atom; the actual
π-system atoms are attached to it as ordinary bonds.

.. code-block:: text

    Metal ← [*:N](haptic-atom₁)(haptic-atom₂)...
             ↑ centroid (dummy, map number N)

The force field treats the centroid as the donor atom:

- the **M–centroid distance** is set from tabulated values in ``HapticDist``
  (keyed by hapticity η, i.e. the number of coordinating atoms);
- the **centroid–ring-atom distances** are constrained via ``HapticCentR``
  (geometric radius of the coordinating fragment), keeping the ligand correctly
  positioned during MM optimisation.

No hybridisation-based angle constraints are applied to centroid atoms, since
the ring or bond geometry is fully described by the centroid–atom distance
constraints.


SMILES encoding
---------------

The centroid atom is written as ``[*:N]`` where ``N`` is the coordination
position (1–6 for ``OH``, 1–4 for ``SP``). The dative bond goes **from the
centroid to the metal**, just like any other donor atom. The haptic atoms are
attached to the centroid with ordinary single bonds.

General pattern::

    [Metal]([other-ligands...])<-[*:N](haptic-atom₁)(haptic-atom₂)...

The hapticity η is usually determined automatically by counting the non-metal
neighbours of the centroid atom. For Cp/Cp*/arene-style ligands, a centroid
bonded to a single ring carbon is automatically expanded to the full 5- or
6-membered pi ring internally.


Examples
--------

σ-H₂ (Kubas complex, η²)
~~~~~~~~~~~~~~~~~~~~~~~~~

A W(0) Kubas-type complex with five CO/PR₃ ligands and one coordinated H₂
molecule. The H₂ centroid occupies position 6:

.. code-block:: python

    from mace import Complex

    smiles = "[W]([CO:1])([CO:2])([CO:3])([PH3:4])([PH3:5])<-[*:6]([H])[H]"
    c = Complex(smiles, "OH")
    c.AddConformer()

The two H atoms are attached to the centroid; the force field pins the
centroid 1.85 Å from W and each H atom 0.55 Å from the centroid (≈ half the
free H–H bond length).

η²-alkene (Zeise's salt type, SP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Square-planar Pt(II) complex with one η²-ethylene ligand:

.. code-block:: python

    smiles = "[Pt]([Cl:1])([Cl:2])([Cl:3])<-[*:4]([CH2])[CH2]"
    c = Complex(smiles, "SP")
    c.AddConformer()

η³-allyl
~~~~~~~~~

.. code-block:: python

    smiles = "[Pd]([Cl:1])([Cl:2])<-[*:3]([CH2])[CH][CH2]"
    c = Complex(smiles, "SP")

η⁵-cyclopentadienyl (half-sandwich, OH)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An octahedral Rh complex with a Cp ring occupying one coordination position
and five classical ligands in the remaining positions:

.. code-block:: python

    smiles = "[Rh]([CO:1])([CO:2])([CO:3])([Cl:4])([Cl:5])<-[*:6]1[cH][cH][cH][cH][cH]1"
    c = Complex(smiles, "OH")
    c.AddConformer()

η⁶-arene (SP)
~~~~~~~~~~~~~~

.. code-block:: python

    smiles = "[Ru]([Cl:1])([Cl:2])([Cl:3])<-[*:4]1[cH][cH][cH][cH][cH][cH]1"
    c = Complex(smiles, "SP")
    c.AddConformer()


Stereomer generation
--------------------

Stereomer search works exactly as for classical complexes. The centroid dummy
atom is a regular donor atom with a map number, so all permutation and
equivalence logic is unchanged. For example:

.. code-block:: python

    smiles = "[Rh]([CO:1])([CO:2])([CO:3])<-[*:4][cH][cH][cH][cH][cH]"
    c = Complex(smiles, "OH")
    stereomers = c.GetStereomers(regime="CA", dropEnantiomers=True)


Default distance parameters
-----------------------------

The tables below show the default distances used by the force field. They are
geometry- and metal-independent approximations; real M–centroid distances vary
with metal identity and oxidation state.

**M–centroid distance** (``HapticDist``, Å):

+-------------+-------+-----------------------------+
| Hapticity η | Value | Typical systems             |
+=============+=======+=============================+
| η²          | 1.85  | σ-H₂, alkene, alkyne        |
+-------------+-------+-----------------------------+
| η³          | 2.00  | allyl, open-ring Cp         |
+-------------+-------+-----------------------------+
| η⁴          | 2.00  | butadiene, COD fragment     |
+-------------+-------+-----------------------------+
| η⁵          | 1.90  | Cp, indenyl                 |
+-------------+-------+-----------------------------+
| η⁶          | 1.75  | benzene, arene              |
+-------------+-------+-----------------------------+

**Centroid–haptic-atom distance** (``HapticCentR``, Å):

+-------------+-------+----------------------------------------------+
| Hapticity η | Value | Geometric basis                              |
+=============+=======+==============================================+
| η²          | 0.55  | ≈ half H–H or half C=C bond                  |
+-------------+-------+----------------------------------------------+
| η³          | 1.20  | allyl half-width                             |
+-------------+-------+----------------------------------------------+
| η⁴          | 1.15  | diene half-width                             |
+-------------+-------+----------------------------------------------+
| η⁵          | 1.21  | Cp ring radius (R = 0.851 × C–C, C–C = 1.40)|
+-------------+-------+----------------------------------------------+
| η⁶          | 1.40  | arene ring radius (R = C–C = 1.40)           |
+-------------+-------+----------------------------------------------+

To override a default for a specific complex, modify the class attribute
before creating the object:

.. code-block:: python

    from mace import Complex

    # Tighten M-centroid distance for a W(0) σ-H₂ complex
    Complex._HapticDist[2] = 1.78
    c = Complex("[W]([CO:1])([CO:2])([CO:3])([PH3:4])([PH3:5])<-[*:6]([H])[H]", "OH")


Known limitations
-----------------

- **Distance defaults are metal-independent.** For quantitatively accurate
  starting geometries, always follow **epic-mace** output with a DFT
  optimisation. The MM structures are intended as reasonable starting points
  only.

- **Full sandwiches (ferrocene-type)** are supported through the ``SAN``
  geometry, with two opposed haptic centroids around the metal.

- **Piano-stool [(η⁵-Cp)ML₃] complexes** are handled with the ``OH`` geometry
  by placing the Cp centroid at one position and the three L ligands at three
  of the remaining five positions (with two dummy-filled positions). A dedicated
  ``HS`` (half-sandwich) geometry type is planned.

- **Substituents on haptic ligands** (e.g. methylated Cp) are fully supported —
  just attach them to the ring carbon atoms in the normal way.
