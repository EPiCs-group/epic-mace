Command line interface
======================

.. note::
    Before you start reading this section, we recommend to look through the package :ref:`cookbook <cookbook link>` to get a better understanding of the **epic-mace**'s features. This should not be a problem even if you have no programming experience, just concentrate on the input and output data.


There are two ways to specify the task for **epic-mace**: via an input file or via command line arguments. Since they can easily be transformed into each other, we will start by describing syntax of input files as a more convenient and reproducible method.

.. warning::
    Do not forget to activate the environment with the installed **epic-mace** package before using it!


Preparing input
---------------

To avoid the inconvenience of creating input files from scratch, one can use **epic-mace-quickstart** command, which creates templates of mace input file **mace_input.yaml** and substituents' file **substituents.yaml** in the working directory:

.. code-block:: bash

    >> epic-mace-quickstart
    >> ls
    mace_input.yaml  substituents.yaml

If you want to create these file in another directory, specify the corresponding path as the first argument:

.. code-block:: bash

    >> mkdir subdir
    >> epic-mace-quickstart subdir
    >> ls subdir
    mace_input.yaml  substituents.yaml

To launch **epic-mace** with the given input file, use the following command:

.. code-block:: bash
    
    >> epic-mace --input mace_input.yaml


Syntax of input file
--------------------

**epic-mace**'s input file is a `YAML`_-formatted text file, containing key-value pairs for the command line tool's arguments. They can be provided in any order, however, it is more convenient to divide them into five groups:
  
  1. Output
  2. Structure
  3. Stereomer search
  4. Conformer generation
  5. Filtering conformers
  6. Substituents


Output
^^^^^^

Path to an output directory:

.. code-block:: yaml
    
    # output directory
    out_dir: ./

- **out_dir** (default value: *./*): specifies the path to the output directory. If you use a relative path, do not forget that it is relative to the working directory from which the script is run, not to the directory where the input file is located.


Structure
^^^^^^^^^

Parameters describing structure of a complex:

.. code-block:: yaml
    
    # structure
    name: RhCl_MeCN_bipy
    geom: SP
    res-structs: 1
    # via complex
    complex: "[Cl-:1][Rh+]1([N:2]#CC)[N:3]2=CC([*])=CC=C2C2=[N:4]1C=C([*])C=C2 |$;;;;;;;;_R2;;;;;;;;_R1;;$,c:8,10,19,t:5,16,C:0.0,2.1,5.4,13.14|"

- **name**: a name of a system which will be used in names of output files.

- **geom**: molecular geometry of a metal complex and can take two values:
  
  - "OH" for octahedral geometry;
  - "SP" for square-planar geometry.

- **res-structs** (default value: *1*): a number of resonance structures that is using to compare complex stereomers. This parameter is required for symmetric polydentate ligands with asymmetric standard resonance structures only (see the "Issue of resonance structures" subsection in the "Complex initialization" section of :ref:`cookbook <cookbook link>`).

- **complex**: ChemAxon/RDKit SMILES of the complex to generate (see the :ref:`Input SMILES<input link>` section for the details).

The other way to define a complex is to provide SMILES of ligands and central atom:

.. code-block:: yaml

    # define structure via ligands & CA
    ligands:
      - "[*]C1=C[N:4]=C(C=C1)C1=[N:3]C=C([*])C=C1 |$_R1;;;;;;;;;;;_R2;;$,c:3,5,13,t:1,8,10|"
      - "[N:2]#CC"
      - "[Cl-:1]"
    CA: "[Rh+]"

- **ligands**: a list of ligands' SMILES.

- **CA**: SMILES of the central atom.

.. warning::
    System structure must be specified either via **complex** or via **ligands** and **CA**.

`Example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/01_complex>`_ of both approaches to specify a structure of a complex.


Stereomer search
^^^^^^^^^^^^^^^^

Parameters of a stereomer search:

.. code-block:: yaml

    # stereomer-search
    regime: all # all, CA, ligands, none
    get-enantiomers: false # true
    trans-cycle: no # if no, trans-position for DA-DA donor atoms not allowed
    mer-rule: true # false

- **regime** (default value: *all*): type of the stereomer search:
    
    - *all*: iterates over all stereocenters (`example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/02_stereomers_all>`_);
    - *ligands*: iterates over ligand's stereocenters only (`example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/03_stereomers_ligands>`_);
    - *CA*: changes configuration of central atom only (`example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/04_stereomers_CA>`_;
    - *none*: do not search for other stereomers (`example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/01_complex>`_).

- **get-enantiomers** (default value: *false*): if *true*, generates all possible enantiomers; otherwise leaves only one complex for an enantiomeric pair.

- **trans-cycle** (default value: *no*): minimal number of bonds between neighboring donor atoms required for the trans- spatial arrangement of the donor atoms. If *no*, such arrangement is considered as impossible (`example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/06_trans_pos>`_).

- **mer-rule** (default value: *true*): if *true*, applies empirical rule forbidding fac- configuration for the "rigid" DA-DA-DA fragments of the ligand (`example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/05_bad_mer_rule>`_).


Conformer generation
^^^^^^^^^^^^^^^^^^^^

Parameters of generation of 3D atomic coordinates:

.. code-block:: yaml

    # conformer-generation
    num-confs: 3
    rms-thresh: 1.0

- **num-confs** (default value: *10*): number of conformers to generate.

- **rms-thresh** (default value: *0.0*): drops one of two conformers if their RMSD is less than this threshold.


Filtering conformers
^^^^^^^^^^^^^^^^^^^^

Parameters of representative selection of low-energy conformers:

.. code-block:: yaml

    # conformer post-processing
    num-repr-confs: no # or positive integer
    e-rel-max: 25.0 # kJ/mol
    drop-close-energy: true # false

- **num-repr-confs** (default value: *no*): maximal number of representative conformers to return after filtering all generated conformers. If *no*, no post-processing are applied to the generated conformers.

- **e-rel-max** (default value: *25.0*): maximal relative energy of conformer not to be dropped from consideration.

- **drop-close-energy** (default value: *true*): if *true*, drops one of two conformers with difference in energy less than 0.1 kJ/mol.

`Example <https://github.com/EPiCs-group/epic-mace/tree/master/examples/08_Mn_CNP_chemilab>`_ of post-filtering of a large number of generated conformers.


Substituents
^^^^^^^^^^^^

Substituent info for modifying the core structure:

.. code-block:: yaml

    # substituents
    substituents-file: substituents.yaml # default
    R1: # name: SMILES must be defined in substituents file
      - H
      - NMe2
      - OMe
    R2:
      - H
      - CN
      - NO2

- **substituents-file** (default value: *./substituents.yaml*): path to a YAML-formatted file containing substituents' name-SMILES mapping. Each substituent must contain exactly one dummy atom designating the place of attachment. For the moment only monovalent substituents are supported. If an input system does not contain substituent, the substituents' file will be ignored.

  .. code-block:: yaml
      
      # Alk/Ar
      H: "[*][H]"
      Me: "[*]C"
      Ph: "[*]c1ccccc1"
      # *-oxy
      OH: "[*]O"
      OMe: "[*]OC"
      OAc: "[*]OC(=O)C"
      # amino
      NH2: "[*]N"
      NMe2: "[*]N(C)C"
      # halogens
      F: "[*]F"
      Cl: "[*]Cl"
      Br: "[*]Br"
      I: "[*]I"
      # acceptors
      CN: "[*]C#N"
      NO2: "[*]N(=O)=O"

- **R1**, **R2**, etc.: lists of substituent names specified in the substituents' file. Indices of indicated substituents must correspond to those in an input structure.

Some `examples <https://github.com/EPiCs-group/epic-mace/tree/cli/examples/07_subs>`_ of using substituents.


Command line arguments
----------------------

In addition to using an input file, the program can be run using the appropriate command line arguments. Thus, the input file discussed above would correspond to the following command:

.. code-block:: bash

    >> epic-mace ./ --name RhCl_MeCN_bipy --geom SP --ligands "[*]C1=C[N:4]=C(C=C1)C1=[N:3]C=C([*])C=C1 |$_R1;;;;;;;;;;;_R2;;$,c:3,5,13,t:1,8,10|" "[N:2]#CC" "[Cl-:1]" --CA "[Rh+]" --num-confs 3 --rms-thresh 1.0 --R1 H NMe2 OMe --R2 H CN NO2

We omitted **res-structs**, **substituents-file**, and all arguments of **stereomer search** and **filtering conformers** groups since they have default values. For more details see the help message:

.. code-block:: bash
    
    >> epic-mace -h


.. _YAML: https://yaml.org/


