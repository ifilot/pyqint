.. index:: basis-sets-and-molecules

Basis sets and molecules
========================

.. contents:: Table of Contents
    :depth: 3

Building molecules
------------------

Molecules can be efficiently built from the :code:`Molecule` class. For example,
to build the H\ :sub:`2` molecule, one can run the script below.

.. code-block:: python

    from pyqint import PyQInt, Molecule
    import numpy as np

    # construct integrator object
    integrator = PyQInt()

    # build hydrogen molecule
    mol = Molecule('H2')
    mol.add_atom('H', 0.0, 0.0, 0.0)
    mol.add_atom('H', 0.0, 0.0, 1.4)
    print(mol)

The output of the above script is::

    Molecule: H2
     H (0.000000,0.000000,0.000000)
     H (0.000000,0.000000,1.400000)


Using the MoleculeBuilder class
-------------------------------

Next to constructing molecules from scratch, one can also use the
:code:`MoleculeBuilder` class which contains a number of pre-generated molecules.

To see which molecules are available in the MoleculeBuilder class, run

.. code-block:: python

    from pyqint import MoleculeBuilder
    mol = MoleculeBuilder.print_list_molecules()

To load any of these molecules, one uses the :code:`from_name` function
as shown in the script below

.. code-block:: python

    from pyqint import MoleculeBuilder
    mol = MoleculeBuilder.from_name('ch4')
    print(mol)

The output of the above script shows the elements and the atom positions::

    Molecule: methane
     C     0.000000      0.000000      0.000000
     H     1.195756      1.195756      1.195756
     H    -1.195756     -1.195756      1.195756
     H    -1.195756      1.195756     -1.195756
     H     1.195756     -1.195756     -1.195756

.. note::
    Naming a molecule is completely optional and has no further implications
    on any of the calculations. To name a molecule, populate the :code:`name`
    member of the :code:`Molecule` class.

Alternatively, one can load molecules from a :code:`.xyz` file via the
:code:`from_file` routine.

.. code-block:: python

    mol = MoleculeBuilder.from_file('ch4.xyz')

.. warning::
    It is assumed that the positions inside the `.xyz` file are stored in
    **angstroms**. Internally, :program:`PyQInt` uses Bohr distances and the
    distances as reported in the :code:`.xyz` file are automatically converted.

Constructing basis functions for a molecule
-------------------------------------------

To construct the basis functions for a given molecule, one first needs to
construct the molecule after which the :code:`build_basis` function can be used
to construct a basis.

The following basis sets are supported. For each basis set, the range of atoms
that are supported are given:

* :code:`sto3g` (H-I)
* :code:`sto6g` (H-Kr)
* :code:`p321` (H-Cs)
* :code:`p631` (H-Zn)

The example code below builds the basis functions for the H\ :sub:`2` molecule:

.. code-block:: python

    from pyqint import PyQInt, Molecule
    import numpy as np

    # construct integrator object
    integrator = PyQInt()

    # build hydrogen molecule
    mol = Molecule('H2')
    mol.add_atom('H', 0.0, 0.0, 0.0)
    mol.add_atom('H', 0.0, 0.0, 1.4)
    cgfs, nuclei = mol.build_basis('sto3g')

    for cgf in cgfs:
        print(cgfs)

    for nucleus in nuclei:
        print(nucleus)

The output of the above script is::

    [CGF; R=(0.000000,0.000000,0.000000)
    01 | GTO : c=0.154329, alpha=3.425251, l=0, m=0, n=0, R=(0.000000,0.000000,0.000000)
    02 | GTO : c=0.535328, alpha=0.623914, l=0, m=0, n=0, R=(0.000000,0.000000,0.000000)
    03 | GTO : c=0.444635, alpha=0.168855, l=0, m=0, n=0, R=(0.000000,0.000000,0.000000)
    , CGF; R=(0.000000,0.000000,1.400000)
    01 | GTO : c=0.154329, alpha=3.425251, l=0, m=0, n=0, R=(0.000000,0.000000,1.400000)
    02 | GTO : c=0.535328, alpha=0.623914, l=0, m=0, n=0, R=(0.000000,0.000000,1.400000)
    03 | GTO : c=0.444635, alpha=0.168855, l=0, m=0, n=0, R=(0.000000,0.000000,1.400000)
    ]
    [CGF; R=(0.000000,0.000000,0.000000)
    01 | GTO : c=0.154329, alpha=3.425251, l=0, m=0, n=0, R=(0.000000,0.000000,0.000000)
    02 | GTO : c=0.535328, alpha=0.623914, l=0, m=0, n=0, R=(0.000000,0.000000,0.000000)
    03 | GTO : c=0.444635, alpha=0.168855, l=0, m=0, n=0, R=(0.000000,0.000000,0.000000)
    , CGF; R=(0.000000,0.000000,1.400000)
    01 | GTO : c=0.154329, alpha=3.425251, l=0, m=0, n=0, R=(0.000000,0.000000,1.400000)
    02 | GTO : c=0.535328, alpha=0.623914, l=0, m=0, n=0, R=(0.000000,0.000000,1.400000)
    03 | GTO : c=0.444635, alpha=0.168855, l=0, m=0, n=0, R=(0.000000,0.000000,1.400000)
    ]
    (array([0., 0., 0.]), 1)
    (array([0. , 0. , 1.4]), 1)

Parallel evaluation of integrals
--------------------------------

From a collection of Contracted Gaussian Functions, the complete set of overlap,
kinetic, nuclear attraction and two-electron integrals can be quickly evaluated
using the `build_integrals_openmp` function. The function will automatically determine
the number of available cores to allocate for this process.

.. code-block:: python

    from pyqint import PyQInt, Molecule
    import numpy as np

    # construct integrator object
    integrator = PyQInt()

    # build hydrogen molecule
    mol = Molecule()
    mol.add_atom('H', 0.0, 0.0, 0.0)
    mol.add_atom('H', 0.0, 0.0, 1.4)
    cgfs, nuclei = mol.build_basis('sto3g')

    # evaluate all integrals
    S, T, V, teint = integrator.build_integrals_openmp(cgfs, nuclei)

    print(S)
    print(T)
    print(V)
    print(teint)

The output of the above script is given by::

    [[1.00000011 0.6593185 ]
     [0.6593185  1.00000011]]
    [[0.76003161 0.23645446]
     [0.23645446 0.76003161]]
    [[-1.88044067 -1.19483464]
     [-1.19483464 -1.88044067]]
    [0.7746057639733748, 0.4441076656879813, 0.29702859983423036, 0.5696758530951017, 0.44410766568798105, 0.7746057639733748]