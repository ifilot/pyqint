.. _user-interface:
.. index:: userinterface

User Interface
##############

.. contents:: Table of Contents
    :depth: 3

Basis functions
===============

Gaussian type orbitals (GTO)
----------------------------

:program:`PyQInt` uses cartesian Gaussian type orbitals as given by

.. math::

    \Phi(\alpha,l,m,n,\vec{R}) = N (x - X)^{l} (y - Y)^{m} (z - Z)^{n} \exp \left(- \alpha |\vec{r} - \vec{R}|^{2} \right)

wherein :math:`\alpha` is the exponent, :math:`\vec{R} = \left(X,Y,Z\right)` the
position of the orbital, :math:`(l,m,n)` the orders of the pre-exponential
polynomial, and :math:`N` a normalization constant such that

.. math::

    \left< \Phi | \Phi \right> = 1

.. note::
    The normalization constant is automatically calculated by `PyQInt` and does not have
    to be supplied by the user.

GTOs are a fundamental building block of CGF (see below) and typically a user would
not directly work with them. Nevertheless, GTO objects can be constructed as follows::

    from pyqint import PyQInt, cgf, gto

    alpha = 0.5
    l,m,n = 0,0,0
    p = (0,0,0)
    G = gto(1.0, p, alpha, l, m, n)

.. note::
    If you work with individual GTOs, the first parameter to construct the GTO
    should have a value of 1.0. This first parameter corresponds to the linear
    expansion coefficient used in the formulation of Contracted Gaussian Functions
    (see below).

Contracted Gaussian Functions (CGF)
-----------------------------------

Several GTOs can be combined to produce a so-called Contracted Gaussian Functional which
is esentially a linear combination of GTOs as given by

.. math::

    \phi = \sum_{i} c_{i} \Phi_{i}(\alpha,l,m,n,\vec{R})

To build a CGF, we first have to produce the CGF object and then
add GTOs to it::

    from pyqint import PyQInt, cgf

    # build cgf for hydrogen separated by 1.4 a.u.
    cgf1 = cgf([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

Integral evaluation
===================

Electronic structure calculations require the construction of molecular
integrals. Here, an overview is given of the integrals involved and how these
can be evaluated using :program:`PyQInt`.

Overlap integrals
-----------------

Overlap integrals effectively probe the overlap between two CGFs and are given by

.. math::

    S_{ij} = \left< \phi_{i} | \phi_{j} \right>

CGFs should be normalized and as such, their self-overlap should be equal to
1. In the code snippet below, the overlap matrix :math:`\mathbf{S}` is
calculated for a basis set composed of the two :math:`1s` atomic orbitals on H which
are separated by a distance of 1.4 Bohr.

.. code-block:: python

    from pyqint import PyQInt, cgf
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build cgf for hydrogen separated by 1.4 a.u.
    cgf1 = cgf([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF
    cgf2 = deepcopy(cgf1)
    cgf2.p[2] = 1.4

    # construct empty matrix
    S = np.zeros((2,2))
    S[0,0] = integrator.overlap(cgf1, cgf1)
    S[0,1] = S[1,0] = integrator.overlap(cgf1, cgf2)
    S[1,1] = integrator.overlap(cgf2, cgf2)

    # output result
    print(S)

The result of this script is::

    [[1.00000011 0.6593185 ]
     [0.6593185  1.00000011]]

Kinetic integrals
-----------------

Kinetic integrals determine the kinetic energy of a given orbital and are given
by

.. math::

    T_{ij} = \left< \phi_{i} \left| -\frac{1}{2} \nabla^{2} \right| \phi_{j} \right>

In the code snippet below, the kinetic energy matrix :math:`\mathbf{T}` is
calculated for a basis set composed of the two :math:`1s` atomic orbitals on H which
are separated by a distance of 1.4 Bohr.

.. code-block:: python

    from pyqint import PyQInt, cgf, gto
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build cgf for hydrogen separated by 1.4 a.u.
    cgf1 = cgf([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF
    cgf2 = deepcopy(cgf1)
    cgf2.p[2] = 1.4

    # construct empty matrix
    T = np.zeros((2,2))
    T[0,0] = integrator.kinetic(cgf1, cgf1)
    T[0,1] = T[1,0] = integrator.kinetic(cgf1, cgf2)
    T[1,1] = integrator.kinetic(cgf2, cgf2)

    # output result
    print(T)

The result of the above script is::

    [[0.76003161 0.23645446]
     [0.23645446 0.76003161]]

Nuclear attraction integrals
----------------------------

Nuclear attraction integrals determine the attraction between a given nucleus
and the atomic orbital and are given by

.. math::

    V_{ij} = \left< \phi_{i} \left| -\frac{Z_{c}}{r_{i,c}} \right| \phi_{j} \right>

In the code snippet below, the nuclear attraction energy matrices :math:`\mathbf{V}_{1}`
and :math:`\mathbf{V}_{2}` are calculated for a basis set composed of the
two :math:`1s` atomic orbitals on H which are separated by a distance of 1.4 Bohr.
Due to the symmetry of the system, the nuclear attraction matrices for each of
the nuclei are the same.

.. code-block:: python

    from pyqint import PyQInt, cgf, gto
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build cgf for hydrogen separated by 1.4 a.u.
    cgf1 = cgf([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF
    cgf2 = deepcopy(cgf1)
    cgf2.p[2] = 1.4

    # Build nuclear attraction integrals
    V1 = np.zeros((2,2))
    V1[0,0] = integrator.nuclear(cgf1, cgf1, cgf1.p, 1)
    V1[0,1] = V1[1,0] = integrator.nuclear(cgf1, cgf2, cgf1.p, 1)
    V1[1,1] = integrator.nuclear(cgf2, cgf2, cgf1.p, 1)

    V2 = np.zeros((2,2))
    V2[0,0] = integrator.nuclear(cgf1, cgf1, cgf2.p, 1)
    V2[0,1] = V2[1,0] = integrator.nuclear(cgf1, cgf2, cgf2.p, 1)
    V2[1,1] = integrator.nuclear(cgf2, cgf2, cgf2.p, 1)

    # print result
    print(V1)
    print(V2)

The result of the above script is::

    [[-1.22661358 -0.59741732]
     [-0.59741732 -0.6538271 ]]
    [[-0.6538271  -0.59741732]
     [-0.59741732 -1.22661358]]

Two-electron integrals
----------------------

Two electron integrals capture electron-electron interactions, specifically
electron-electron repulsion and electron exchange. They are defined as

.. math::

    (i,j,k,l) = \left< \phi_{i}(x_{1})\phi_{j}(x_{2}) \left| r_{12}^{-1} \right| \phi_{k}(x_{1})\phi_{l}(x_{2}) \right>

The two-electron integrals are the most expensive terms to calculate in any
electronic structure calculation due to their :math:`N^{4}` scaling where
:math:`N` is the number of basis functions.

.. note::
    :program:`PyQInt` offers a `separate routine <#parallel-evaluation-of-integrals>`_
    for the efficient evaluation of all the integrals including the two electron integrals.

Although there are essentially :math:`N^{4}` different two-electron integrals,
due to certain symmetries the number of unique two-electron integrals is smaller.
In the script below, the six unique two-electron integrals for the H\ :sub:`2`
system are calculated.

.. code-block:: python

    from pyqint import PyQInt, cgf, gto
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build cgf for hydrogen separated by 1.4 a.u.
    cgf1 = cgf([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF
    cgf2 = deepcopy(cgf1)
    cgf2.p[2] = 1.4

    T1111 = integrator.repulsion(cgf1, cgf1, cgf1, cgf1)
    T1122 = integrator.repulsion(cgf1, cgf1, cgf2, cgf2)
    T1112 = integrator.repulsion(cgf1, cgf1, cgf1, cgf2)
    T2121 = integrator.repulsion(cgf2, cgf1, cgf2, cgf1)
    T1222 = integrator.repulsion(cgf1, cgf2, cgf2, cgf2)
    T2211 = integrator.repulsion(cgf2, cgf2, cgf1, cgf1)

    print(T1111)
    print(T1122)
    print(T1112)
    print(T2121)
    print(T1222)
    print(T2211)

The output of the above script is given by::

    0.7746057639733748
    0.5696758530951017
    0.44410766568798127
    0.29702859983423036
    0.4441076656879813
    0.5696758530951017

Dipole-moment integrals
-----------------------

Dipole-moment integrals are defined as

.. math::

    \mu_{x,i,j} = \left< \phi_{i}(x_{1}) \left| x \right| \phi_{j}(x_{1}) \right>

.. math::
    \mu_{y,i,j} = \left< \phi_{i}(x_{1}) \left| y \right| \phi_{j}(x_{1}) \right>

.. math::
    \mu_{z,i,j} = \left< \phi_{i}(x_{1}) \left| z \right| \phi_{j}(x_{1}) \right>

and are evaluated with respect to the coordinate center of the system. Dipole moments
are vector quantities, but in this implementation the dipoles are evaluated
in the :math:`x`, :math:`y`, :math:`z` separately.

In the script below, the dipole integrals are evaluated for the H\ :sub:`2`\ O
molecule using a :code:`sto3g` basis set and in each cartesian direction. The result
is collected in a three-dimensional array.

.. code-block:: python

    from pyqint import PyQInt, Molecule
    import numpy as np

    # construct integrator object
    integrator = PyQInt()

    # build water molecule
    mol = Molecule("H2O")
    mol.add_atom('O',  0.00000, -0.07579, 0.0000, unit='angstrom')
    mol.add_atom('H',  0.86681,  0.60144, 0.0000, unit='angstrom')
    mol.add_atom('H', -0.86681,  0.60144, 0.0000, unit='angstrom')
    cgfs, nuclei = mol.build_basis('sto3g')

    N = len(cgfs)
    D = np.zeros((N,N,3))
    for i in range(N):
        for j in range(i,N):
            for k in range(0,3): # loop over directions
                D[i,j,k] = integrator.dipole(cgfs[i], cgfs[j], k)

    print(D)

The result of the above script is::

    [[[ 0.00000000e+00 -1.43222417e-01  0.00000000e+00]
      [ 0.00000000e+00 -3.39013356e-02  0.00000000e+00]
      [ 5.07919476e-02  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  5.07919476e-02  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  5.07919476e-02]
      [ 2.22964944e-03 -3.75854187e-03  0.00000000e+00]
      [-2.22964944e-03 -3.75854187e-03  0.00000000e+00]]

     [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00 -1.43222278e-01  0.00000000e+00]
      [ 6.41172506e-01  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  6.41172506e-01  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  6.41172506e-01]
      [ 2.62741706e-01  1.49973767e-01  0.00000000e+00]
      [-2.62741706e-01  1.49973767e-01  0.00000000e+00]]

     [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00 -1.43222278e-01  0.00000000e+00]
      [-9.08620418e-18  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 4.37629746e-01  1.08953250e-01  0.00000000e+00]
      [ 4.37629746e-01 -1.08953250e-01  0.00000000e+00]]

     [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00 -1.43222278e-01  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00 -9.08620418e-18]
      [ 1.47399486e-01  3.34092154e-01  0.00000000e+00]
      [-1.47399486e-01  3.34092154e-01  0.00000000e+00]]

     [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00 -1.43222278e-01  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  2.48968067e-01]
      [ 0.00000000e+00  0.00000000e+00  2.48968067e-01]]

     [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 1.63803356e+00  1.13655692e+00  0.00000000e+00]
      [-1.38777878e-17  2.06582174e-01  0.00000000e+00]]

     [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
      [-1.63803356e+00  1.13655692e+00  0.00000000e+00]]]

.. note::
    Each row in the above output corresponds to the dipole moment **vector**.
    There are in total 7 blocks to be observed and each block contains 7
    rows. Each block corresponds to a different basis function in the *bra*
    and each row inside a block loops over the different basis functions in the
    *ket*.

Basis sets and molecules
========================

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

The following molecules are available:

* benzene
* bf3
* ch4
* co
* co2
* ethylene
* h2
* h2o
* he
* lih
* nh3

To load any of these molecules, one uses the :code:`from_name` function
as shown in the script below

.. code-block:: python

    from pyqint import MoleculeBuilder

    mol = MoleculeBuilder().from_name('ch4')
    mol.name = 'CH4'

    print(mol)

The output of the above script shows the elements and the atom positions::

    Molecule: CH4
     C (0.000000,0.000000,0.000000)
     H (1.195756,1.195756,1.195756)
     H (-1.195756,-1.195756,1.195756)
     H (-1.195756,1.195756,-1.195756)
     H (1.195756,-1.195756,-1.195756)

.. note::
    Naming a molecule is completely optional and has no further implications
    on any of the calculations. To name a molecule, populate the :code:`name`
    member of the :code:`Molecule` class.

Alternatively, one can load molecules from a :code:`.xyz` file via the
:code:`from_file` routine.

.. code-block:: python

    mol = MoleculeBuilder().from_file('ch4.xyz')

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

    [<pyqint.cgf.cgf object at 0x000001BDEDB37430>, <pyqint.cgf.cgf object at 0x000001BDEDB37F10>]
    [<pyqint.cgf.cgf object at 0x000001BDEDB37430>, <pyqint.cgf.cgf object at 0x000001BDEDB37F10>]
    [array([0., 0., 0.]), 1]
    [array([0. , 0. , 1.4]), 1]

Parallel evaluation of integrals
--------------------------------

From a collection of Contracted Gaussian Functions, the complete set of overlap,
kinetic, nuclear attraction and two-electron integrals can be quickly evaluated
using the `build_integrals` function. Using the `npar` argument, the number of
threads to be spawned can be set.

.. code-block:: python

    from pyqint import PyQInt, Molecule
    import numpy as np
    import multiprocessing

    # construct integrator object
    integrator = PyQInt()

    # build hydrogen molecule
    mol = Molecule()
    mol.add_atom('H', 0.0, 0.0, 0.0)
    mol.add_atom('H', 0.0, 0.0, 1.4)
    cgfs, nuclei = mol.build_basis('sto3g')

    # evaluate all integrals
    ncpu = multiprocessing.cpu_count()
    S, T, V, teint = integrator.build_integrals(cgfs, nuclei, npar=ncpu, verbose=False)

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

Electronic structure calculations
=================================

Hartree-Fock
------------

The Hartree-Fock procedure is readily available as a separate class in the
:program:`PyQInt` package. It gives rich output allowing the user to investigate
the Hartree-Fock coefficient optimization procedure in detail.

.. code-block:: python

    from pyqint import PyQInt, Molecule, HF
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    def main():
        # calculate sto3g coefficients for h2o
        cgfs, coeff = calculate_co()

        # visualize orbitals
        fig, ax = plt.subplots(2,3, figsize=(18,10))
        for i in range(0,2):
            for j in range(0,3):
                dens = plot_wavefunction(cgfs, coeff[:,i*3+j])
                limit = max(abs(np.min(dens)), abs(np.max(dens)) )
                im = ax[i,j].imshow(dens, origin='lower', interpolation='bilinear',
                  extent=[-2,2,-2,2], cmap='PiYG', vmin=-limit, vmax=limit)
                ax[i,j].set_xlabel('Distance a.u.')
                ax[i,j].set_ylabel('Distance a.u.')
                divider = make_axes_locatable(ax[i,j])
                cax = divider.append_axes('right', size='5%', pad=0.05)
                fig.colorbar(im, cax=cax, orientation='vertical')

    def calculate_co():
        mol = Molecule()
        mol.add_atom('C', 0.0, -0.5, 0.0)
        mol.add_atom('O', 0.0, 0.5, 0.0)

        result = HF().rhf(mol, 'sto3g')

        return result['cgfs'], result['orbc']

    def plot_wavefunction(cgfs, coeff):
        # build integrator
        integrator = PyQInt()

        # build grid
        x = np.linspace(-2, 2, 100)
        y = np.linspace(-2, 2, 100)
        xx, yy = np.meshgrid(x,y)
        zz = np.zeros(len(x) * len(y))
        grid = np.vstack([xx.flatten(), yy.flatten(), zz]).reshape(3,-1).T
        res = integrator.plot_wavefunction(grid, coeff, cgfs).reshape((len(y), len(x)))

        return res

    if __name__ == '__main__':
        main()

.. figure:: _static/img/co.jpg

    Canonical molecular orbitals of CO visualized using contour plots.

Result dictionary
-----------------

The result of a Hartree-Fock calculation is captured inside a dictionary
object. This dictionary objects contains the following keys

.. list-table:: Description of the data contained in the result library
   :widths: 25 75
   :header-rows: 1

   * - Key
     - Description
   * - :code:`energy`
     - Final energy of the electronic structure calculation
   * - :code:`nuclei`
     - List of elements and their position in Bohr units
   * - :code:`cgfs`
     - List of contracted Gaussian functional objects
   * - :code:`energies`
     - List of energies during the self-convergence procedure
   * - :code:`orbe`
     - Orbital energies (converged) (array of N element)
   * - :code:`orbc`
     - Orbital coefficients (converted) (matrix of N x N elements)
   * - :code:`density`
     - Density matrix :math:`\mathbf{P}`
   * - :code:`fock`
     - Fock matrix :math:`\mathbf{F}`
   * - :code:`transform`
     - Unitary transformation matrix :math:`\mathbf{X}`
   * - :code:`overlap`
     - Overlap matrix :math:`\mathbf{S}`
   * - :code:`kinetic`
     - Kinetic energy matrix :math:`\mathbf{T}`
   * - :code:`nuclear`
     - Nuclear attraction matrix :math:`\mathbf{V}`
   * - :code:`hcore`
     - Core Hamiltonian matrix :math:`\mathbf{H_\textrm{core}}`
   * - :code:`tetensor`
     - Two-electron tensor object :math:`(i,j,k,l)`
   * - :code:`time_stats`
     - Time statistics object
   * - :code:`ecore`
     - Sum of kinetic and nuclear attraction energy
   * - :code:`ekin`
     - Total kinetic energy
   * - :code:`enuc`
     - Total nuclear attraction energy
   * - :code:`erep`
     - Total electron-electron repulsion energy
   * - :code:`ex`
     - Total exchange energy
   * - :code:`enucrep`
     - Electrostatic repulsion energy of the nuclei
   * - :code:`nelec`
     - Total number of electrons
   * - :code:`forces`
     - Forces on the atoms (if calculated, else :code:`None`)

To provide an example how one can use the above data, let us consider the
situation wherein the user wants to decompose the individual components of the
total energy as given by

.. math::

    E_{\textrm{total}} = E_{\textrm{kin}} + E_{\textrm{nuc}} + E_{\textrm{e-e}} + E_{\textrm{ex}} + E_{\textrm{nuc,rep}}

Via the script below, one can easily verify that the above equation holds and
that the total energy is indeed the sum of the kinetic, nuclear attraction,
electron-electron repulsion, exchange and nuclear repulsion energies within a
Hartree-Fock calculation.

.. code-block:: python

    from pyqint import MoleculeBuilder,HF

    mol = MoleculeBuilder().from_name('ch4')
    mol.name = 'CH4'

    res = HF().rhf(mol, 'sto3g')
    print()
    print('Kinetic energy: ', res['ekin'])
    print('Nuclear attraction energy: ', res['enuc'])
    print('Electron-electron repulsion: ', res['erep'])
    print('Exchange energy: ', res['ex'])
    print('Repulsion between nuclei: ', res['enucrep'])
    print()
    print('Total energy: ', res['energy'])
    print('Sum of the individual terms: ',
          res['ekin'] + res['enuc'] + res['erep'] + res['ex'] + res['enucrep'])

The output of the above script yields::

    Kinetic energy:  39.42613774982387
    Nuclear attraction energy:  -118.63789179775034
    Electron-electron repulsion:  32.7324270326041
    Exchange energy:  -6.609004673631048
    Repulsion between nuclei:  13.362026647057352

    Total energy:  -39.72630504189621
    Sum of the individual terms:  -39.726305041896055

Orbital visualization
=====================

Since orbitals are essentially three-dimensional scalar fields, there are two
useful procedures to visualize them. The scalar field can either be projected
onto a plane, creating so-called contour plots. Alternatively, a specific
value (i.e. the isovalue) of the scalar field can be chosen and all points in
space that have this value can be tied together creating a so-called isosurface.

Contour plots can be easily created using `matplotlib <https://matplotlib.org/>`_.
For the creation of isosurfaces, we use `PyTessel <https://pytessel.imc-tue.nl.>`_.

Contour plots
-------------

.. code-block:: python

    from pyqint import PyQInt, Molecule
    import matplotlib.pyplot as plt
    import numpy as np

    # coefficients (calculated by Hartree-Fock using a sto3g basis set)
    coeff = [8.37612e-17, -2.73592e-16,  -0.713011, -1.8627e-17, 9.53496e-17, -0.379323,  0.379323]

    # construct integrator object
    integrator = PyQInt()

    # build water molecule
    mol = Molecule('H2O')
    mol.add_atom('O', 0.0, 0.0, 0.0)
    mol.add_atom('H', 0.7570, 0.5860, 0.0)
    mol.add_atom('H', -0.7570, 0.5860, 0.0)
    cgfs, nuclei = mol.build_basis('sto3g')

    # build grid
    x = np.linspace(-2, 2, 50)
    y = np.linspace(-2, 2, 50)
    xx, yy = np.meshgrid(x,y)
    zz = np.zeros(len(x) * len(y))
    grid = np.vstack([xx.flatten(), yy.flatten(), zz]).reshape(3,-1).T
    res = integrator.plot_wavefunction(grid, coeff, cgfs).reshape((len(y), len(x)))

    # plot wave function
    plt.imshow(res, origin='lower', extent=[-2,2,-2,2], cmap='PiYG')
    plt.colorbar()
    plt.title('1b$_{2}$ Molecular orbital of H$_{2}$O')


Constructing isosurfaces
------------------------

.. note::
    Isosurface generation requires the :program:`PyTessel` package to be
    installed. More information can be found `here <https://pytessel.imc-tue.nl>`_.

.. code-block:: python

    from pyqint import PyQInt, Molecule, HF
    import numpy as np
    from pytessel import PyTessel

    def main():
        # calculate sto3g coefficients for h2o
        cgfs, coeff = calculate_co()

        # build isosurface of the fifth MO
        # isovalue = 0.1
        # store result as .ply file
        build_isosurface('co_04.ply', cgfs, coeff[:,4], 0.1)

    def build_isosurface(filename, cgfs, coeff, isovalue):
        # generate some data
        sz = 100
        integrator = PyQInt()
        grid = integrator.build_rectgrid3d(-5, 5, sz)
        scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff, cgfs), (sz, sz, sz))
        unitcell = np.diag(np.ones(3) * 10.0)

        pytessel = PyTessel()
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue)
        pytessel.write_ply(filename, vertices, normals, indices)

    def calculate_co():
        mol = Molecule()
        mol.add_atom('C', 0.0, -0.5, 0.0)
        mol.add_atom('O', 0.0, 0.5, 0.0)

        result = HF().rhf(mol, 'sto3g')

        return result['cgfs'], result['orbc']

    if __name__ == '__main__':
        main()

Orbital localization: Foster-Boys
=================================

The code below first performs a Hartree-Fock calculation on the CO molecule
after which the localized molecular orbitals are calculated using the
`Foster-Boys method <https://en.wikipedia.org/wiki/Localized_molecular_orbitals#Foster-Boys>`_.
The Foster-Boys localization procedure is present as a separate class in the
:program:`PyQInt` package. It takes the output of a Hartree-Fock calculation
as its input.

.. note::
    The code below uses the PyTessel package for constructing the isosurfaces.
    PyTessel is an external package for easy construction of isosurfaces from
    scalar fields. More information is given `in the corresponding section below <#constructing-isosurfaces>`_.

.. code-block:: python

    from pyqint import Molecule, HF, PyQInt, FosterBoys
    import pyqint
    import numpy as np
    from pytessel import PyTessel

    def main():
        res = calculate_co(1.145414)
        resfb = FosterBoys(res).run()

        for i in range(len(res['cgfs'])):
            build_isosurface('MO_%03i' % (i+1),
                             res['cgfs'],
                             resfb['orbc'][:,i],
                             0.1)

    def calculate_co(d):
        """
        Full function for evaluation
        """
        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, -d/2, unit='angstrom')
        mol.add_atom('O', 0.0, 0.0,  d/2, unit='angstrom')

        result = HF().rhf(mol, 'sto3g')

        return result

    def build_isosurface(filename, cgfs, coeff, isovalue, sz=5, npts=100):
        # generate some data
        isovalue = np.abs(isovalue)
        integrator = PyQInt()
        grid = integrator.build_rectgrid3d(-sz, sz, npts)
        scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff, cgfs), (npts, npts, npts))
        unitcell = np.diag(np.ones(3) * 2 * sz)

        pytessel = PyTessel()
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue)
        fname = filename + '_pos.ply'
        pytessel.write_ply(fname, vertices, normals, indices)

        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), -isovalue)
        fname = filename + '_neg.ply'
        pytessel.write_ply(fname, vertices, normals, indices)

    if __name__ == '__main__':
        main()

.. figure:: _static/img/co_canonical_isosurfaces.jpg

    Canonical molecular orbitals of CO visualized using isosurfaces with an
    isovalue of +/-0.03.

.. figure:: _static/img/co_fosterboys_isosurfaces.jpg

    Localized molecular orbitals of CO visualized using isosurfaces with an
    isovalue of +/-0.03. Note that the localization procedure has only been
    applied to the occupied molecular orbitals. Observe that the localized
    orbitals contain a triple-degenerate state corresponding to the triple
    bond and two lone pairs for C and O.

Geometry optimization
=====================

Crystal Orbital Hamilton Population Analysis
============================================
