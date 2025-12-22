.. index:: integral-evaluation

Integral evaluation
###################

.. contents:: Table of Contents
    :depth: 3

Electronic structure calculations require the construction of molecular
integrals. Here, an overview is given of the integrals involved and how these
can be evaluated using :program:`PyQInt`.

Overlap integrals
=================

Overlap integrals effectively probe the overlap between two CGFs and are given by

.. math::

    S_{ij} = \left< \phi_{i} | \phi_{j} \right>

CGFs should be normalized and as such, their self-overlap should be equal to
1. In the code snippet below, the overlap matrix :math:`\mathbf{S}` is
calculated for a basis set composed of the two :math:`1s` atomic orbitals on H which
are separated by a distance of 1.4 Bohr.

.. code-block:: python

    from pyqint import PyQInt, CGF
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build CGF for a H atom located at the origin
    cgf1 = CGF([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF located 1.4 a.u. separated from CGF1
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
=================

Kinetic integrals determine the kinetic energy of a given orbital and are given
by

.. math::

    T_{ij} = \left< \phi_{i} \left| -\frac{1}{2} \nabla^{2} \right| \phi_{j} \right>

In the code snippet below, the kinetic energy matrix :math:`\mathbf{T}` is
calculated for a basis set composed of the two :math:`1s` atomic orbitals on H which
are separated by a distance of 1.4 Bohr.

.. code-block:: python

    from pyqint import PyQInt, CGF
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build CGF for a H atom located at the origin
    cgf1 = CGF([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF located 1.4 a.u. separated from CGF1
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
============================

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

    from pyqint import PyQInt, CGF
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build CGF for a H atom located at the origin
    cgf1 = CGF([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF located 1.4 a.u. separated from CGF1
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
======================

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

    from pyqint import PyQInt, CGF
    import numpy as np
    from copy import deepcopy

    # construct integrator object
    integrator = PyQInt()

    # build CGF for a H atom located at the origin
    cgf1 = CGF([0.0, 0.0, 0.0])

    cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

    # create a copy of the CGF located 1.4 a.u. separated from CGF1
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
=======================

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