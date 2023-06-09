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
In the script below, the six unique two-electron integrals for the H:sub:`2`
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

Basis sets and molecules
========================

Building molecules
------------------

Molecules can be efficiently built from the `Molecule` class. For example,
to build the H:sub:`2` molecule, one can run the script below.

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


Constructing basis functions for a molecule
-------------------------------------------

To construct the basis functions for a given molecule, one first needs to
construct the molecule after which the `build_basis` function can be used
to construct a basis.

The following basis sets are supported. For each basis set, the range of atoms
that are supported are given:

* `sto3g` (H-I)
* `sto6g` (H-Kr)
* `p321` (H-Cs)
* `p631` (H-Zn)

The example code below builds the basis functions for the H:sub:`2` molecule:

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

.. code-block:: python

    from pyqint import PyQInt, Molecule, HF
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    def main():
        # calculate sto-3g coefficients for h2o
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

    # coefficients (calculated by Hartree-Fock using a sto-3g basis set)
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
        # calculate sto-3g coefficients for h2o
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
