.. index:: orbital_visualization

Orbital visualization
=====================

.. contents:: Table of Contents
    :depth: 3

Since orbitals are essentially three-dimensional scalar fields, there are two
useful procedures to visualize them. The scalar field can either be projected
onto a plane, creating so-called contour plots. Alternatively, a specific
value (i.e. the isovalue) of the scalar field can be chosen and all points in
space that have this value can be tied together creating a so-called isosurface.

Contour plots can be easily created using `matplotlib <https://matplotlib.org/>`_.
For the creation of isosurfaces, we use `PyTessel <https://ifilot.github.io/pytessel/>`_.

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
    * Isosurface generation requires the :program:`PyTessel` package to be
      installed. Make sure you have installed :program:`PyTessel` alongside 
      :program:`PyQInt`. For more details, see the :ref:`installation`.
    * Optionally, have a look at `PyTessel's documentation <https://ifilot.github.io/pytessel/>`_.

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

        result = HF(mol, 'sto3g').rhf()

        return result['cgfs'], result['orbc']

    if __name__ == '__main__':
        main()