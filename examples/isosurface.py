from pyqint import PyQInt, Molecule, HF
import numpy as np
from pytessel import PyTessel
from purepytessel.api import isosurface_to_ply

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

    isosurface_to_ply(scalarfield, unitcell, isovalue, filename.replace('04','05'))

def calculate_co():
    mol = Molecule()
    mol.add_atom('C', 0.0, -0.5, 0.0)
    mol.add_atom('O', 0.0, 0.5, 0.0)

    result = HF(mol, 'sto3g').rhf()

    return result['cgfs'], result['orbc']

if __name__ == '__main__':
    main()