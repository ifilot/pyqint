from pyqint import PyQInt, Molecule, HF
import numpy as np
from pytessel import PyTessel

def main():
    # calculate the sto3g molecular orbitals of CO
    cgfs, coeff = calculate_co()

    # build the isosurfaces of the fifth molecular orbital at an isovalue
    # of 0.1 and store these as .ply files
    build_isosurface('co_mo5', cgfs, coeff[:,4], 0.1)

def build_isosurface(filename, cgfs, coeff, isovalue):
    """
    Construct the isosurfaces of a molecular orbital at +isovalue and
    -isovalue, i.e. its positive and negative lobes, and store these as
    <filename>_pos.ply and <filename>_neg.ply
    """
    # evaluate the molecular orbital on a cubic grid
    sz = 100
    integrator = PyQInt()
    grid = integrator.build_rectgrid3d(-5, 5, sz)
    scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff, cgfs), (sz, sz, sz))
    unitcell = np.diag(np.ones(3) * 10.0)

    pytessel = PyTessel()
    for suffix, value in (('pos', isovalue), ('neg', -isovalue)):
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), value)
        pytessel.write_ply('%s_%s.ply' % (filename, suffix), vertices, normals, indices)

def calculate_co():
    mol = Molecule()
    mol.add_atom('C', 0.0, -0.5, 0.0)
    mol.add_atom('O', 0.0, 0.5, 0.0)

    result = HF(mol, 'sto3g').rhf()

    return result['cgfs'], result['orbc']

if __name__ == '__main__':
    main()
