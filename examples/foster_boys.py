from pyqint import Molecule, HF, PyQInt, FosterBoys
import pyqint
import numpy as np
from pytessel import PyTessel

#
# Plot the isosurfaces for the CO molecule
#

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
