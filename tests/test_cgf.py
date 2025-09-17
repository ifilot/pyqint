import unittest
from pyqint import PyQInt, Molecule, CGF
import numpy as np
import os

class TestCGF(unittest.TestCase):

    def testFunctionsCGF(self):
        """
        Test getting amplitude from CGF
        """

        # construct integrator object
        integrator = PyQInt()

        # get compile info
        compile_info = integrator.get_compile_info()

        # build hydrogen molecule
        mol = Molecule()
        mol.add_atom('H', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')

        # test values at these coordinates
        coords = []
        for x in np.linspace(0, 10, 10):
            coords.append([x, x, x])

        # results
        ans = [6.2825e-01, 7.1229e-02, 6.8672e-03, 3.0000e-04, 3.7662e-06,
               1.3536e-08, 1.3927e-11, 4.1021e-15, 3.4591e-19, 8.3505e-24]

        # test for each coord
        amps = []
        for i,coord in enumerate(coords):
            amp = cgfs[0].get_amp(coord)
            amps.append(amp)

        np.testing.assert_almost_equal(amps, ans, 4)

    def testPlotGrid(self):
        """
        Test plotting of 1b2 molecular orbital of H2O
        """
        # coefficients
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

        ans = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_orb_1b2.npy'))
        np.testing.assert_almost_equal(res, ans, 6)

    def testSphericalHarmonicsOnSite(self):
        """
        Check if the overlap matrix of all spherical harmonics up to l=6 is the identity matrix
        """
        basis_functions = []
        p0 = [0.0, 0.0, 0.0]
        for l in range(7):
            for m in range(-l, l+1):
                orb = CGF(p0)
                orb.add_spherical_gto(1.0, 1.0, l, m)
                basis_functions.append(orb)
        # build overlap matrix
        om = np.zeros((len(basis_functions), len(basis_functions)))
        integrator = PyQInt()
        for i,orb1 in enumerate(basis_functions):
            for j,orb2 in enumerate(basis_functions):
                om[i,j] = integrator.overlap(orb1, orb2)
        np.testing.assert_almost_equal(om, np.eye(om.shape[0]), 6)

if __name__ == '__main__':
    unittest.main()
