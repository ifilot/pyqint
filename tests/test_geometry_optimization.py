import unittest
from pyqint import Molecule,GeometryOptimization
import numpy as np

class TestGeometryOptimization(unittest.TestCase):
    """
    Test the calculation of analytical gradients by comparing these
    with the finite difference results
    """

    def test_optimization_h2(self):
        # create new H2 molecule with slight adjustment in geometry and
        # seed the calculation with the previous converged result
        mol = Molecule()
        mol.add_atom('H', 0.9, 0.0, 0.0)
        mol.add_atom('H', -0.9, 0.0, 0.0)

        res = GeometryOptimization(mol, 'sto3g').run()
        np.testing.assert_almost_equal(res['energies'][-1], -1.117506)

    def test_optimization_ch4(self):
        # create new CH4 molecule with slight adjustment in geometry and
        # seed the calculation with the previous converged result
        mol = Molecule()
        dist = 2.0/2
        mol.add_atom('C', 0.0, 0.0, 0.0, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

        res = GeometryOptimization(mol, 'sto3g', verbose=False).run()
        np.testing.assert_almost_equal(res['energies'][-1], -39.7268635)

if __name__ == '__main__':
    unittest.main()
