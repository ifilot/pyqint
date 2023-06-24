import unittest
from pyqint import Molecule,GeometryOptimization
import numpy as np

class TestGeometryOptimization(unittest.TestCase):
    """
    Test the calculation of analytical gradients by comparing these
    with the finite difference results
    """

    def test_optimization_h2_sto3g(self):
        # create new H2 molecule with slight adjustment in geometry
        mol = Molecule()
        mol.add_atom('H', 0.9, 0.0, 0.0)
        mol.add_atom('H', -0.9, 0.0, 0.0)

        res = GeometryOptimization(verbose=False).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res.fun, -1.117506)

    def test_optimization_h2_p321(self):
        # create new H2 molecule with slight adjustment in geometry
        mol = Molecule()
        mol.add_atom('H', 0.9, 0.0, 0.0)
        mol.add_atom('H', -0.9, 0.0, 0.0)

        res = GeometryOptimization(verbose=False).run(mol, 'p321')
        np.testing.assert_almost_equal(res.fun, -1.1229598312004625, decimal=4)

    def test_optimization_ch4(self):
        # create new CH4 molecule with slight adjustment in geometry
        mol = Molecule()
        dist = 2.0/2
        mol.add_atom('C', 0.1, 0.0, 0.1, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

        res = GeometryOptimization(verbose=False).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res.fun, -39.7268635)

    def test_optimization_h2o(self):
        # create new H2O molecule with slight adjustment in geometry
        mol = Molecule()
        mol.add_atom('O', 0.0,  0.0, -0.2271310707, unit='angstrom')
        mol.add_atom('H', 0.0,  -0.8580158822, 0.5085242828, unit='angstrom')
        mol.add_atom('H', 0.0,  0.8580158822, 0.5085242828, unit='angstrom')

        res = GeometryOptimization(verbose=False).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res.fun, -74.96590013836217)

    def test_optimization_c2h4(self):
        # create new C2H4 molecule with slight adjustment in geometry
        mol = Molecule()
        mol.add_atom('C', -0.5530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
        mol.add_atom('C',  0.4530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
        mol.add_atom('H', -1.3288875372, -0.9556191261 ,0.1000000000, unit='angstrom')
        mol.add_atom('H', -1.1288875372,  0.9356191261 ,0.0000000000, unit='angstrom')
        mol.add_atom('H',  1.3288875372,  0.9556191261 ,0.0000000000, unit='angstrom')
        mol.add_atom('H',  1.1288875372, -0.9156191261 ,0.1000000000, unit='angstrom')

        res = GeometryOptimization(verbose=True).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res.fun, -77.0739544)

if __name__ == '__main__':
    unittest.main()
