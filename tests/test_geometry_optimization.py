import unittest
from pyqint import Molecule,GeometryOptimization
import numpy as np
import sys

class TestGeometryOptimization(unittest.TestCase):
    """
    Test the calculation of analytical gradients by comparing these
    with the finite difference results
    """

    def test_optimization_h2_sto3g(self):
        """
        Optimize dihydrogen molecule and assess that the energy corresponds to
        an optimized structure
        """
        # create new H2 molecule with perturbed geometry
        mol = Molecule()
        mol.add_atom('H', 0.9, 0.0, 0.0)
        mol.add_atom('H', -0.9, 0.0, 0.0)

        res = GeometryOptimization(verbose=False).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res['opt'].fun, -1.117506)

        self.assertEqual(len(res['energies']), len(res['forces']))
        self.assertEqual(len(res['energies']), len(res['coordinates']))
        self.assertEqual(res['coordinates'][0].shape, (len(mol.get_atoms()),3))

        # test existence of data object
        N = len(res['data']['cgfs'])
        self.assertEqual(N,2)
        self.assertEqual(res['energies'][-1], res['data']['energies'][-1])
        self.assertEqual(res['data']['overlap'].shape, (N, N))
        self.assertEqual(res['data']['kinetic'].shape, (N, N))
        self.assertEqual(res['data']['nuclear'].shape, (N, N))
        self.assertEqual(res['data']['tetensor'].shape, (N, N, N, N))

    @unittest.skip("Skipping H2 P321 test for receiving inconsistent results")
    def test_optimization_h2_p321(self):
        """
        Optimize dihydrogen molecule and assess that the energy corresponds to
        an optimized structure
        """
        # create new H2 molecule with perturbed geometry
        mol = Molecule()
        mol.add_atom('H', 0.9, 0.0, 0.0)
        mol.add_atom('H', -0.9, 0.0, 0.0)

        res = GeometryOptimization(verbose=False).run(mol, 'p321')
        np.testing.assert_almost_equal(res['opt'].fun, -1.1230, decimal=3)

        self.assertEqual(len(res['energies']), len(res['forces']))
        self.assertEqual(len(res['energies']), len(res['coordinates']))
        self.assertEqual(res['coordinates'][0].shape, (len(mol.get_atoms()),3))

        # test existence of data object
        N = len(res['data']['cgfs'])
        self.assertEqual(N,4)
        self.assertEqual(res['energies'][-1], res['data']['energies'][-1])
        self.assertEqual(res['data']['overlap'].shape, (N, N))
        self.assertEqual(res['data']['kinetic'].shape, (N, N))
        self.assertEqual(res['data']['nuclear'].shape, (N, N))
        self.assertEqual(res['data']['tetensor'].shape, (N, N, N, N))

    def test_optimization_ch4(self):
        """
        Optimize methane molecule and assess that the energy corresponds to
        an optimized structure
        """
        # create new CH4 molecule with perturbed geometry
        mol = Molecule()
        dist = 2.0/2
        mol.add_atom('C', 0.1, 0.0, 0.1, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

        res = GeometryOptimization(verbose=False).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res['opt'].fun, -39.72691085946399, decimal=4)

        self.assertEqual(len(res['energies']), len(res['forces']))
        self.assertEqual(len(res['energies']), len(res['coordinates']))
        self.assertEqual(res['coordinates'][0].shape, (len(mol.get_atoms()),3))

        # test existence of data object
        N = len(res['data']['cgfs'])
        self.assertEqual(N,9)
        self.assertEqual(res['energies'][-1], res['data']['energies'][-1])
        self.assertEqual(res['data']['overlap'].shape, (N, N))
        self.assertEqual(res['data']['kinetic'].shape, (N, N))
        self.assertEqual(res['data']['nuclear'].shape, (N, N))
        self.assertEqual(res['data']['tetensor'].shape, (N, N, N, N))

    def test_optimization_h2o(self):
        """
        Optimize the water molecule and assess that the energy corresponds to
        an optimized structure
        """
        # create new H2O molecule with perturbed geometry
        mol = Molecule()
        mol.add_atom('O', 0.0,  0.0, -0.2271310707, unit='angstrom')
        mol.add_atom('H', 0.0,  -0.8580158822, 0.5085242828, unit='angstrom')
        mol.add_atom('H', 0.0,  0.8580158822, 0.5085242828, unit='angstrom')

        res = GeometryOptimization(verbose=False).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res['opt'].fun, -74.96590347517174, decimal=4)

        self.assertEqual(len(res['energies']), len(res['forces']))
        self.assertEqual(len(res['energies']), len(res['coordinates']))
        self.assertEqual(res['coordinates'][0].shape, (len(mol.get_atoms()),3))

        # test existence of data object
        N = len(res['data']['cgfs'])
        self.assertEqual(N,7)
        self.assertEqual(res['energies'][-1], res['data']['energies'][-1])
        self.assertEqual(res['data']['overlap'].shape, (N, N))
        self.assertEqual(res['data']['kinetic'].shape, (N, N))
        self.assertEqual(res['data']['nuclear'].shape, (N, N))
        self.assertEqual(res['data']['tetensor'].shape, (N, N, N, N))

    @unittest.skipIf(sys.platform == "darwin",
                     "skipping test for MacOS")
    def test_optimization_c2h4(self):
        """
        Optimize ethylene molecule and assess that the energy corresponds to
        an optimized structure
        """
        # create new C2H4 molecule with perturbed geometry
        mol = Molecule()
        mol.add_atom('C', -0.5530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
        mol.add_atom('C',  0.4530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
        mol.add_atom('H', -1.3288875372, -0.9556191261 ,0.1000000000, unit='angstrom')
        mol.add_atom('H', -1.1288875372,  0.9356191261 ,0.0000000000, unit='angstrom')
        mol.add_atom('H',  1.3288875372,  0.9556191261 ,0.0000000000, unit='angstrom')
        mol.add_atom('H',  1.1288875372, -0.9156191261 ,0.1000000000, unit='angstrom')

        res = GeometryOptimization(verbose=False).run(mol, 'sto3g')
        np.testing.assert_almost_equal(res['opt'].fun, -77.07396213047552, decimal=4)

        self.assertEqual(len(res['energies']), len(res['forces']))
        self.assertEqual(len(res['energies']), len(res['coordinates']))
        self.assertEqual(res['coordinates'][0].shape, (len(mol.get_atoms()),3))

        # test existence of data object
        N = len(res['data']['cgfs'])
        self.assertEqual(N,14)
        self.assertEqual(res['energies'][-1], res['data']['energies'][-1])
        self.assertEqual(res['data']['overlap'].shape, (N, N))
        self.assertEqual(res['data']['kinetic'].shape, (N, N))
        self.assertEqual(res['data']['nuclear'].shape, (N, N))
        self.assertEqual(res['data']['tetensor'].shape, (N, N, N, N))

if __name__ == '__main__':
    unittest.main()
