import unittest
from pyqint import PyQInt, Molecule
import numpy as np
import os

class TestPlotData(unittest.TestCase):

    def test_plot_wavefunction(self):
        """
        Test generation of wave function
        """
        # build hydrogen molecule
        mol = Molecule("H2")
        mol.add_atom('H', 0.0, 0.00, -0.7)  # distances in Bohr lengths
        mol.add_atom('H', 0.0, 0.00, 0.7)   # distances in Bohr lengths
        cgfs, nuclei = mol.build_basis('sto3g')
        
        # construct integrator object
        integrator = PyQInt()
        
        # build grid points
        x = np.linspace(-2, 2, 6, endpoint=True)
        coord = np.flipud(np.vstack(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1)).T
        
        c = np.array([0.54893397, 0.54893397])
        wf = integrator.plot_wavefunction(coord, c, cgfs)
        
        res = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2_wf.npy'))
        np.testing.assert_almost_equal(wf, res, 4)
        
    def test_plot_gradient(self):
        """
        Test generation of wave function
        """
        # build hydrogen molecule
        mol = Molecule("H2")
        mol.add_atom('H', 0.0, 0.00, -0.7)  # distances in Bohr lengths
        mol.add_atom('H', 0.0, 0.00, 0.7)   # distances in Bohr lengths
        cgfs, nuclei = mol.build_basis('sto3g')
        
        # construct integrator object
        integrator = PyQInt()
                       
        c = np.array([0.54893397, 0.54893397])
        
        # test couple of single points
        sc = np.array([[-2.,-2.,-2.]])
        grad = integrator.plot_gradient(sc, c, cgfs)
        self.assertEqual(grad.shape[0], 1)
        self.assertEqual(grad.shape[1], 3)
        res = np.array([[0.00926244,0.00926244,0.0076777]])
        np.testing.assert_almost_equal(grad, res, 4)
        
        sc = np.array([[-2., 2.,-2.]])
        grad = integrator.plot_gradient(sc, c, cgfs)
        self.assertEqual(grad.shape[0], 1)
        self.assertEqual(grad.shape[1], 3)
        res = np.array([[0.00926244, -0.00926244, 0.0076777]])
        np.testing.assert_almost_equal(grad, res, 4)
        
        # test two points
        tc = np.array([[-2.,-2.,-2.],[-2., 2.,-2.]])
        grad = integrator.plot_gradient(tc, c, cgfs)
        self.assertEqual(grad.shape[0], 2)
        self.assertEqual(grad.shape[1], 3)
        res = np.array([[0.00926244,0.00926244,0.0076777],[0.00926244, -0.00926244, 0.0076777]])
        np.testing.assert_almost_equal(grad, res, 4)
        
        # test grid of points
        # x = np.linspace(-2, 2, 6, endpoint=True)
        # coord = np.flipud(np.vstack(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1)).T
        coord = integrator.build_rectgrid3d(-2, 2, 6, endpoint=True)
        grad = integrator.plot_gradient(coord, c, cgfs)
        self.assertEqual(grad.shape[0], 6*6*6)
        self.assertEqual(grad.shape[1], 3)        
        res = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2_grad.npy'))
        np.testing.assert_almost_equal(grad, res, 4)

if __name__ == '__main__':
    unittest.main()
