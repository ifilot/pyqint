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
        np.testing.assert_almost_equal(wf, res, 7)

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
        res = np.array([[0.009262441, 0.009262441, 0.007677702]])
        np.testing.assert_almost_equal(grad, res, 7)

        sc = np.array([[-2., 2.,-2.]])
        grad = integrator.plot_gradient(sc, c, cgfs)
        self.assertEqual(grad.shape[0], 1)
        self.assertEqual(grad.shape[1], 3)
        res = np.array([[ 0.009262441, -0.009262441, 0.007677702]])
        np.testing.assert_almost_equal(grad, res, 7)

        # test two points
        tc = np.array([[-2.,-2.,-2.],[-2., 2.,-2.]])
        grad = integrator.plot_gradient(tc, c, cgfs)
        self.assertEqual(grad.shape[0], 2)
        self.assertEqual(grad.shape[1], 3)
        res = np.array(
            [
                [0.009262441, 0.009262441, 0.007677702],
                [0.009262441, -0.009262441, 0.007677702],
            ]
        )
        np.testing.assert_almost_equal(grad, res, 7)

        # test grid of points
        # x = np.linspace(-2, 2, 6, endpoint=True)
        # coord = np.flipud(np.vstack(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1)).T
        coord = integrator.build_rectgrid3d(-2, 2, 6, endpoint=True)
        grad = integrator.plot_gradient(coord, c, cgfs)
        self.assertEqual(grad.shape[0], 6*6*6)
        self.assertEqual(grad.shape[1], 3)
        res = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2_grad.npy'))
        np.testing.assert_almost_equal(grad, res, 7)

    def test_plot_basis_functions(self):
        """
        Test that all basis functions evaluated in a single call match the
        basis functions evaluated one at a time
        """
        # build water molecule; 6-31G contains s- and p-functions on
        # several centers
        mol = Molecule("H2O")
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 1.43, 1.11)
        mol.add_atom('H', 0.0, -1.43, 1.11)
        cgfs, nuclei = mol.build_basis('p631')

        integrator = PyQInt()
        coord = integrator.build_rectgrid3d(-2, 2, 7, endpoint=True)

        amps = integrator.plot_basis_functions(coord, cgfs)
        self.assertEqual(amps.shape, (len(cgfs), len(coord)))
        for i, cgf in enumerate(cgfs):
            np.testing.assert_allclose(amps[i], integrator.plot_basis_function(coord, cgf),
                                       rtol=1e-14, atol=1e-14)

    def test_plot_basis_function_gradients(self):
        """
        Test that the gradients of all basis functions evaluated in a single
        call match the gradients evaluated one at a time
        """
        mol = Molecule("H2O")
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 1.43, 1.11)
        mol.add_atom('H', 0.0, -1.43, 1.11)
        cgfs, nuclei = mol.build_basis('p631')

        integrator = PyQInt()
        coord = integrator.build_rectgrid3d(-2, 2, 7, endpoint=True)

        grads = integrator.plot_basis_function_gradients(coord, cgfs)
        self.assertEqual(grads.shape, (len(cgfs), len(coord), 3))
        for i, cgf in enumerate(cgfs):
            np.testing.assert_allclose(grads[i], integrator.plot_gradient(coord, [1.0], [cgf]),
                                       rtol=1e-14, atol=1e-14)

    def test_plot_basis_functions_single_point(self):
        """
        Test that a single grid point given as a flat vector is accepted
        """
        mol = Molecule("H2")
        mol.add_atom('H', 0.0, 0.00, -0.7)
        mol.add_atom('H', 0.0, 0.00, 0.7)
        cgfs, nuclei = mol.build_basis('sto3g')

        integrator = PyQInt()
        point = np.array([0.1, 0.2, 0.3])
        amps = integrator.plot_basis_functions(point, cgfs)
        grads = integrator.plot_basis_function_gradients(point, cgfs)
        self.assertEqual(amps.shape, (2, 1))
        self.assertEqual(grads.shape, (2, 1, 3))
        for i, cgf in enumerate(cgfs):
            np.testing.assert_allclose(amps[i, 0], cgf.get_amp(point), rtol=1e-14)
            np.testing.assert_allclose(grads[i, 0], cgf.get_grad(point), rtol=1e-14)

if __name__ == '__main__':
    unittest.main()
