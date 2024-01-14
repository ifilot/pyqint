import unittest
from pyqint import PyQInt, Molecule
import numpy as np

class TestSafetyGuards(unittest.TestCase):
    """
    Tests certain safety features preventing user errors
    """

    def test_safety_plot_wavefunction(self):
        """
        Test for correct array sizes
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
        grid = np.flipud(np.vstack(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1)).T
        
        # put in matrix object -> should yield error
        c = np.identity(2)
        with self.assertRaises(TypeError) as raises_cm:
            integrator.plot_wavefunction(grid, c, cgfs)

        exception = raises_cm.exception
        self.assertTrue('arrays can be converted to Python scalars' in exception.args[0])
        
        # put in matrix object -> should yield error
        c = np.identity(3)
        with self.assertRaises(Exception) as raises_cm:
            integrator.plot_wavefunction(grid, c, cgfs)

        exception = raises_cm.exception
        self.assertEqual(exception.args, ('Dimensions of cgf list and coefficient matrix do not match (2 != 3)',))
        
        # put in matrix object -> should yield error
        c = np.ones(3)
        with self.assertRaises(Exception) as raises_cm:
            integrator.plot_wavefunction(grid, c, cgfs)
            
        exception = raises_cm.exception
        self.assertEqual(exception.args, ('Dimensions of cgf list and coefficient matrix do not match (2 != 3)',))
        
    def test_safety_plot_gradient(self):
        """
        Test for correct array sizes
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
        grid = np.flipud(np.vstack(np.meshgrid(x, x, x, indexing='ij')).reshape(3,-1)).T
        
        # put in matrix object -> should yield error
        c = np.identity(2)
        with self.assertRaises(TypeError) as raises_cm:
            integrator.plot_gradient(grid, c, cgfs)

        exception = raises_cm.exception
        self.assertEqual(type(exception), TypeError)
        
        # put in matrix object -> should yield error
        c = np.identity(3)
        with self.assertRaises(Exception) as raises_cm:
            integrator.plot_gradient(grid, c, cgfs)

        exception = raises_cm.exception
        self.assertEqual(exception.args, ('Dimensions of cgf list and coefficient matrix do not match (2 != 3)',))
        
        # put in matrix object -> should yield error
        c = np.ones(3)
        with self.assertRaises(Exception) as raises_cm:
            integrator.plot_gradient(grid, c, cgfs)
            
        exception = raises_cm.exception
        self.assertEqual(exception.args, ('Dimensions of cgf list and coefficient matrix do not match (2 != 3)',))

if __name__ == '__main__':
    unittest.main()
