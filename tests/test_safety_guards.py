import unittest
from pyqint import PyQInt, Molecule, HF
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
        cgfs, _ = mol.build_basis('sto3g')

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

    def test_safety_uhf_multiplicity(self):
        """
        Test that impossible spin states are rejected by UHF
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        # multiplicity out of range (10 electrons)
        for multiplicity in (0, -1, 12):
            with self.assertRaises(ValueError):
                HF(mol, 'sto3g').uhf(multiplicity=multiplicity)

        # parity of the multiplicity does not match the number of electrons
        for multiplicity, nelec in ((2, None), (4, None), (1, 9), (3, 11)):
            with self.assertRaises(ValueError) as raises_cm:
                HF(mol, 'sto3g').uhf(multiplicity=multiplicity, nelec=nelec)
            self.assertIn('incompatible', raises_cm.exception.args[0])

    def test_safety_uhf_basis_size(self):
        """
        Test that UHF rejects spin states that do not fit in the basis set
        """
        # triplet He would require two alpha electrons in a single 1s
        # basis function
        mol = Molecule()
        mol.add_atom('He', 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError) as raises_cm:
            HF(mol, 'sto3g').uhf(multiplicity=3)
        self.assertIn('basis functions', raises_cm.exception.args[0])

if __name__ == '__main__':
    unittest.main()
