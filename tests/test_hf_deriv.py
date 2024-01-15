import unittest
from pyqint import Molecule, HF
from copy import deepcopy
import numpy as np

class TestHFDeriv(unittest.TestCase):

    def test_hartree_fock_forces_h2(self):
        """
        Test Hartree-Fock calculation on H2 using STO-3G basis set
        """
        for i in range(0,5):
            mol = Molecule()
            mol.add_atom('H', 0.0000, 0.0000, 0.30 + i * 0.50, unit='angstrom')
            mol.add_atom('H', 0.0000, 0.0000, -0.30 - i * 0.50, unit='angstrom')

            # calculate forces using analytical derivatives
            solver = HF()
            res = solver.rhf(mol, 'sto3g', calc_forces=True)

            # calculate forces using finite difference
            forces = calculate_forces_finite_difference(mol)

            np.testing.assert_almost_equal(res['forces'], forces, decimal=4)

    def test_hartree_fock_forces_h2o(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.1, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        # calculate forces using analytical derivatives
        solver = HF()
        res = solver.rhf(mol, 'sto3g', calc_forces=True)

        # calculate forces using finite difference
        forces = calculate_forces_finite_difference(mol)

        np.testing.assert_almost_equal(res['forces'], forces, decimal=4)

    def test_hartree_fock_forces_ch4(self):
        """
        Test Hartree-Fock calculation on CH4 using STO-3G basis set

        Note that the CH4 molecule is slightly perturbed
        """
        mol = Molecule()
        mol.add_atom('C', 0.0,  0.1, 0.0, unit='angstrom')
        mol.add_atom('H', 0.1,  0.15, 1.0830098121, unit='angstrom')
        mol.add_atom('H', 0.0,  -1.0210714424, -0.3610032723, unit='angstrom')
        mol.add_atom('H', -0.8842738057,  0.5105357237, -0.3610032748, unit='angstrom')
        mol.add_atom('H', 0.8842738117,  0.5105357203, -0.3610032651, unit='angstrom')

        # calculate forces using analytical derivatives
        solver = HF()
        res = solver.rhf(mol, 'sto3g', calc_forces=True, tolerance=1e-12)

        # calculate forces using finite difference
        forces = calculate_forces_finite_difference(mol)

        np.testing.assert_almost_equal(res['forces'], forces, decimal=3)

    def test_hartree_fock_forces_co2(self):
        """
        Test Hartree-Fock calculation on CH4 using STO-3G basis set

        Note that the CO2 molecule is slightly perturbed
        """
        mol = Molecule()
        mol.add_atom('C', 0.0,  0.0, 0.0, unit='angstrom')
        mol.add_atom('O', 0.0,  0.0, 1.2879700928, unit='angstrom')
        mol.add_atom('O', 0.0,  0.0, -1.2879700928, unit='angstrom')

        # calculate forces using analytical derivatives
        solver = HF()
        res = solver.rhf(mol, 'sto3g', calc_forces=True)

        # calculate forces using finite difference
        forces = calculate_forces_finite_difference(mol)

        np.testing.assert_almost_equal(res['forces'], forces, decimal=3)

def perform_hf(mol):
    sol = HF().rhf(mol, 'sto3g', tolerance=1e-12)
    return sol

def calculate_forces_finite_difference(mol):
    """
    Calculates the forces on each of the atoms using a finite difference
    approach.
    """
    forces = np.zeros((len(mol.get_atoms()),3))

    sz = 1e-3

    for i in range(0, len(mol.get_atoms())): # loop over nuclei
        for j in range(0, 3): # loop over directions
            mol1 = deepcopy(mol)
            mol1.get_atoms()[i][1][j] -= sz / 2
            mol2 = deepcopy(mol)
            mol2.get_atoms()[i][1][j] += sz / 2

            energy1 = perform_hf(mol1)['energy']
            energy2 = perform_hf(mol2)['energy']

            forces[i,j] = (energy2 - energy1) / sz

    return forces

if __name__ == '__main__':
    unittest.main()
