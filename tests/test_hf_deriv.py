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
            res = HF(mol, 'sto3g').rhf(calc_forces=True)

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
        res = HF(mol, 'sto3g').rhf(calc_forces=True)

        # calculate forces using finite difference
        forces = calculate_forces_finite_difference(mol)

        np.testing.assert_allclose(res['forces'], forces, rtol=1e-5, atol=1e-3)

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
        res = HF(mol, 'sto3g').rhf(calc_forces=True, tolerance=1e-12)

        # calculate forces using finite difference
        forces = calculate_forces_finite_difference(mol)

        np.testing.assert_allclose(res['forces'], forces, rtol=1e-5, atol=1e-3)

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
        res = HF(mol, 'sto3g').rhf(calc_forces=True)

        # calculate forces using finite difference
        forces = calculate_forces_finite_difference(mol)

        np.testing.assert_allclose(res['forces'], forces, rtol=1e-5, atol=1e-3)

def perform_hf(mol):
    sol = HF(mol, 'sto3g').rhf(tolerance=1e-15)
    return sol

def calculate_forces_finite_difference(mol):
    """
    Finite-difference reference forces using a 5-point central stencil.

    Forces are defined as:
        F = - dE / dR
    """
    forces = np.zeros((len(mol.get_atoms()), 3), dtype=float)

    h = 2e-2   # larger step is fine with 5-point stencil

    for i in range(len(mol.get_atoms())):
        for j in range(3):

            # R - 2h
            mol_m2 = deepcopy(mol)
            mol_m2.get_atoms()[i][1][j] -= 2*h
            e_m2 = perform_hf(mol_m2)['energy']

            # R - h
            mol_m1 = deepcopy(mol)
            mol_m1.get_atoms()[i][1][j] -= h
            e_m1 = perform_hf(mol_m1)['energy']

            # R + h
            mol_p1 = deepcopy(mol)
            mol_p1.get_atoms()[i][1][j] += h
            e_p1 = perform_hf(mol_p1)['energy']

            # R + 2h
            mol_p2 = deepcopy(mol)
            mol_p2.get_atoms()[i][1][j] += 2*h
            e_p2 = perform_hf(mol_p2)['energy']

            # 5-point central difference for gradient
            dE = (
                -e_p2
                + 8.0*e_p1
                - 8.0*e_m1
                + e_m2
            ) / (12.0*h)

            # Convert gradient → force
            forces[i, j] = dE

    return forces

if __name__ == '__main__':
    unittest.main()
