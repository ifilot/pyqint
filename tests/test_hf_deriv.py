import unittest
from pyqint import CGF, HF, Molecule, MoleculeBuilder
from copy import deepcopy
import numpy as np

class TestHFDeriv(unittest.TestCase):

    def test_closed_shell_uhf_gradient_matches_rhf(self):
        solver = HF(MoleculeBuilder.from_name("h2o"), "sto3g")
        rhf = solver.rhf(calc_forces=True, tolerance_grad=1e-8)
        uhf = solver.uhf(multiplicity=1, calc_forces=True, tolerance_grad=1e-8)

        self.assertTrue(rhf["converged"] and uhf["converged"])
        np.testing.assert_allclose(uhf["forces"], rhf["forces"], atol=1e-8)
        self.assertIsNone(solver.uhf(multiplicity=1)["forces"])

    def test_open_shell_uhf_gradient_matches_energy_difference(self):
        mol = MoleculeBuilder.from_name("h2o")
        mol.set_charge(1)
        mol.get_atoms()[1][1][:] += np.array([0.12, -0.04, 0.08])
        options = dict(multiplicity=2, use_diis=False, tolerance=1e-13, tolerance_grad=1e-9)
        result = HF(mol, "sto3g", eri_storage="packed", eri_engine="auto").uhf(calc_forces=True, **options)
        self.assertTrue(result["converged"])
        orbitals = {"alpha": result["orbc_alpha"], "beta": result["orbc_beta"]}

        step = 1e-4
        expected = np.zeros_like(result["forces"])
        for i, j in np.ndindex(expected.shape):
            energies = []
            for sign in (-1, 1):
                displaced = deepcopy(mol)
                displaced.get_atoms()[i][1][j] += sign * step
                state = HF(displaced, "sto3g").uhf(orbc_init=orbitals, **options)
                self.assertTrue(state["converged"])
                energies.append(state["energy"])
            expected[i, j] = (energies[1] - energies[0]) / (2 * step)

        np.testing.assert_allclose(result["forces"], expected, rtol=1e-6, atol=5e-7)
        np.testing.assert_allclose(result["forces"].sum(axis=0), 0, atol=1e-10)

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

    def test_hartree_fock_forces_cgf_list(self):
        """
        Test that forces are also available when the basis set is supplied
        as a list of CGFs rather than by name
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0, unit='angstrom')
        mol.add_atom('H', -0.7570, 0.5860, 0.0, unit='angstrom')

        ref = HF(mol, 'sto3g').rhf(calc_forces=True)
        cgfs, _ = mol.build_basis('sto3g')
        res = HF(mol, cgfs).rhf(calc_forces=True)

        np.testing.assert_almost_equal(res['forces'], ref['forces'], decimal=10)

    def test_hartree_fock_forces_multicenter_salcs(self):
        """H2 gradients agree in the atom-centered and symmetry-adapted bases."""
        def calculate_salc_hf(mol, calc_forces=False, eri_storage="dense"):
            aos, _ = mol.build_basis("sto3g")
            salcs = []
            for phase in (1.0, -1.0):
                salc = CGF()
                for ao, weight in zip(aos, (1.0, phase)):
                    for gto in ao.gtos:
                        salc.add_gto_with_position(weight * gto.c, gto.p, gto.alpha, gto.l, gto.m, gto.n)
                salcs.append(salc)
            return HF(mol, salcs, eri_storage=eri_storage, eri_engine="auto").rhf(
                calc_forces=calc_forces, tolerance=1e-13)

        mol = Molecule()
        mol.add_atom("H", 0.0, 0.0, -0.7)
        mol.add_atom("H", 0.0, 0.0, 0.9)
        result = calculate_salc_hf(mol, calc_forces=True)
        reference = HF(mol, "sto3g").rhf(calc_forces=True)
        self.assertAlmostEqual(result["energy"], reference["energy"], places=12)
        np.testing.assert_allclose(result["forces"], reference["forces"], atol=1e-10)

        packed = calculate_salc_hf(mol, calc_forces=True, eri_storage="packed")
        self.assertTrue(packed["converged"])
        self.assertIsNone(packed["tetensor"])
        self.assertEqual(packed["packed_tei"].size, 6)
        for key in ("energy", "fock", "density", "erep", "ex", "forces"):
            np.testing.assert_allclose(packed[key], result[key], atol=1e-12, rtol=1e-11, err_msg=key)

        step = 1e-4
        plus, minus = deepcopy(mol), deepcopy(mol)
        plus.get_atoms()[0][1][2] += step
        minus.get_atoms()[0][1][2] -= step
        expected = (calculate_salc_hf(plus)["energy"] - calculate_salc_hf(minus)["energy"]) / (2.0 * step)
        self.assertAlmostEqual(result["forces"][0, 2], expected, places=8)

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
    """Reference energy gradient from a five-point central stencil."""
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

            forces[i, j] = dE

    return forces

if __name__ == '__main__':
    unittest.main()
