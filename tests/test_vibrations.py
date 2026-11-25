import unittest

import numpy as np

from pyqint import CGF, HF, Molecule, MoleculeBuilder
from pyqint.methods.vibrations import WAVENUMBER_FACTOR, harmonic_analysis


class TestVibrations(unittest.TestCase):

    def test_h2_stretch_matches_energy_curvature(self):
        def molecule(distance):
            mol = Molecule()
            mol.add_atom("H", 0.0, 0.0, -distance / 2)
            mol.add_atom("H", 0.0, 0.0, distance / 2)
            return mol

        distance = 1.346
        step = 1e-3
        mol = molecule(distance)
        original = np.array([position.copy() for _, position in mol])
        masses = np.array([1.0, 2.0])
        result = HF(mol, "sto3g", eri_storage="packed", eri_engine="auto").vibrational_frequencies(
            step=step, masses=masses,
        )
        energy = [HF(molecule(distance + shift), "sto3g").rhf()["energy"] for shift in (-step, 0, step)]
        spring = (energy[0] - 2 * energy[1] + energy[2]) / step ** 2
        reduced_mass = 2 / 3

        np.testing.assert_allclose(result["frequencies"], [WAVENUMBER_FACTOR * np.sqrt(spring / reduced_mass)],
                                   rtol=2e-5)
        np.testing.assert_allclose(result["hessian"], result["hessian"].T, atol=0)
        np.testing.assert_allclose(np.linalg.norm(result["normal_modes"]), 1)
        np.testing.assert_allclose(masses @ result["normal_modes"][0], 0, atol=1e-12)
        np.testing.assert_array_equal([position for _, position in mol], original)

    def test_uhf_h2_cation_stretch_matches_energy_curvature(self):
        def solver(distance):
            mol = Molecule()
            mol.add_atom("H", 0, 0, -distance / 2)
            mol.add_atom("H", 0, 0, distance / 2)
            mol.set_charge(1)
            return HF(mol, "sto3g")

        distance, step = 2.0, 1e-3
        result = solver(distance).vibrational_frequencies(multiplicity=2, step=step)
        energies = [solver(distance + shift).uhf(multiplicity=2)["energy"] for shift in (-step, 0, step)]
        spring = (energies[0] - 2 * energies[1] + energies[2]) / step ** 2
        reduced_mass = result["masses"][0] / 2
        np.testing.assert_allclose(result["frequencies"], [WAVENUMBER_FACTOR * np.sqrt(spring / reduced_mass)],
                                   rtol=2e-5)

    def test_bent_uhf_frequencies_against_reference(self):
        mol = Molecule("water cation")
        mol.add_atom("O", 0.0, 0.0, 0.259)
        mol.add_atom("H", 0.12, 1.415, -0.954)
        mol.add_atom("H", 0.0, -1.455, -1.034)
        mol.set_charge(1)
        result = HF(mol, "sto3g", eri_storage="packed", eri_engine="auto").vibrational_frequencies(
            multiplicity=2, use_diis=False,
        )

        # PySCF 2.15.1 analytic UHF Hessian using PyQInt's basis and masses.
        expected = [1746.908169, 3909.969941, 4409.345016]
        np.testing.assert_allclose(result["frequencies"], expected, rtol=0, atol=0.05)
        self.assertEqual(result["normal_modes"].shape, (3, 3, 3))

    def test_multicenter_salc_stretch(self):
        mol = MoleculeBuilder.from_name("h2")
        cgfs, _ = mol.build_basis("sto3g")
        basis = []
        for phase in (1, -1):
            salc = CGF()
            for cgf, coefficient in zip(cgfs, (1, phase)):
                for gto in cgf.gtos:
                    salc.add_gto_with_position(coefficient * gto.c, gto.p, gto.alpha, gto.l, gto.m, gto.n)
            basis.append(salc)

        reference = HF(mol, "sto3g").vibrational_frequencies()
        reduced = HF(mol, basis).vibrational_frequencies()
        np.testing.assert_allclose(reduced["hessian"], reference["hessian"], atol=1e-8)
        np.testing.assert_allclose(reduced["frequencies"], reference["frequencies"], atol=1e-5)

    def test_nonlinear_projection_and_imaginary_modes(self):
        coordinates = np.array([[0., 0., 0.], [0., 0., 1.], [0., 2., 0.]])
        masses = np.array([16., 1., 1.])
        hessian = -np.diag(np.repeat(masses, 3))
        frequencies, modes = harmonic_analysis(hessian, coordinates, masses)
        np.testing.assert_allclose(frequencies, [-WAVENUMBER_FACTOR] * 3)
        self.assertEqual(modes.shape, (3, 3, 3))
        np.testing.assert_allclose(np.einsum("a,mac->mc", masses, modes), 0, atol=1e-12)
