import unittest
from unittest.mock import patch

from pyqint import HF, MoleculeBuilder, PyQInt


class TestMemoryPreflight(unittest.TestCase):

    def setUp(self):
        self.mol = MoleculeBuilder.from_name("h2")
        self.cgfs, self.nuclei = self.mol.build_basis("sto3g")

    def test_rejects_integral_builds(self):
        integrator = PyQInt()
        with patch("pyqint._memory._available_memory", return_value=0):
            with self.assertRaisesRegex(MemoryError, "Dense ERI build.*available.*packed"):
                integrator.build_tei_openmp(iter(self.cgfs), profile=True)
            with self.assertRaisesRegex(MemoryError, "Packed ERI build.*available"):
                integrator.build_tei_openmp(self.cgfs, packed=True)
            with self.assertRaisesRegex(MemoryError, "Packed ERI build.*available"):
                integrator.build_integrals_openmp(self.cgfs, self.nuclei, packed=True)
            with self.assertRaisesRegex(MemoryError, "Geometric derivative build.*available"):
                integrator.build_geometric_derivatives_openmp(self.cgfs, self.nuclei)
            with self.assertRaisesRegex(MemoryError, "Dipole integral build.*available"):
                integrator.build_dipole_openmp(self.cgfs)

    def test_checks_force_memory_before_scf(self):
        solver = HF(self.mol, "sto3g", eri_storage="packed")
        with patch("pyqint._memory._available_memory", return_value=0):
            with self.assertRaisesRegex(MemoryError, "Geometric derivative build"):
                solver.rhf(calc_forces=True)
            with self.assertRaisesRegex(MemoryError, "Geometric derivative build"):
                solver.uhf(multiplicity=1, calc_forces=True)

    def test_bounds_automatic_quartet_planning(self):
        integrator = PyQInt(eri_engine="auto")
        with patch("pyqint.pyqint_core.check_eri_memory", return_value=0):
            with self.assertRaisesRegex(MemoryError, "quartet planning.*workspace"):
                integrator.build_tei_openmp(self.cgfs, packed=True)
        with patch("pyqint._memory._available_memory", return_value=None):
            tensor = integrator.build_tei_openmp(self.cgfs, packed=True)
        self.assertEqual(tensor.shape, (6,))
