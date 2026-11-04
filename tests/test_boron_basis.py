import unittest

import numpy as np

from pyqint import HF, Molecule


class TestBoronBasis(unittest.TestCase):
    def test_boron_sto6g(self):
        """Check the core and valence contractions against an atomic reference."""
        mol = Molecule('B')
        mol.add_atom('B', 0.0, 0.0, 0.0)
        result = HF(mol, 'sto6g').uhf(multiplicity=2, tolerance=1e-11)
        # PySCF 2.14, UHF/STO-6G. PyQInt stores rounded coefficients.
        self.assertAlmostEqual(result['energy'], -24.3942945593718, places=6)

    def test_borane_sto6g(self):
        """Exercise the corrected boron entry in a closed-shell molecule."""
        mol = Molecule('BH3')
        mol.add_atom('B', 0.0, 0.0, 0.0)
        for angle in (0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0):
            mol.add_atom('H', 2.19 * np.cos(angle), 2.19 * np.sin(angle), 0.0)
        result = HF(mol, 'sto6g').rhf(tolerance=1e-11)
        # PySCF 2.14, RHF/STO-6G, at the same geometry in bohr.
        self.assertAlmostEqual(result['energy'], -26.3311297683041, places=4)
