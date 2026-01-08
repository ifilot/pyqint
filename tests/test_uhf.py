import unittest
from pyqint import Molecule, HF
import numpy as np

class TestHF(unittest.TestCase):

    def test_unrestricted_hartree_fock_h2o(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        results_rhf = HF(mol, 'sto3g').rhf(tolerance=1e-12)
        results_uhf = HF(mol, 'sto3g').uhf(multiplicity=1, tolerance=1e-12)

        # check that energy matches
        np.testing.assert_almost_equal(results_rhf['energy'], results_uhf['energy'], 9)

        # verify that terms are being calculated
        np.testing.assert_almost_equal(results_uhf['orbe_alpha'], results_uhf['orbe_beta'], decimal=5)
        np.testing.assert_almost_equal(results_uhf['orbe_alpha'], results_rhf['orbe'], decimal=5)

    def test_unrestricted_hartree_fock_ch3(self):
        """
        Test unrestricted Hartree-Fock calculation on the methyl radical (CH3)
        using STO-3G basis set.

        Geometry:
            - Planar CH3
            - C-H bond length = 2.039 a.u.
            - H-C-H angles = 120 degrees
        """
        R = 2.039
        sqrt3 = np.sqrt(3.0)

        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, 0.0)
        mol.add_atom('H',  R, 0.0, 0.0)
        mol.add_atom('H', -0.5 * R,  0.5 * sqrt3 * R, 0.0)
        mol.add_atom('H', -0.5 * R, -0.5 * sqrt3 * R, 0.0)

        # CH3 is a doublet → multiplicity = 2
        results_uhf = HF(mol, 'sto3g').uhf(multiplicity=2, tolerance=1e-12)

        # basic sanity checks
        self.assertIn('energy', results_uhf)
        self.assertIn('orbe_alpha', results_uhf)
        self.assertIn('orbe_beta', results_uhf)

        # alpha and beta orbital energies should NOT be identical for
        # an open-shell system
        with self.assertRaises(AssertionError):
            np.testing.assert_almost_equal(
                results_uhf['orbe_alpha'],
                results_uhf['orbe_beta'],
                decimal=6
            )

        # energies should be finite and negative
        self.assertTrue(np.isfinite(results_uhf['energy']))
        self.assertLess(results_uhf['energy'], 0.0)

    def test_unrestricted_hartree_fock_n2(self):
        """
        Test unrestricted Hartree-Fock calculation on the N2 molecule
        using p631 basis set.
        """
        R = 2.074

        mol = Molecule()
        mol.add_atom('N', 0.0, 0.0, 0.0)
        mol.add_atom('N', 0.0, 0.0, R)

        # N2 is a doublet → multiplicity = 2
        results_uhf = HF(mol, 'p631').uhf(nelec=15, multiplicity=2, tolerance=1e-14)
        np.testing.assert_almost_equal(results_uhf['energy'], -108.75042559940775, decimal=4)

if __name__ == '__main__':
    unittest.main()
