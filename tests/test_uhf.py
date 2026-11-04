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
        np.testing.assert_almost_equal(results_rhf['energy'], results_uhf['energy'], 7)

        # verify that terms are being calculated
        np.testing.assert_almost_equal(results_uhf['orbe_alpha'], results_uhf['orbe_beta'], decimal=5)
        np.testing.assert_almost_equal(results_uhf['orbe_alpha'], results_rhf['orbe'], decimal=5)

        # a closed-shell singlet is free of spin contamination
        self.assertEqual(results_uhf['s2_exact'], 0.0)
        np.testing.assert_almost_equal(results_uhf['s2'], 0.0, decimal=8)

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

        # reference values obtained with PySCF 2.14 (UHF/STO-3G, Cartesian);
        # the small energy difference stems from the basis set coefficients
        np.testing.assert_almost_equal(results_uhf['energy'], -39.0767088890, decimal=4)
        self.assertEqual(results_uhf['s2_exact'], 0.75)
        np.testing.assert_almost_equal(results_uhf['s2'], 0.765224, decimal=4)

    def test_unrestricted_hartree_fock_n2_anion(self):
        """
        Test unrestricted Hartree-Fock calculation on the N2- anion
        using p631 basis set.
        """
        R = 2.074

        mol = Molecule()
        mol.add_atom('N', 0.0, 0.0, 0.0)
        mol.add_atom('N', 0.0, 0.0, R)

        # N2- (15 electrons) is a doublet → multiplicity = 2
        results_uhf = HF(mol, 'p631').uhf(nelec=15, multiplicity=2, tolerance=1e-14)
        self.assertEqual(results_uhf['nalpha'], 8)
        self.assertEqual(results_uhf['nbeta'], 7)

        # reference values obtained with PySCF 2.14 (UHF/6-31G, Cartesian)
        np.testing.assert_almost_equal(results_uhf['energy'], -108.7504624176, decimal=4)
        np.testing.assert_almost_equal(results_uhf['s2'], 0.757046, decimal=4)

    def test_unrestricted_hartree_fock_atoms(self):
        """
        Test unrestricted Hartree-Fock calculations on doublet atoms against
        reference values obtained with PySCF 2.14 (UHF/STO-3G)
        """
        for element, energy_ref in (('H', -0.4665818496), ('Li', -7.3155259813)):
            mol = Molecule()
            mol.add_atom(element, 0.0, 0.0, 0.0)

            results_uhf = HF(mol, 'sto3g').uhf(multiplicity=2, tolerance=1e-12)
            np.testing.assert_almost_equal(results_uhf['energy'], energy_ref, decimal=5)

            # a single unpaired electron outside a closed shell gives a
            # pure doublet
            np.testing.assert_almost_equal(results_uhf['s2'], 0.75, decimal=8)

    def test_unrestricted_hartree_fock_electron_count(self):
        """
        Test the distribution of electrons over the alpha and beta channels
        for various multiplicities and electron counts
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        for nelec, multiplicity, nalpha, nbeta in ((10, 1, 5, 5),
                                                   (10, 3, 6, 4),
                                                   (10, 5, 7, 3),
                                                   ( 9, 2, 5, 4),
                                                   (11, 2, 6, 5)):
            res = HF(mol, 'sto3g').uhf(multiplicity=multiplicity, nelec=nelec)
            self.assertEqual(res['nelec'], nelec)
            self.assertEqual(res['nalpha'], nalpha)
            self.assertEqual(res['nbeta'], nbeta)
            self.assertEqual(res['multiplicity'], multiplicity)

            # the total and spin densities integrate to the electron counts
            S = res['overlap']
            np.testing.assert_almost_equal(np.trace(res['density'] @ S), nelec)
            np.testing.assert_almost_equal(
                np.trace((res['density_alpha'] - res['density_beta']) @ S),
                nalpha - nbeta)

            # <S^2> cannot be lower than the exact value S(S+1)
            S_ = 0.5 * (multiplicity - 1)
            self.assertEqual(res['s2_exact'], S_ * (S_ + 1))
            self.assertGreaterEqual(res['s2'], res['s2_exact'] - 1e-8)

if __name__ == '__main__':
    unittest.main()
