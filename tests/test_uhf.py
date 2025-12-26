import unittest
from pyqint import Molecule, HF
import numpy as np

class TestHF(unittest.TestCase):

    def test_hartree_fock_h2o(self):
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
        np.testing.assert_almost_equal(results_rhf['energy'], results_uhf['energy'], 5)

        # verify that terms are being calculated
        np.testing.assert_almost_equal(results_uhf['orbe_alpha'], results_uhf['orbe_beta'], decimal=5)
        np.testing.assert_almost_equal(results_uhf['orbe_alpha'], results_rhf['orbe'], decimal=5)

if __name__ == '__main__':
    unittest.main()
