import unittest
from pyqint import Molecule, HF
import numpy as np

class TestHF(unittest.TestCase):

    def testHartreeFockWater(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        results = perform_hf(mol)

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -73.21447132, 4)

def perform_hf(mol):
    results = HF().rhf(mol, 'sto3g', tolerance=1e-9)
    return results

if __name__ == '__main__':
    unittest.main()
