import unittest
from pyqint import MoleculeBuilder, HF
import numpy as np

class TestEnergyDecomposition(unittest.TestCase):

    def test_hartree_fock_h2o(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = MoleculeBuilder().from_name("h2o")

        res = HF().rhf(mol, 'sto3g')
        P = res['density']
        T = res['kinetic']
        V = res['nuclear']
        H = res['fock']
        enucrep = res['enucrep']
        energy = res['energy']
        
        np.testing.assert_almost_equal(0.5 * np.einsum('ji,ij', P, T+V+H) + enucrep, energy, decimal=16)
        
if __name__ == '__main__':
    unittest.main()