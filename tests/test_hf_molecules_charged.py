import unittest
from pyqint import Molecule, GeometryOptimization
import numpy as np

class TestHFMoleculeCharged(unittest.TestCase):

    def testHCO(self):
        """
        Test Hartree-Fock calculation for CO

        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('O', -0.12770367,  -0.00000000,   1.23419284)
        mol.add_atom('C', -0.50625665,   0.00000000,  -1.09149431)
        mol.add_atom('H',  1.57882331,  -0.00000000,  -2.06681794)
        mol.set_charge(-1)

        results = GeometryOptimization().run(mol, 'sto3g')

        # check that energy matches
        self.assertEqual(results['data']['nelec'], 16)    
        np.testing.assert_almost_equal(results['energies'][-1], -111.523147758727, 5)

if __name__ == '__main__':
    unittest.main()
