import unittest
from pyqint import Molecule, GeometryOptimization
import numpy as np

class TestHFMoleculeCharged(unittest.TestCase):

    def testHCO(self):
        """
        Test geometry optimization of the formyl anion HCO-

        At the optimized geometry, the energy agrees with PySCF, using the
        same basis set, to within 1e-10 Ht
        """
        mol = Molecule()
        mol.add_atom('O', -0.12770367,  -0.00000000,   1.23419284)
        mol.add_atom('C', -0.50625665,   0.00000000,  -1.09149431)
        mol.add_atom('H',  1.57882331,  -0.00000000,  -2.06681794)
        mol.set_charge(-1)

        results = GeometryOptimization(mol, 'sto3g').run()

        # check that energy matches
        self.assertEqual(results['data']['nelec'], 16)
        np.testing.assert_almost_equal(results['energies'][-1], -111.52314846187623, 5)

if __name__ == '__main__':
    unittest.main()
