import unittest
from pyqint import Molecule, HF, COHP
import numpy as np

class TestCOHP(unittest.TestCase):

    def testCO(self):
        """
        Test construction of localized orbitals using Foster-Boys procedure
        for the CO molecule
        """
        d = 1.145414
        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, -d/2, unit='angstrom')
        mol.add_atom('O', 0.0, 0.0,  d/2, unit='angstrom')

        res = HF().rhf(mol, 'sto3g')
        cohp = COHP(res).run(res['orbc'], 0, 1)

        cohp_ref = np.array([
            0.0399,
            0.0104,
           -0.4365,
            0.2051,
           -0.2918,
           -0.2918,
            0.1098,
            0.5029,
            0.5029,
            6.4827
        ])

        # note that Foster-Boys optimization is somewhat random and thus
        # we use relatively loose testing criteria
        np.testing.assert_almost_equal(cohp,
                                       cohp_ref,
                                       decimal=4)

if __name__ == '__main__':
    unittest.main()
