import unittest
from pyqint import Molecule, HF, FosterBoys
import numpy as np

class TestFosterBoys(unittest.TestCase):

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

        # note that a seed is given here for reproducibility purposes
        res_fb = FosterBoys(res, seed=0).run(nr_runners=5)

        orbe_ref = np.array([
            -20.30750217,
            -11.0370294,
            -0.83093927,
            -0.8309353,
            -0.83084896,
            -0.81363734,
            -0.52411525,
            res['orbe'][7],
            res['orbe'][8],
            res['orbe'][9]
        ])

        # note that Foster-Boys optimization is somewhat random and thus
        # we use relatively loose testing criteria
        np.testing.assert_almost_equal(res_fb['orbe'],
                                       orbe_ref,
                                       decimal=2)

    def testCH4(self):
        """
        Test construction of localized orbitals using Foster-Boys procedure
        for the CH4 molecule
        """
        mol = Molecule()
        dist = 1.78/2
        mol.add_atom('C', 0.0, 0.0, 0.0, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

        res = HF().rhf(mol, 'sto3g')

        # note that a seed is given here for reproducibility purposes
        res_fb = FosterBoys(res, seed=0).run(nr_runners=5)

        orbe_ref = np.array([
            -11.050113,
            -0.47136919,
            -0.47136899,
            -0.47136892,
            -0.47136892,
            res['orbe'][5],
            res['orbe'][6],
            res['orbe'][7],
            res['orbe'][8]
        ])

        # assert orbital energies
        np.testing.assert_almost_equal(res_fb['orbe'], orbe_ref, decimal=2)

        # specifically test for quadruple degenerate orbital
        for i in range(0,4):
            for j in range(i+1,4):
                # note that Foster-Boys optimization is somewhat random and thus
                # we use relatively loose testing criteria
                np.testing.assert_almost_equal(res_fb['orbe'][i+1],
                                               res_fb['orbe'][j+1],
                                               decimal=2)

if __name__ == '__main__':
    unittest.main()
