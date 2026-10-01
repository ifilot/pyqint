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

        res = HF(mol, 'sto3g').rhf()

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

        res = HF(mol, 'sto3g').rhf()

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

    def testCOUHF(self):
        """
        For a closed-shell system, localization of the UHF orbitals should
        yield the same result for both spin channels as for RHF
        """
        d = 1.145414
        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, -d/2, unit='angstrom')
        mol.add_atom('O', 0.0, 0.0,  d/2, unit='angstrom')

        res_rhf = HF(mol, 'sto3g').rhf(tolerance=1e-12)
        res_uhf = HF(mol, 'sto3g').uhf(multiplicity=1, tolerance=1e-12)

        res_fb_rhf = FosterBoys(res_rhf, seed=0).run(nr_runners=5)
        res_fb_uhf = FosterBoys(res_uhf, seed=0).run(nr_runners=5)

        for spin in ('alpha', 'beta'):
            np.testing.assert_almost_equal(res_fb_uhf['r2final_' + spin],
                                           res_fb_rhf['r2final'],
                                           decimal=3)
            np.testing.assert_almost_equal(res_fb_uhf['orbe_' + spin],
                                           res_fb_rhf['orbe'],
                                           decimal=2)

    def testCH3UHF(self):
        """
        Test localization of the alpha and beta orbitals of the methyl
        radical
        """
        R = 2.039
        sqrt3 = np.sqrt(3.0)

        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, 0.0)
        mol.add_atom('H',  R, 0.0, 0.0)
        mol.add_atom('H', -0.5 * R,  0.5 * sqrt3 * R, 0.0)
        mol.add_atom('H', -0.5 * R, -0.5 * sqrt3 * R, 0.0)

        res = HF(mol, 'sto3g').uhf(multiplicity=2, tolerance=1e-12)
        res_fb = FosterBoys(res, seed=0).run(nr_runners=3)

        for spin in ('alpha', 'beta'):
            nocc = res['n' + spin]
            C = res_fb['orbc_' + spin]

            # localization increases the Boys functional
            self.assertGreater(res_fb['r2final_' + spin],
                               res_fb['r2start_' + spin])

            # orbitals remain orthonormal
            np.testing.assert_almost_equal(C.T @ res['overlap'] @ C,
                                           np.identity(C.shape[1]))

            # rotations among occupied orbitals preserve the density
            Cocc = C[:, :nocc]
            np.testing.assert_almost_equal(Cocc @ Cocc.T,
                                           res['density_' + spin])

            # the three C-H bonds are equivalent
            np.testing.assert_almost_equal(res_fb['orbe_' + spin][nocc-3:nocc],
                                           [res_fb['orbe_' + spin][nocc-1]] * 3,
                                           decimal=3)

    def testSingleOccupiedUHF(self):
        """
        Spin channels with fewer than two occupied orbitals are left
        untouched
        """
        for element, nalpha, nbeta in (('H', 1, 0), ('Li', 2, 1)):
            mol = Molecule()
            mol.add_atom(element, 0.0, 0.0, 0.0)

            res = HF(mol, 'sto3g').uhf(multiplicity=2)
            res_fb = FosterBoys(res, seed=0).run()

            self.assertEqual(res_fb['nalpha'], nalpha)
            self.assertEqual(res_fb['nbeta'], nbeta)
            self.assertEqual(res_fb['nriter_beta'], 0)
            np.testing.assert_almost_equal(res_fb['orbe_beta'],
                                           res['orbe_beta'], decimal=5)

if __name__ == '__main__':
    unittest.main()
