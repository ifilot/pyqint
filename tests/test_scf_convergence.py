import unittest
import copy
import warnings
from pyqint import Molecule, HF
import numpy as np

def water():
    mol = Molecule()
    mol.add_atom('O', 0.0, 0.0, 0.0)
    mol.add_atom('H', 0.7570, 0.5860, 0.0)
    mol.add_atom('H', -0.7570, 0.5860, 0.0)
    return mol

class TestSCFConvergence(unittest.TestCase):

    def test_converged(self):
        """
        Test that a converged RHF calculation is reported as such and that its
        Fock and density matrices commute (FPS = SPF)
        """
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            res = HF(water(), 'sto3g').rhf()

        self.assertTrue(res['converged'])
        self.assertEqual(res['niter'], len(res['energies']))
        F, P, S = res['fock'], res['density'], res['overlap']
        self.assertLess(np.max(np.abs(F @ P @ S - S @ P @ F)), 1e-5)

        # the final energy belongs to the final density and Fock matrix
        H = res['hcore']
        np.testing.assert_almost_equal(res['energy'],
                                       0.5 * np.einsum('ij,ji', P, H + F) + res['enucrep'], 12)

    def test_not_converged(self):
        """
        Test that a warning is issued when the SCF does not converge
        """
        with self.assertWarns(RuntimeWarning):
            res = HF(water(), 'sto3g').rhf(itermax=3)
        self.assertFalse(res['converged'])

        with self.assertWarns(RuntimeWarning):
            res = HF(water(), 'sto3g').uhf(multiplicity=1, itermax=3)
        self.assertFalse(res['converged'])

    def test_uhf_converged(self):
        """
        Test that a converged UHF calculation is reported as such
        """
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            res = HF(water(), 'sto3g').uhf(multiplicity=1)

        self.assertTrue(res['converged'])
        ref = HF(water(), 'sto3g').rhf()
        np.testing.assert_almost_equal(res['energy'], ref['energy'], 8)

    def test_linear_dependency(self):
        """
        Test that canonical orthogonalization removes linear dependencies:
        duplicating every basis function leaves the energy unchanged
        """
        mol = Molecule()
        mol.add_atom('H', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        ref = HF(mol, cgfs).rhf()
        res = HF(mol, cgfs + [copy.deepcopy(c) for c in cgfs]).rhf()

        self.assertTrue(res['converged'])
        self.assertEqual(res['orbc'].shape, (4, 2))
        np.testing.assert_almost_equal(res['energy'], ref['energy'], 10)

    def test_odd_electrons(self):
        """
        Test that RHF rejects an odd number of electrons
        """
        mol = Molecule()
        mol.add_atom('H', 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            HF(mol, 'sto3g').rhf()

        res = HF(mol, 'sto3g').uhf(multiplicity=2)
        self.assertTrue(res['converged'])

if __name__ == '__main__':
    unittest.main()
