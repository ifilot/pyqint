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

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -73.21447132, 4)

        # verify that terms are being calculated
        np.testing.assert_almost_equal(results['density'], np.einsum('ik,jk,k->ij', results['orbc'], results['orbc'], [2,2,2,2,2,0,0]), decimal=5)
        np.testing.assert_almost_equal(results['ekin'] + results['enuc'] + results['erep'] + results['ex'] + results['enucrep'], results['energy'], decimal=5)

    def test_hartree_fock_ch4(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        dist = 1.78/2
        mol.add_atom('C', 0.0, 0.0, 0.0, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that orbital energies are correctly approximated
        ans = np.array([-11.0707,
                         -0.7392,
                         -0.3752,
                         -0.3752,
                         -0.3752,
                          0.2865,
                          0.4092,
                          0.4092,
                          0.4092])
        np.testing.assert_almost_equal(results['orbe'], ans, 4)

        en = -39.35007843284954
        np.testing.assert_almost_equal(results['energies'][-1], en, 4)
        
    def test_hartree_fock_ch4_symmetric(self):
        """
        Test Hartree-Fock calculation on CH4 using an STO-3g basis set and
        symmetric orthogonalization
        """
        mol = Molecule()
        dist = 1.78/2
        mol.add_atom('C', 0.0, 0.0, 0.0, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

        results = HF().rhf(mol, 'sto3g', ortho='symmetric')

        # check that orbital energies are correctly approximated
        ans = np.array([-11.0707,
                         -0.7392,
                         -0.3752,
                         -0.3752,
                         -0.3752,
                          0.2865,
                          0.4092,
                          0.4092,
                          0.4092])
        np.testing.assert_almost_equal(results['orbe'], ans, 4)

        en = -39.35007843284954
        np.testing.assert_almost_equal(results['energies'][-1], en, 4)

    def test_hartree_fock_restart(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        dist = 1.78/2
        mol.add_atom('C', 0.0, 0.0, 0.0, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')
        results1 = HF().rhf(mol, 'sto3g')

        # check that orbital energies are correctly approximated
        ans = np.array([-11.0707,
                         -0.7392,
                         -0.3752,
                         -0.3752,
                         -0.3752,
                          0.2865,
                          0.4092,
                          0.4092,
                          0.4092])
        np.testing.assert_almost_equal(results1['orbe'], ans, 4)

        en = -39.35007843284954
        np.testing.assert_almost_equal(results1['energies'][-1], en, 4)

        # create new CH4 molecule with slight adjustment in geometry and
        # seed the calculation with the previous converged result
        mol = Molecule()
        dist = 1.78/2
        mol.add_atom('C', 0.0, 0.0, 0.1, unit='angstrom')
        mol.add_atom('H', dist, dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
        mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
        mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

        # perform HF calculation with coefficient matrix from previous
        # calculation to speed up the convergence
        results2 = HF().rhf(mol, 'sto3g', orbc_init=results1['orbc'])

        # assess that the energy of the perturbed result is different
        # (and also higher)
        en = -39.34538546003782
        np.testing.assert_almost_equal(results2['energies'][-1], en, 3)

        # check that the convergence is quicker
        self.assertTrue(len(results1['energies']) > len(results2['energies']))

if __name__ == '__main__':
    unittest.main()
