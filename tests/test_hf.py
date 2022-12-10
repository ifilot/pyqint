import unittest
from pyqint import PyQInt, cgf, gto, Molecule, HF
from copy import deepcopy
import numpy as np
import multiprocessing
import os

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

        # verify that time statistics are being recorded
        self.assertTrue(results['time_stats']['integral_evaluation'] > 0)
        self.assertTrue(results['time_stats']['self_convergence'] > 0)

    def testHartreeFockCH4Properties(self):
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

        results = perform_hf(mol)

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

def perform_hf(mol):
    results = HF().rhf(mol, 'sto3g')
    return results

if __name__ == '__main__':
    unittest.main()
