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

        energy = perform_hf(mol)
        np.testing.assert_almost_equal(energy, -73.21447132, 4)

def perform_hf(mol):
    sol = HF().rhf(mol, 'sto3g')
    return sol['energy']

if __name__ == '__main__':
    unittest.main()
