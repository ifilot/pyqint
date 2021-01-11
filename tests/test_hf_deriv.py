import unittest
from pyqint import PyQInt, cgf, gto, Molecule, HF
from copy import deepcopy
import numpy as np
import multiprocessing
import os
from nose.tools import nottest

class TestHFDeriv(unittest.TestCase):

    def testHartreeFockForces(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        # calculate forces using analytical derivatives
        solver = HF()
        res = solver.rhf(mol, 'sto3g', calc_forces=True)

        # calculate forces using finite difference
        forces = calculate_forces_finite_difference(mol)

        np.testing.assert_almost_equal(res['forces'], forces)

def perform_hf(mol):
    sol = HF().rhf(mol, 'sto3g')
    return sol['energy']

def calculate_forces_finite_difference(mol):
    forces = np.zeros((3,3))

    sz = 0.000001

    for i in range(0, len(mol.atoms)): # loop over nuclei
        for j in range(0, 3): # loop over directions
            mol1 = deepcopy(mol)
            mol1.atoms[i][1][j] -= sz / 2
            mol2 = deepcopy(mol)
            mol2.atoms[i][1][j] += sz / 2

            energy1 = perform_hf(mol1)
            energy2 = perform_hf(mol2)

            forces[i,j] = (energy2 - energy1) / sz

    return forces

if __name__ == '__main__':
    unittest.main()
