import unittest
from pyqint import PyQInt, cgf, gto, Molecule, HF
from copy import deepcopy
import numpy as np
import multiprocessing
import os

class TestFerrocene(unittest.TestCase):

    @unittest.skip("Skipping expensive test")
    def testHartreeFockFerrocene(self):
        """
        Test Hartree-Fock calculation on Ferrocene
        """

        # build Integrator object and check if it has more than 1 thread
        integrator = PyQInt()
        numthreads = integrator.get_num_threads()
        self.assertTrue(numthreads > 1)

        atoms = read_xyz(os.path.join(os.path.dirname(__file__), 'input', 'ferrocene.xyz'))
        mol = Molecule('Ferrocene')
        for atom in atoms:
            mol.add_atom(atom[0], atom[1], atom[2], atom[3], unit='angstrom')

        results = perform_hf(mol)
        print(results['energy'])

def perform_hf(mol):
    return HF().rhf(mol, 'sto3g')

def read_xyz(filename):
    f = open(filename, 'r')
    nratoms = int(f.readline())
    f.readline()

    # read all atoms
    atoms = []
    for i in range(0, nratoms):
        pieces = f.readline().split()
        atoms.append([pieces[0], float(pieces[1]), float(pieces[2]), float(pieces[3])])

    # center structure
    ctr = [np.average([at[1] for at in atoms]),
           np.average([at[2] for at in atoms]),
           np.average([at[3] for at in atoms])]

    for at in atoms:
        for i in range(0,3):
            at[i+1] -= ctr[i]

    return atoms

if __name__ == '__main__':
    unittest.main()
