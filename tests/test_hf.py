import unittest
from pyqint import PyQInt, cgf, gto, Molecule
from copy import deepcopy
import numpy as np
import multiprocessing
import os

class TestCGF(unittest.TestCase):

    def testHartreeFockWater(self):
        """
        Test Hartree-Fock calculation on water
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        energy = perform_hf(mol)
        np.testing.assert_almost_equal(energy, -74.099410, 4)

def perform_hf(mol):
    # build cgfs, nuclei and calculate nr of electrons
    cgfs, nuclei = mol.build_basis('p321')
    nelec = int(np.sum([at[1] for at in nuclei]))

    # build integrals
    integrator = PyQInt()
    S, T, V, teint = integrator.build_integrals(cgfs, nuclei)

    # diagonalize S
    s, U = np.linalg.eigh(S)

    # construct transformation matrix X
    X = U.dot(np.diag(1.0/np.sqrt(s)))

    # create empty P matrix as initial guess
    P = np.zeros(S.shape)

    # start iterative procedure
    energy_old = 1.0
    for niter in range(0,100):

        # calculate G
        G = np.zeros(S.shape)
        for i in range(S.shape[0]):
            for j in range(S.shape[0]):
                for k in range(S.shape[0]):
                    for l in range(S.shape[0]):
                        idx_rep = integrator.teindex(i,j,l,k)
                        idx_exc = integrator.teindex(i,k,l,j)
                        G[i,j] += P[k,l] * (teint[idx_rep] - 0.5 * teint[idx_exc])

        # build Fock matrix
        F = T + V + G

        # transform Fock matrix
        Fprime = X.transpose().dot(F).dot(X)

        # diagonalize F
        e, Cprime = np.linalg.eigh(Fprime)

        # back-transform
        C = X.dot(Cprime)

        # calculate energy E
        energy = 0.0
        M = T + V + F
        for i in range(S.shape[0]):
            for j in range(S.shape[0]):
                energy += 0.5 * P[j,i] * M[i,j]

        # calculate repulsion of the nuclei
        for i in range(0, len(nuclei)):
            for j in range(i+1, len(nuclei)):
                r = np.linalg.norm(np.array(nuclei[i][0]) - np.array(nuclei[j][0]))
                energy += nuclei[i][1] * nuclei[j][1] / r

        # print info for this iteration
        # print("Iteration: %i Energy: %f" % (niter, energy))

        # calculate energy difference between this and the previous
        # iteration; terminate the loop when energy difference is less
        # than threshold
        ediff = np.abs(energy - energy_old)
        if ediff < 1e-5:
            print("Stopping SCF cycle, convergence reached.")
            break

        # store energy for next iteration
        energy_old = energy

        # calculate a new P
        P = np.zeros(S.shape)
        for i in range(S.shape[0]):
            for j in range(S.shape[0]):
                for k in range(0,int(nelec/2)):
                    P[i,j] += 2.0 * C[i,k] * C[j,k]

    return energy_old

if __name__ == '__main__':
    unittest.main()
