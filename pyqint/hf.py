# -*- coding: utf-8 -*-

import json
import os
from .cgf import cgf
import numpy as np
from . import PyQInt

class HF:
    """
    Routines to perform Hartree-Fock calculations
    """

    def rhf(mol, basis, verbose=False):
        """
        Performs a Hartree-Fock type calculation
        """

        # build cgfs, nuclei and calculate nr of electrons
        cgfs, nuclei = mol.build_basis(basis)
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
        energies = []
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
            if verbose:
                print("Iteration: %i Energy: %f" % (niter, energy))

            # calculate energy difference between this and the previous
            # iteration; terminate the loop when energy difference is less
            # than threshold
            if niter > 1:
                ediff = np.abs(energy - energies[-1])
                if ediff < 1e-5:
                    print("Stopping SCF cycle, convergence reached.")
                    break

            # store energy for next iteration
            energies.append(energy)

            # calculate a new P
            P = np.zeros(S.shape)
            for i in range(S.shape[0]):
                for j in range(S.shape[0]):
                    for k in range(0,int(nelec/2)):
                        P[i,j] += 2.0 * C[i,k] * C[j,k]

        # build solution dictionary
        sol = {
            "energy": energies[-1],
            "energies": energies,
            "coefficients": C,
            "overlap": S,
            "kinetic": T,
            "nuclear": V,
            "teint": teint
        }

        return sol
