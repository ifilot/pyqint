# -*- coding: utf-8 -*-

import json
import os
from .cgf import cgf
import numpy as np
from . import PyQInt
import time

# couple of hardcoded variables for the DIIS algorithm
SUBSPACE_LENGTH = 8
SUBSPACE_START = 1

class HF:
    """
    Routines to perform a restricted Hartree-Fock calculations
    """
    def rhf(self, mol, basis, calc_forces=False, itermax=100, verbose=False):
        """
        Performs a Hartree-Fock type calculation

        mol:            molecule
        basis:          basis set
        calc_forces:    whether first derivatives need to be calculed
        verbose:        whether verbose output is given
        """

        # create empty dictionary for time tracking statistics
        time_stats = {}

        # build cgfs, nuclei and calculate nr of electrons
        cgfs, nuclei = mol.build_basis(basis)
        nelec = int(np.sum([at[1] for at in nuclei]))

        # build integrals
        integrator = PyQInt()
        start = time.time()
        S, T, V, teint = integrator.build_integrals_openmp(cgfs, nuclei)
        end = time.time()
        time_stats['integral_evaluation'] = end - start

        # diagonalize S
        s, U = np.linalg.eigh(S)

        # construct transformation matrix X
        X = U.dot(np.diag(1.0/np.sqrt(s)))

        # create empty P matrix as initial guess
        P = np.zeros(S.shape)

        # keep track of time
        start = time.time()

        # build containers to store per-iteration data
        energies = []
        time_stats['iterations'] = []
        fmats_diis = []
        pmat_diis = []
        evs_diis = []

        # start iterations
        for niter in range(0,itermax):
            # keep track of iterations
            iterstart = time.time()

            if niter > SUBSPACE_START:
                diis_coeff = self.calculate_diis_coefficients(evs_diis)

                F = self.extrapolate_fock_from_diis_coefficients(fmats_diis, diis_coeff)
                Fprime = X.transpose().dot(F).dot(X)
                e, Cprime = np.linalg.eigh(Fprime)
                C = X.dot(Cprime)
                P = np.zeros(S.shape)
                for i in range(S.shape[0]):
                    for j in range(S.shape[0]):
                        for k in range(0,int(nelec/2)):
                            P[i,j] += 2.0 * C[i,k] * C[j,k]

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
                    # store iteration time
                    iterend = time.time()
                    time_stats['iterations'].append(iterend - iterstart)

                    # terminate self-convergence cycle
                    if verbose:
                        print("Stopping SCF cycle, convergence reached.")
                    break

            # store energy for next iteration
            energies.append(energy)

            # for the first few iterations, build a new density
            # matrix from the coefficients, else, resort to the DIIS
            # algorithm
            if niter <= SUBSPACE_START:
                P = np.zeros(S.shape)
                for i in range(S.shape[0]):
                    for j in range(S.shape[0]):
                        for k in range(0,int(nelec/2)):
                            P[i,j] += 2.0 * C[i,k] * C[j,k]

            # calculate DIIS coefficients
            e = (F.dot(P.dot(S)) - S.dot(P.dot(F))).flatten()   # calculate error vector
            enorm = np.linalg.norm(e)                           # store error vector norm
            fmats_diis.append(F)                                # add Fock matrix to list
            pmat_diis.append(P)                                 # add density matrix to list
            evs_diis.append(e)

            # prune size of the old Fock, density and error vector lists
            # only SUBSPACE_LENGTH iterations are used to guess the new
            # solution
            if len(fmats_diis) > SUBSPACE_LENGTH:
                fmats_diis = fmats_diis[-SUBSPACE_LENGTH:]
                pmat_diis = pmat_diis[-SUBSPACE_LENGTH:]
                evs_diis = evs_diis[-SUBSPACE_LENGTH:]

            # store iteration time
            iterend = time.time()
            time_stats['iterations'].append(iterend - iterstart)

        # store time for self-converging iterations
        end = time.time()
        time_stats['self_convergence'] = end - start

        # build solution dictionary
        sol = {
            "energy": energies[-1],
            "nuclei" : nuclei,
            "cgfs": cgfs,
            "energies": energies,
            "orbe": e,
            "orbc": C,
            "density": P,
            "overlap": S,
            "kinetic": T,
            "nuclear": V,
            "time_stats" : time_stats,
            "ecore": np.sum(P * (T + V)),
            "teint": teint,
            "forces": self.rhf_forces(mol, basis, C, P, e) if calc_forces else None
        }

        return sol

    def calculate_diis_coefficients(self, evs_diis):
        """
        Calculate the DIIS coefficients
        """
        B = np.zeros((len(evs_diis)+1, len(evs_diis)+1))
        B[-1,:] = -1
        B[:,-1] = -1
        B[-1,-1]=  0

        rhs = np.zeros((len(evs_diis)+1, 1))
        rhs[-1,-1] = -1

        for i in range(len(evs_diis)):
            for j in range(i+1):
                B[i,j] = np.dot(evs_diis[i].transpose(), evs_diis[j])
                B[j,i] = B[i,j]

        *diis_coeff, _ = np.linalg.solve(B,rhs)

        return diis_coeff

    def extrapolate_fock_from_diis_coefficients(self, fmats_diis, diis_coeff):
        """
        Extrapolate the Fock matrix from the DIIS coefficients
        """
        norbs = fmats_diis[-1].shape[0]
        fguess = np.zeros((norbs,norbs))

        for i in range(len(fmats_diis)):
            fguess += fmats_diis[i]*diis_coeff[i]

        return fguess

    # def rhf_forces(self, mol, basis, C, P, e):
    #     forces = np.zeros((len(mol.atoms), 3))

    #     for i in range(0, len(mol.atoms)): # loop over nuclei
    #         for j in range(0, 3): # loop over directions
    #             forces[i,j] = self.rhf_force_nuc_dir(mol, basis, C, P, e, i, j)

    #     return forces

    # def rhf_force_core(self, mol, basis, nucid, direction):
    #     # build cgfs, nuclei and calculate nr of electrons
    #     cgfs, nuclei = mol.build_basis(basis)

    #     # build integrator object
    #     integrator = PyQInt()

    #     # build overlap and kinetic derivatives
    #     T = np.zeros((len(cgfs), len(cgfs)))
    #     for i in range(0, len(cgfs)):
    #         for j in range(0, len(cgfs)):
    #             T[i,j] = integrator.kinetic_deriv(cgfs[i], cgfs[j], nuclei[nucid][0], direction)

    #     # build nuclear derivatives
    #     V = np.zeros((len(cgfs), len(cgfs)))
    #     for i in range(0, len(cgfs)):
    #         for j in range(0, len(cgfs)):
    #             for k in range(0, len(nuclei)):
    #                 V[i,j] += integrator.nuclear_deriv(cgfs[i], cgfs[j], nuclei[k][0], nuclei[k][1], nuclei[nucid][0], direction)

    #     return T + V

    # def rhf_force_nuc_dir(self, mol, basis, C, P, e, nucleus, direction):
    #     # build cgfs, nuclei and calculate nr of electrons
    #     cgfs, nuclei = mol.build_basis(basis)
    #     nelec = int(np.sum([at[1] for at in nuclei]))
    #     nratoms = len(nuclei)

    #     # build integrator object
    #     integrator = PyQInt()

    #     # build overlap and kinetic derivatives
    #     S = np.zeros((len(cgfs), len(cgfs)))
    #     for i in range(0, len(cgfs)):
    #         for j in range(i, len(cgfs)):
    #             S[i,j] = S[j,i] = integrator.overlap_deriv(cgfs[i], cgfs[j], nuclei[nucleus][0], direction)

    #     # build Q matrix
    #     Q = np.zeros(S.shape)
    #     for i in range(S.shape[0]):
    #         for j in range(S.shape[0]):
    #             for k in range(0,int(nelec/2)):
    #                 Q[i,j] += 2.0 * e[k] * C[i,k] * C[j,k]

    #     # build two-electron derivatives
    #     N = len(cgfs)
    #     teint_calc = np.multiply(np.ones(integrator.teindex(N,N,N,N)), -1.0)
    #     teint = np.zeros(integrator.teindex(N,N,N,N))
    #     for i, cgf1 in enumerate(cgfs):
    #         for j, cgf2 in enumerate(cgfs):
    #             ij = i*(i+1)/2 + j
    #             for k, cgf3 in enumerate(cgfs):
    #                 for l, cgf4 in enumerate(cgfs):
    #                     kl = k * (k+1)/2 + l
    #                     if ij <= kl:
    #                         idx = integrator.teindex(i,j,k,l)
    #                         if teint_calc[idx] < 0:
    #                             teint_calc[idx] = 1
    #                             teint[idx] = integrator.repulsion_deriv(cgfs[i], cgfs[j], cgfs[k], cgfs[l], nuclei[nucleus][0], direction)

    #     # build H-core derivatives
    #     Hcore = self.rhf_force_core(mol, basis, nucleus, direction)

    #     # calculate electronic derivate
    #     deriv = 0.0
    #     for i in range(0, len(cgfs)):
    #         for j in range(0, len(cgfs)):
    #             deriv += P[j,i] * Hcore[i,j]
    #             for k in range(0, len(cgfs)):
    #                 for l in range(0, len(cgfs)):
    #                     idx_rep = integrator.teindex(i,j,k,l)
    #                     idx_exc = integrator.teindex(i,l,k,j)
    #                     deriv += 0.5 * P[i,j] * P[k,l] * (teint[idx_rep] - 0.5 * teint[idx_exc])
    #             deriv -= Q[i,j] * S[i,j]

    #     # calculate nuclear derivative
    #     Vnn = 0.0
    #     pc = nuclei[nucleus][0]
    #     for i in range(0, len(nuclei)):
    #         if nucleus != i:
    #             pi = nuclei[i][0]
    #             Vnn += nuclei[nucleus][1] * nuclei[i][1] * (pi[direction] - pc[direction]) / np.linalg.norm(pi - pc)**3

    #     return deriv + Vnn

    # def check_symmetric(self, a, rtol=1e-05, atol=1e-08):
    #     return np.allclose(a, a.T, rtol=rtol, atol=atol)
