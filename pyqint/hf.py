# -*- coding: utf-8 -*-

import json
import os
import numpy as np
from . import PyQInt
import time

# couple of hardcoded variables for the DIIS algorithm
SUBSPACE_LENGTH = 4
SUBSPACE_START = 1

class HF:
    """
    Routines to perform a restricted Hartree-Fock calculations
    """
    def rhf(self, mol, basis, calc_forces=False, itermax=100,
            use_diis=True, verbose=False, tolerance=1e-7,
            orbc_init=None):
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
        N = len(cgfs)
        occ = [2 if i < nelec//2 else 0 for i in range(N)]

        # build integrals
        integrator = PyQInt()
        start = time.time()
        S, T, V, teints = integrator.build_integrals_openmp(cgfs, nuclei)
        end = time.time()
        time_stats['integral_evaluation'] = end - start

        # calculate nuclear repulsion
        nuc_rep = 0.0
        for i in range(0, len(nuclei)):
            for j in range(i+1, len(nuclei)):
                r = np.linalg.norm(np.array(nuclei[i][0]) - np.array(nuclei[j][0]))
                nuc_rep += nuclei[i][1] * nuclei[j][1] / r

        # diagonalize S
        s, U = np.linalg.eigh(S)

        # construct transformation matrix X
        X = U.dot(np.diag(1.0/np.sqrt(s)))

        # create empty P matrix as initial guess
        if orbc_init is None:
            P = np.zeros(S.shape)
        else:
            P = np.einsum('ik,jk,k->ij', orbc_init, orbc_init, occ)

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

            if niter > SUBSPACE_START and use_diis:
                try:
                    diis_coeff = self.calculate_diis_coefficients(evs_diis)
                except np.linalg.LinAlgError:
                    # stop diis procedure and revert to linear stepping
                    use_diis = False
                    continue

                F = self.extrapolate_fock_from_diis_coefficients(fmats_diis, diis_coeff)
                Fprime = X.transpose().dot(F).dot(X)
                e, Cprime = np.linalg.eigh(Fprime)
                C = X.dot(Cprime)
                P = np.einsum('ik,jk,k->ij', C, C, occ)

            # calculate G
            G = np.zeros_like(S)
            for i in range(N):
                for j in range(N):
                    for k in range(N):
                        for l in range(N):
                            idx_rep = integrator.teindex(i,j,l,k)
                            idx_exc = integrator.teindex(i,k,l,j)
                            G[i,j] += P[k,l] * (teints[idx_rep] - 0.5 * teints[idx_exc])

            # build Fock matrix
            F = T + V + G

            # transform Fock matrix
            Fprime = X.transpose().dot(F).dot(X)

            # diagonalize F
            orbe, Cprime = np.linalg.eigh(Fprime)

            # back-transform
            C = X.dot(Cprime)

            # calculate energy E
            energy = 0.0
            M = T + V + F
            for i in range(S.shape[0]):
                for j in range(S.shape[0]):
                    energy += 0.5 * P[j,i] * M[i,j]

            # add nuclear repulsion
            energy += nuc_rep

            # store energy for next iteration
            energies.append(energy)

            # for the first few iterations, build a new density
            # matrix from the coefficients, else, resort to the DIIS
            # algorithm
            if niter <= SUBSPACE_START or not use_diis:
                Pold = P.copy()
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

            # print info for this iteration
            if verbose:
                print("Iteration: %i Energy: %f Time: %f" % (niter, energy, time_stats['iterations'][-1]))

            # calculate energy difference between this and the previous
            # iteration; terminate the loop when energy difference is less
            # than threshold
            if niter > 1:
                ediff = np.abs(energies[-2] - energies[-1])

                if ediff < tolerance: # convergence criterion needs to be at least 1e-7!
                    # store iteration time
                    iterend = time.time()
                    time_stats['iterations'].append(iterend - iterstart)

                    # terminate self-convergence cycle
                    if verbose:
                        print("Stopping SCF cycle, convergence reached.")
                    break

        # store time for self-converging iterations
        end = time.time()
        time_stats['self_convergence'] = end - start

        # build two-electron integral tensor
        tetens = np.zeros((N,N,N,N))
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        idx = integrator.teindex(i,j,k,l)
                        tetens[i,j,k,l] = teints[idx]

        # build solution dictionary
        sol = {
            "energy": energies[-1],
            "nuclei" : nuclei,
            "cgfs": cgfs,
            "energies": energies,
            "orbe": orbe,
            "orbc": C,
            "density": P,
            "fock": F,
            "transform": X,
            "overlap": S,
            "kinetic": T,
            "nuclear": V,
            'hcore': T+V,
            'tetens': tetens,
            "time_stats" : time_stats,
            "ecore": np.sum(P * (T + V)),
            "ekin": np.einsum('ij,ji', T, P),
            "enuc": np.einsum('ij,ji', V, P),
            "erep": np.einsum('ijlk,ij,kl', tetens, P, P),
            "ex": np.einsum('ijlk,ij,kl', tetens, P, P),
            "enucrep": nuc_rep,
            "teint": teints,
            "nelec": nelec,
            "forces": self.rhf_forces(mol, basis, C, P, orbe) if calc_forces else None
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

    def rhf_forces(self, mol, basis, C, P, e):
        """
        Calculate derivatives of nuclear coordinates to obtain forces.
        
        Parameters:
        mol:    Molecule object to be evaluated
        basis:  string identifier of the basis
        C:      matrix of orbital coefficients
        P:      density matrix
        e:      orbital energies

        Returns array of forces with shape [Natoms, 3]
        """
        # intialization
        integrator = PyQInt()
        cgfs, nuclei = mol.build_basis(basis)
        nelec = np.sum([nucleus[1] for nucleus in nuclei])
        forces = np.zeros((len(nuclei),3))
        N = len(cgfs)
        occ = [2 if i < nelec//2 else 0 for i in range(N)]
    
        # calculate energy weighted density matrix
        ew_density = np.einsum('ik,jk,k,k->ij', C, C, occ, e)

        # collect derivatives
        S, T, V, teints = integrator.build_geometric_derivatives_openmp(cgfs, nuclei)

        # Loop over all cartesian direction for every nucleus
        # This could be made more efficient by incorporating symmetry
        for n, deriv_nucleus in enumerate(nuclei):
            for d in range(3):

                # derivate nucleus-nucleus repulsion
                term_nn = 0
                for nucleus in nuclei:
                    if np.linalg.norm(nucleus[0] - deriv_nucleus[0]) > 0.0001:
                        term_nn += deriv_nucleus[1] * nucleus[1] * \
                                   (nucleus[0][d] - deriv_nucleus[0][d]) / \
                                   (np.linalg.norm(nucleus[0]- deriv_nucleus[0])**3)

                hcore = T[n,d] + V[n,d]
                term_hcore = np.sum(np.multiply(P, hcore))

                term_repulsion = 0
                for i in range(N):
                    for j in range(N):
                        for k in range(N):
                            for l in range(N):
                                idx1 = integrator.teindex(i,j,k,l)
                                idx2 = integrator.teindex(i,k,j,l)

                                term_repulsion += 0.5 * P[i,j] * P[k,l] * (teints[n,d,idx1] - 0.5 * teints[n,d,idx2])
                
                term_overlap = - np.sum(np.multiply(ew_density, S[n,d]))
           
                # F = - d/dR E
                forces[n,d] = term_hcore + term_repulsion + term_overlap + term_nn

        return forces
