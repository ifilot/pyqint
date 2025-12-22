# -*- coding: utf-8 -*-

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt

from .pyqint_core import PyQInt
from .molecule import Molecule

# These are assumed to exist elsewhere in your module
SUBSPACE_START: int
SUBSPACE_LENGTH: int


Vec = npt.NDArray[np.float64]
Mat = npt.NDArray[np.float64]

# nucleus entry: (position, charge)
Nucleus = Tuple[Vec, int]

# couple of hardcoded variables for the DIIS algorithm
SUBSPACE_LENGTH = 4
SUBSPACE_START = 1

class HF:
    """
    Perform Hartree-Fock calculations
    """
    def __init__(self, mol: Molecule, basis: Union[str, Sequence[Any]]) -> None:
        """
        Initialize the Hartree-Fock solver.

        Parameters
        ----------
        mol
            Molecule object.
        basis
            Basis set name or a list of CGF objects.
        """
        self._mol: Molecule = mol
        self._basis: Union[str, Sequence[Any]] = basis

        # build cgfs, nuclei and calculate nr of electrons
        if issubclass(type(basis), str): # if a basis set name is given
            self._cgfs, self._nuclei = self._mol.build_basis(basis)
        else: # either assume a list of CGFs objects is given
            self._cgfs = basis
            self._nuclei = mol.get_nuclei()
            
        self._nelec = mol.get_nelec()

    def rhf(
        self,
        calc_forces: bool = False,
        itermax: int = 100,
        use_diis: bool = True,
        verbose: bool = False,
        tolerance: float = 1e-9,
        orbc_init: Optional[Mat] = None,
        ortho: str = "canonical",
    ) -> Dict[str, Any]:
        """
        Perform a restricted Hartree-Fock calculation.

        Parameters
        ----------
        calc_forces
            Whether analytic nuclear forces are computed.
        itermax
            Maximum number of SCF iterations.
        use_diis
            Enable DIIS acceleration.
        verbose
            Print per-iteration diagnostics.
        tolerance
            Energy convergence threshold.
        orbc_init
            Optional initial MO coefficient matrix.
        ortho
            Orthogonalization scheme ("canonical" or "symmetric").

        Returns
        -------
        dict
            Result dictionary containing energies, orbitals, matrices, etc.
        """

        # crank up the tolerance when calculating forces
        if calc_forces and tolerance > 1e-12:
            tolerance = 1e-12

        # create empty dictionary for time tracking statistics
        time_stats = {}
        N = len(self._cgfs)
        occ = [2 if i < self._nelec//2 else 0 for i in range(N)]

        # build integrals
        integrator = PyQInt()
        start = time.time()
        S, T, V, tetensor = integrator.build_integrals_openmp(self._cgfs, self._nuclei)
        end = time.time()
        time_stats['integral_evaluation'] = end - start

        # calculate nuclear repulsion
        nuc_rep = 0.0
        for i in range(0, len(self._nuclei)):
            for j in range(i+1, len(self._nuclei)):
                r = np.linalg.norm(np.array(self._nuclei[i][0]) - np.array(self._nuclei[j][0]))
                nuc_rep += self._nuclei[i][1] * self._nuclei[j][1] / r

        # diagonalize S
        s, U = np.linalg.eigh(S)

        # construct transformation matrix X
        if ortho == 'canonical': # perform canonical orthogonalization
            X = U @ np.diag(1.0/np.sqrt(s))
        elif ortho == 'symmetric':
            X = U @ np.diag(1.0/np.sqrt(s)) @ U.transpose()
        else:
            raise Exception("Invalid orthogonalization option selected: ", ortho)
        
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
                    diis_coeff = self.__calculate_diis_coefficients(evs_diis)
                except np.linalg.LinAlgError:
                    # stop diis procedure and revert to linear stepping
                    use_diis = False
                    continue

                F = self.__extrapolate_fock_from_diis_coefficients(fmats_diis, diis_coeff)
                Fprime = X.transpose() @ F @ X
                e, Cprime = np.linalg.eigh(Fprime)
                C = X @ Cprime
                P = np.einsum('ik,jk,k->ij', C, C, occ)

            # calculate G
            G = (np.einsum('kl,ijlk->ij', P, tetensor) - 0.5 * np.einsum('kl,iklj->ij', P, tetensor))

            # build Fock matrix
            F = T + V + G

            # transform Fock matrix
            Fprime = X.transpose() @ F @ X

            # diagonalize F
            orbe, Cprime = np.linalg.eigh(Fprime)

            # back-transform
            C = X @ Cprime

            # calculate energy E
            energy = 0.0
            M = T + V + F
            energy = 0.5 * np.einsum('ij,ji', P, M)

            # add nuclear repulsion
            energy += nuc_rep

            # store energy for next iteration
            energies.append(energy)

            # for the first few iterations, build a new density
            # matrix from the coefficients, else, resort to the DIIS
            # algorithm
            if niter <= SUBSPACE_START or not use_diis:
                P = np.einsum('ik,jk,k->ij', C, C, occ)

            # calculate DIIS coefficients
            e = (F.dot(P.dot(S)) - S.dot(P.dot(F))).flatten()   # calculate error vector
            #enorm = np.linalg.norm(e)                          # store error vector norm
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
                print("Iteration: %2i | Energy: %12.6f Ht | Time: %6.4f s" % (niter, energy, time_stats['iterations'][-1]))

            # calculate energy difference between this and the previous
            # iteration; terminate the loop when energy difference is less
            # than threshold
            if niter > 4: # require at least 4 steps
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

        # update density matrix and final energy
        P = np.einsum('ik,jk,k->ij', C, C, occ)
        energies[-1] = 0.5 * np.einsum('ji,ij', P, T+V+F) + nuc_rep

        # build solution dictionary
        sol = {
            "energy": energies[-1],
            "nuclei" : self._nuclei,
            "cgfs": self._cgfs,
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
            'tetensor': tetensor,
            "time_stats" : time_stats,
            "ecore": np.sum(P * (T + V)),
            "ekin": np.einsum('ij,ji', T, P),
            "enuc": np.einsum('ij,ji', V, P),
            "erep": 0.5 * np.einsum('ijlk,ij,kl', tetensor, P, P),
            "ex": -0.25 * np.einsum('iklj,ij,kl', tetensor, P, P),
            "enucrep": nuc_rep,
            "nelec": self._nelec,
            "mol": self._mol,
            "forces": self.__rhf_forces(self._mol, self._basis, C, P, orbe) if calc_forces else None
        }

        return sol

    def uhf(
            self,
            multiplicity: int,
            itermax: int = 100,
            use_diis: bool = True,
            verbose: bool = False,
            tolerance: float = 1e-9,
            orbc_init: Optional[Dict[str, "Mat"]] = None,   # expects {"alpha": Ca, "beta": Cb}
            ortho: str = "canonical",
        ) -> Dict[str, Any]:
        """
        Perform an unrestricted Hartree-Fock calculation.

        Parameters
        ----------
        multiplicity
            Spin multiplicity (number of unpaired electrons)
        calc_forces
            Whether analytic nuclear forces are computed. (Not implemented here.)
        itermax
            Maximum number of SCF iterations.
        use_diis
            Enable DIIS acceleration.
        verbose
            Print per-iteration diagnostics.
        tolerance
            Energy convergence threshold.
        orbc_init
            Optional initial MO coefficient matrices: {"alpha": Ca, "beta": Cb}.
        ortho
            Orthogonalization scheme ("canonical" or "symmetric").

        Returns
        -------
        dict
            Result dictionary containing energies, orbitals, matrices, etc.
        """

        # electron counts
        nelec = int(self._nelec)
        if multiplicity < 1 or multiplicity > nelec + 1:
            raise ValueError(f"Invalid multiplicity={multiplicity} for nelec={nelec}.")

        # Nα - Nβ = 2S = multiplicity-1 ; Nα + Nβ = nelec
        nalpha = (nelec + (multiplicity - 1)) // 2
        nbeta  = nelec - nalpha
        if nalpha < 0 or nbeta < 0:
            raise ValueError("Computed negative Nalpha/Nbeta; check nelec/multiplicity.")

        # basis size / occupations
        time_stats = {}
        N = len(self._cgfs)
        occ_a = np.array([1 if i < nalpha else 0 for i in range(N)], dtype=float)
        occ_b = np.array([1 if i < nbeta  else 0 for i in range(N)], dtype=float)

        # build integrals
        integrator = PyQInt()
        start = time.time()
        S, T, V, tetensor = integrator.build_integrals_openmp(self._cgfs, self._nuclei)
        end = time.time()
        time_stats['integral_evaluation'] = end - start

        # calculate nuclear repulsion
        nuc_rep = 0.0
        for i in range(0, len(self._nuclei)):
            for j in range(i+1, len(self._nuclei)):
                r = np.linalg.norm(np.array(self._nuclei[i][0]) - np.array(self._nuclei[j][0]))
                nuc_rep += self._nuclei[i][1] * self._nuclei[j][1] / r

        # diagonalize S
        s, U = np.linalg.eigh(S)

        # construct transformation matrix X
        if ortho == 'canonical':
            X = U @ np.diag(1.0/np.sqrt(s))
        elif ortho == 'symmetric':
            X = U @ np.diag(1.0/np.sqrt(s)) @ U.transpose()
        else:
            raise Exception("Invalid orthogonalization option selected: ", ortho)

        # initial guess densities
        if orbc_init is None:
            Pa = np.zeros_like(S)
            Pb = np.zeros_like(S)
        else:
            Ca0 = orbc_init.get("alpha", None)
            Cb0 = orbc_init.get("beta", None)
            if Ca0 is None or Cb0 is None:
                raise ValueError("orbc_init for UHF must be a dict with keys {'alpha','beta'}.")
            Pa = np.einsum('ik,jk,k->ij', Ca0, Ca0, occ_a)
            Pb = np.einsum('ik,jk,k->ij', Cb0, Cb0, occ_b)

        # time tracking
        start = time.time()

        # per-iteration containers
        energies = []
        time_stats['iterations'] = []

        # DIIS containers (separate DIIS for alpha and beta)
        fmats_diis_a, evs_diis_a = [], []
        fmats_diis_b, evs_diis_b = [], []

        # predefine outputs in case itermax==0
        orbe_a = orbe_b = None
        Ca = Cb = None
        Fa = Fb = None

        def build_J(Ptot):
            """
            Auxiliary matrix for coulombic repulsion
            """
            return np.einsum('kl,ijlk->ij', Ptot, tetensor)

        def build_K(Pspin):
            """
            Auxiliary matrix for exchange
            """
            return np.einsum('kl,iklj->ij', Pspin, tetensor)

        # SCF iterations
        for niter in range(0, itermax):
            iterstart = time.time()

            # Optionally do DIIS extrapolation once we have a subspace
            did_diis = False
            if niter > SUBSPACE_START and use_diis:
                try:
                    coeff_a = self.__calculate_diis_coefficients(evs_diis_a)
                    coeff_b = self.__calculate_diis_coefficients(evs_diis_b)
                except np.linalg.LinAlgError:
                    # stop diis procedure and revert to linear stepping
                    use_diis = False
                    coeff_a = coeff_b = None

                if use_diis:
                    Fa = self.__extrapolate_fock_from_diis_coefficients(fmats_diis_a, coeff_a)
                    Fb = self.__extrapolate_fock_from_diis_coefficients(fmats_diis_b, coeff_b)

                    # diagonalize both and update densities
                    Fa_p = X.transpose() @ Fa @ X
                    Fb_p = X.transpose() @ Fb @ X

                    orbe_a, Ca_p = np.linalg.eigh(Fa_p)
                    orbe_b, Cb_p = np.linalg.eigh(Fb_p)

                    Ca = X @ Ca_p
                    Cb = X @ Cb_p

                    Pa = np.einsum('ik,jk,k->ij', Ca, Ca, occ_a)
                    Pb = np.einsum('ik,jk,k->ij', Cb, Cb, occ_b)

                    did_diis = True

            # Build Coulomb/exchange from current densities
            Ptot = Pa + Pb
            J = build_J(Ptot)
            Ka = build_K(Pa)
            Kb = build_K(Pb)

            # Build Fock matrices
            Fa = T + V + J - Ka
            Fb = T + V + J - Kb

            # Transform and diagonalize
            Fa_p = X.transpose() @ Fa @ X
            Fb_p = X.transpose() @ Fb @ X

            orbe_a, Ca_p = np.linalg.eigh(Fa_p)
            orbe_b, Cb_p = np.linalg.eigh(Fb_p)

            Ca = X @ Ca_p
            Cb = X @ Cb_p

            # Build new densities
            Pa_new = np.einsum('ik,jk,k->ij', Ca, Ca, occ_a)
            Pb_new = np.einsum('ik,jk,k->ij', Cb, Cb, occ_b)

            # Calculate energy
            Ptot_new = Pa_new + Pb_new
            J_new = build_J(Ptot_new)
            Ka_new = build_K(Pa_new)
            Kb_new = build_K(Pb_new)

            E_elec = (
                np.einsum('ij,ji', Ptot_new, (T + V))
                + 0.5 * np.einsum('ij,ji', Ptot_new, J_new)
                - 0.5 * np.einsum('ij,ji', Pa_new, Ka_new)
                - 0.5 * np.einsum('ij,ji', Pb_new, Kb_new)
            )
            energy = E_elec + nuc_rep
            energies.append(energy)

            # update densities for next iteration
            Pa, Pb = Pa_new, Pb_new

            # DIIS error vectors for alpha and beta
            ea = (Fa.dot(Pa.dot(S)) - S.dot(Pa.dot(Fa))).flatten()
            eb = (Fb.dot(Pb.dot(S)) - S.dot(Pb.dot(Fb))).flatten()

            fmats_diis_a.append(Fa); evs_diis_a.append(ea)
            fmats_diis_b.append(Fb); evs_diis_b.append(eb)

            # prune DIIS subspace
            if len(fmats_diis_a) > SUBSPACE_LENGTH:
                fmats_diis_a = fmats_diis_a[-SUBSPACE_LENGTH:]
                evs_diis_a   = evs_diis_a[-SUBSPACE_LENGTH:]
            if len(fmats_diis_b) > SUBSPACE_LENGTH:
                fmats_diis_b = fmats_diis_b[-SUBSPACE_LENGTH:]
                evs_diis_b   = evs_diis_b[-SUBSPACE_LENGTH:]

            iterend = time.time()
            time_stats['iterations'].append(iterend - iterstart)

            if verbose:
                print(
                    "Iteration: %2i | Energy: %16.10f Ht | Time: %6.4f s"
                    % (niter, energy, time_stats['iterations'][-1])
                )

            # convergence check
            if niter > 4:
                ediff = abs(energies[-1] - energies[-2])
                if ediff < tolerance:
                    if verbose:
                        print("Stopping SCF cycle, convergence reached.")
                    break

        end = time.time()
        time_stats['self_convergence'] = end - start

        # Final recompute of common matrices for reporting
        Ptot = Pa + Pb
        J = build_J(Ptot)
        Ka = build_K(Pa)
        Kb = build_K(Pb)
        Fa = T + V + J - Ka
        Fb = T + V + J - Kb

        # calculate energy terms
        E_core = np.einsum('ij,ji', Ptot, (T + V))
        E_J    = 0.5 * np.einsum('ij,ji', Ptot, J)
        E_Ka   = -0.5 * np.einsum('ij,ji', Pa, Ka)
        E_Kb   = -0.5 * np.einsum('ij,ji', Pb, Kb)
        E_elec = E_core + E_J + E_Ka + E_Kb
        E_tot  = E_elec + nuc_rep

        sol = {
            "energy": E_tot,
            "nuclei": self._nuclei,
            "cgfs": self._cgfs,
            "energies": energies,
            "orbe_alpha": orbe_a,
            "orbe_beta": orbe_b,
            "orbc_alpha": Ca,
            "orbc_beta": Cb,
            "density_alpha": Pa,
            "density_beta": Pb,
            "density": Ptot,
            "fock_alpha": Fa,
            "fock_beta": Fb,
            "transform": X,
            "overlap": S,
            "kinetic": T,
            "nuclear": V,
            "hcore": T + V,
            "tetensor": tetensor,
            "time_stats": time_stats,

            # energy breakdown
            "ecore": E_core,
            "ej": E_J,
            "ex_alpha": E_Ka,
            "ex_beta": E_Kb,
            "enucrep": nuc_rep,

            "nelec": nelec,
            "nalpha": nalpha,
            "nbeta": nbeta,
            "multiplicity": multiplicity,

            "mol": self._mol,
            "forces": None, # Forces: not available for UHF
        }

        return sol

    def __calculate_diis_coefficients(self, evs_diis):
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

    def __extrapolate_fock_from_diis_coefficients(self, fmats_diis, diis_coeff):
        """
        Extrapolate the Fock matrix from the DIIS coefficients
        """
        norbs = fmats_diis[-1].shape[0]
        fguess = np.zeros((norbs,norbs))

        for i in range(len(fmats_diis)):
            fguess += fmats_diis[i]*diis_coeff[i]

        return fguess

    def __rhf_forces(self, mol, basis, C, P, e):
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
        nelec = mol.get_nelec()
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
