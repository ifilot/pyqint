# -*- coding: utf-8 -*-

from __future__ import annotations

import time
import warnings
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt

from ..pyqint_core import PyQInt
from ..structure.molecule import Molecule

Vec = npt.NDArray[np.float64]
Mat = npt.NDArray[np.float64]

# nucleus entry: (position, charge)
Nucleus = Tuple[Vec, int]

# couple of hardcoded variables for the DIIS algorithm
SUBSPACE_LENGTH = 4
SUBSPACE_START = 1

# overlap eigenvalues below this threshold are removed by canonical
# orthogonalization (near-linear dependencies in the basis set)
LINDEP_THRESHOLD = 1e-7

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
        tolerance_grad: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Perform a restricted Hartree-Fock calculation.

        The calculation is converged when the change in energy is below
        ``tolerance`` and the largest element of the commutator FPS - SPF,
        expressed in the orthogonalized basis, is below ``tolerance_grad``.
        A RuntimeWarning is issued when this is not reached within
        ``itermax`` iterations; the ``converged`` entry of the result
        reports the outcome.

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
            Orthogonalization scheme ("canonical" or "symmetric"). Canonical
            orthogonalization removes the combinations of basis functions
            whose overlap eigenvalue is below 1e-7, such that the number of
            molecular orbitals can be smaller than the number of basis
            functions.
        tolerance_grad
            Threshold for the commutator FPS - SPF; defaults to
            sqrt(tolerance).

        Returns
        -------
        dict
            Result dictionary containing energies, orbitals, matrices, etc.

        Raises
        ------
        ValueError
            If the number of electrons is odd or exceeds the capacity of the
            molecular orbitals.
        """
        if self._nelec % 2 != 0:
            raise ValueError(
                f"RHF requires an even number of electrons, not {self._nelec}; "
                "use uhf() for open-shell systems."
            )
        if calc_forces and itermax < 1:
            raise ValueError("Calculating forces requires at least one SCF iteration.")

        # crank up the tolerance when calculating forces
        if calc_forces and tolerance > 1e-12:
            tolerance = 1e-12
        if tolerance_grad is None:
            tolerance_grad = np.sqrt(tolerance)

        # create empty dictionary for time tracking statistics
        time_stats = {}

        # build integrals
        integrator = PyQInt()
        start = time.time()
        S, T, V, tetensor = integrator.build_integrals_openmp(self._cgfs, self._nuclei)
        end = time.time()
        time_stats['integral_evaluation'] = end - start

        # calculate nuclear repulsion
        nuc_rep = self.__nuclear_repulsion()

        # construct transformation matrix X; the number of molecular
        # orbitals M is the number of columns of X
        X = self.__transformation_matrix(S, ortho)
        M = X.shape[1]
        if self._nelec // 2 > M:
            raise ValueError(
                f"Cannot place {self._nelec} electrons in {M} molecular orbitals "
                f"({len(self._cgfs)} basis functions)."
            )
        occ = [2 if i < self._nelec//2 else 0 for i in range(M)]

        # create empty P matrix as initial guess
        if orbc_init is None:
            P = np.zeros(S.shape)
        else:
            occ_init = [2 if i < self._nelec//2 else 0 for i in range(orbc_init.shape[1])]
            P = np.einsum('ik,jk,k->ij', orbc_init, orbc_init, occ_init)

        # keep track of time
        start = time.time()

        # build containers to store per-iteration data
        energies = []
        time_stats['iterations'] = []
        fmats_diis = []
        evs_diis = []
        converged = False

        # predefine outputs in case itermax==0
        orbe = C = None

        # start iterations
        for niter in range(0,itermax):
            # keep track of iterations
            iterstart = time.time()

            if niter > SUBSPACE_START and use_diis:
                try:
                    diis_coeff = self.__calculate_diis_coefficients(evs_diis)
                except np.linalg.LinAlgError:
                    # stop diis procedure and revert to linear stepping,
                    # i.e. continue from the orbitals of the previous
                    # iteration
                    use_diis = False
                    P = np.einsum('ik,jk,k->ij', C, C, occ)
                else:
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
            energy = 0.5 * np.einsum('ij,ji', P, T + V + F)

            # add nuclear repulsion
            energy += nuc_rep

            # store energy for next iteration
            energies.append(energy)

            # calculate the error vector from the density matrix that was
            # used to build the Fock matrix; it vanishes at self-consistency
            err = F.dot(P.dot(S)) - S.dot(P.dot(F))
            errmax = np.max(np.abs(X.transpose() @ err @ X))

            # store Fock matrix and error vector for DIIS; the zero initial
            # guess is skipped as its error vector vanishes trivially
            if niter > 0 or orbc_init is not None:
                fmats_diis.append(F)                            # add Fock matrix to list
                evs_diis.append(err.flatten())                  # add error vector to list

            # for the first few iterations, build a new density
            # matrix from the coefficients, else, resort to the DIIS
            # algorithm
            if niter <= SUBSPACE_START or not use_diis:
                P = np.einsum('ik,jk,k->ij', C, C, occ)

            # prune size of the old Fock, density and error vector lists
            # only SUBSPACE_LENGTH iterations are used to guess the new
            # solution
            if len(fmats_diis) > SUBSPACE_LENGTH:
                fmats_diis = fmats_diis[-SUBSPACE_LENGTH:]
                evs_diis = evs_diis[-SUBSPACE_LENGTH:]

            # store iteration time
            iterend = time.time()
            time_stats['iterations'].append(iterend - iterstart)

            # print info for this iteration
            if verbose:
                print("Iteration: %2i | Energy: %12.6f Ht | Error: %8.2e | Time: %6.4f s" % (niter, energy, errmax, time_stats['iterations'][-1]))

            # calculate energy difference between this and the previous
            # iteration; terminate the loop when both the energy difference
            # and the error vector are less than their thresholds
            if niter > 4: # require at least 4 steps
                ediff = np.abs(energies[-2] - energies[-1])

                if ediff < tolerance and errmax < tolerance_grad:
                    converged = True

                    # terminate self-convergence cycle
                    if verbose:
                        print("Stopping SCF cycle, convergence reached.")
                    break

        # store time for self-converging iterations
        end = time.time()
        time_stats['self_convergence'] = end - start

        if itermax > 0 and not converged:
            warnings.warn(
                f"RHF did not converge within {itermax} iterations "
                f"(last energy change {ediff if niter > 4 else float('nan'):.2e} Ht, "
                f"error {errmax:.2e}).",
                RuntimeWarning, stacklevel=2,
            )

        # build the density matrix from the final orbitals and the Fock
        # matrix and energy that belong to it; without any iterations, the
        # initial guess density is retained
        if C is not None:
            P = np.einsum('ik,jk,k->ij', C, C, occ)
        G = (np.einsum('kl,ijlk->ij', P, tetensor) - 0.5 * np.einsum('kl,iklj->ij', P, tetensor))
        F = T + V + G
        energy = 0.5 * np.einsum('ji,ij', P, T+V+F) + nuc_rep
        if energies:
            energies[-1] = energy

        # build solution dictionary
        sol = {
            "energy": energy,
            "converged": converged,
            "niter": len(energies),
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
            "forces": self.__rhf_forces(C, P, orbe) if calc_forces else None
        }

        return sol

    def uhf(
            self,
            multiplicity: int,
            itermax: int = 100,
            use_diis: bool = True,
            verbose: bool = False,
            tolerance: float = 1e-9,
            orbc_init: Optional[Dict[str, Mat]] = None,   # expects {"alpha": Ca, "beta": Cb}
            ortho: str = "canonical",
            nelec: Optional[int] = None,
            tolerance_grad: Optional[float] = None,
        ) -> Dict[str, Any]:
        """
        Perform an unrestricted Hartree-Fock calculation.

        Convergence is assessed as for :meth:`rhf`, using the larger of the
        commutators of the alpha and beta Fock and density matrices.

        Parameters
        ----------
        multiplicity
            Spin multiplicity 2S + 1 (number of unpaired electrons plus one)
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
            Orthogonalization scheme ("canonical" or "symmetric"); see
            :meth:`rhf`.
        nelec
            Number of electrons; defaults to the number of electrons of the
            neutral molecule.
        tolerance_grad
            Threshold for the commutators FPS - SPF; defaults to
            sqrt(tolerance).

        Returns
        -------
        dict
            Result dictionary containing energies, orbitals, matrices, etc.

        Raises
        ------
        ValueError
            If the multiplicity is incompatible with the number of electrons
            or the basis set is too small to hold the alpha electrons.
        """

        # electron counts
        if nelec is None:
            nelec = int(self._nelec)
        if multiplicity < 1 or multiplicity > nelec + 1:
            raise ValueError(f"Invalid multiplicity={multiplicity} for nelec={nelec}.")
        if (nelec + multiplicity - 1) % 2 != 0:
            raise ValueError(
                f"Multiplicity={multiplicity} is incompatible with nelec={nelec}: "
                "an even number of electrons requires an odd multiplicity and "
                "vice versa."
            )

        # Nα - Nβ = 2S = multiplicity-1 ; Nα + Nβ = nelec
        nalpha = (nelec + (multiplicity - 1)) // 2
        nbeta  = nelec - nalpha

        if tolerance_grad is None:
            tolerance_grad = np.sqrt(tolerance)

        time_stats = {}

        # build integrals
        integrator = PyQInt()
        start = time.time()
        S, T, V, tetensor = integrator.build_integrals_openmp(self._cgfs, self._nuclei)
        end = time.time()
        time_stats['integral_evaluation'] = end - start

        # calculate nuclear repulsion
        nuc_rep = self.__nuclear_repulsion()

        # construct transformation matrix X; the number of molecular
        # orbitals M is the number of columns of X
        X = self.__transformation_matrix(S, ortho)
        M = X.shape[1]
        if nalpha > M:
            raise ValueError(
                f"Cannot place {nalpha} alpha electrons in {M} molecular orbitals "
                f"({len(self._cgfs)} basis functions)."
            )

        # occupations
        occ_a = np.array([1 if i < nalpha else 0 for i in range(M)], dtype=float)
        occ_b = np.array([1 if i < nbeta  else 0 for i in range(M)], dtype=float)

        # initial guess densities
        if orbc_init is None:
            Pa = np.zeros_like(S)
            Pb = np.zeros_like(S)
        else:
            Ca0 = orbc_init.get("alpha", None)
            Cb0 = orbc_init.get("beta", None)
            if Ca0 is None or Cb0 is None:
                raise ValueError("orbc_init for UHF must be a dict with keys {'alpha','beta'}.")
            Pa = np.einsum('ik,jk,k->ij', Ca0, Ca0, (np.arange(Ca0.shape[1]) < nalpha).astype(float))
            Pb = np.einsum('ik,jk,k->ij', Cb0, Cb0, (np.arange(Cb0.shape[1]) < nbeta).astype(float))

        # time tracking
        start = time.time()

        # per-iteration containers
        energies = []
        time_stats['iterations'] = []

        # DIIS containers (separate DIIS for alpha and beta)
        fmats_diis_a, evs_diis_a = [], []
        fmats_diis_b, evs_diis_b = [], []
        converged = False

        # predefine outputs in case itermax==0
        orbe_a = orbe_b = None
        Ca = Cb = None
        Fa = Fb = None

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
            J = np.einsum('kl,ijlk->ij', Ptot, tetensor)
            Ka = np.einsum('kl,iklj->ij', Pa, tetensor)
            Kb = np.einsum('kl,iklj->ij', Pb, tetensor)

            # Build Fock matrices
            Fa = T + V + J - Ka
            Fb = T + V + J - Kb

            # DIIS error vectors for alpha and beta, from the densities that
            # were used to build the Fock matrices
            err_a = Fa.dot(Pa.dot(S)) - S.dot(Pa.dot(Fa))
            err_b = Fb.dot(Pb.dot(S)) - S.dot(Pb.dot(Fb))
            errmax = max(np.max(np.abs(X.transpose() @ err_a @ X)),
                         np.max(np.abs(X.transpose() @ err_b @ X)))

            # store them for DIIS; the zero initial guess is skipped as its
            # error vectors vanish trivially
            if niter > 0 or orbc_init is not None:
                fmats_diis_a.append(Fa); evs_diis_a.append(err_a.flatten())
                fmats_diis_b.append(Fb); evs_diis_b.append(err_b.flatten())

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
            J_new = np.einsum('kl,ijlk->ij', Ptot_new, tetensor)
            Ka_new = np.einsum('kl,iklj->ij', Pa_new, tetensor)
            Kb_new = np.einsum('kl,iklj->ij', Pb_new, tetensor)

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
                    "Iteration: %2i | Energy: %16.10f Ht | Error: %8.2e | Time: %6.4f s"
                    % (niter, energy, errmax, time_stats['iterations'][-1])
                )

            # convergence check
            if niter > 4:
                ediff = abs(energies[-1] - energies[-2])
                if ediff < tolerance and errmax < tolerance_grad:
                    converged = True
                    if verbose:
                        print("Stopping SCF cycle, convergence reached.")
                    break

        end = time.time()
        time_stats['self_convergence'] = end - start

        if itermax > 0 and not converged:
            warnings.warn(
                f"UHF did not converge within {itermax} iterations "
                f"(last energy change {ediff if niter > 4 else float('nan'):.2e} Ht, "
                f"error {errmax:.2e}).",
                RuntimeWarning, stacklevel=2,
            )

        # Final recompute of common matrices for reporting
        Ptot = Pa + Pb
        J = np.einsum('kl,ijlk->ij', Ptot, tetensor)
        Ka = np.einsum('kl,iklj->ij', Pa, tetensor)
        Kb = np.einsum('kl,iklj->ij', Pb, tetensor)
        Fa = T + V + J - Ka
        Fb = T + V + J - Kb

        # calculate energy terms
        E_core = np.einsum('ij,ji', Ptot, (T + V))
        E_J    = 0.5 * np.einsum('ij,ji', Ptot, J)
        E_Ka   = -0.5 * np.einsum('ij,ji', Pa, Ka)
        E_Kb   = -0.5 * np.einsum('ij,ji', Pb, Kb)
        E_elec = E_core + E_J + E_Ka + E_Kb
        E_tot  = E_elec + nuc_rep

        # expectation value of the total spin <S^2>; deviation from the exact
        # value S(S+1) measures spin contamination
        Sz = 0.5 * (nalpha - nbeta)
        s2 = Sz * (Sz + 1.0) + nbeta - np.einsum('ij,ji', Pa @ S, Pb @ S)

        sol = {
            "energy": E_tot,
            "converged": converged,
            "niter": len(energies),
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
            "s2": s2,
            "s2_exact": Sz * (Sz + 1.0),

            "mol": self._mol,
            "forces": None, # Forces: not available for UHF
        }

        return sol

    def __nuclear_repulsion(self) -> float:
        """
        Calculate the repulsion energy between the nuclei
        """
        nuc_rep = 0.0
        for i in range(0, len(self._nuclei)):
            for j in range(i+1, len(self._nuclei)):
                r = np.linalg.norm(np.array(self._nuclei[i][0]) - np.array(self._nuclei[j][0]))
                nuc_rep += self._nuclei[i][1] * self._nuclei[j][1] / r

        return nuc_rep

    def __transformation_matrix(self, S: Mat, ortho: str) -> Mat:
        """
        Construct the matrix X that orthogonalizes the basis, X^T S X = 1

        Canonical orthogonalization drops the eigenvectors of S whose
        eigenvalue is below LINDEP_THRESHOLD, such that X can have fewer
        columns than rows; symmetric orthogonalization keeps all of them.
        """
        # diagonalize S
        s, U = np.linalg.eigh(S)

        if ortho == 'canonical': # perform canonical orthogonalization
            keep = s > LINDEP_THRESHOLD
            X = U[:,keep] @ np.diag(1.0/np.sqrt(s[keep]))
        elif ortho == 'symmetric':
            if s[0] < LINDEP_THRESHOLD:
                warnings.warn(
                    f"The basis set is nearly linearly dependent (smallest overlap "
                    f"eigenvalue {s[0]:.2e}); use canonical orthogonalization.",
                    RuntimeWarning, stacklevel=3,
                )
            X = U @ np.diag(1.0/np.sqrt(s)) @ U.transpose()
        else:
            raise Exception("Invalid orthogonalization option selected: ", ortho)

        return X

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

    def __rhf_forces(self, C, P, e):
        """
        Calculate derivatives of nuclear coordinates to obtain forces.

        The basis functions and nuclei of the calculation are used, such that
        forces are also available when a list of CGFs was supplied as basis.

        Parameters:
        C:      matrix of orbital coefficients
        P:      density matrix
        e:      orbital energies

        Returns array of forces with shape [Natoms, 3]
        """
        # intialization
        integrator = PyQInt()
        cgfs, nuclei = self._cgfs, self._nuclei
        nelec = self._nelec
        forces = np.zeros((len(nuclei),3))
        N = len(cgfs)
        occ = [2 if i < nelec//2 else 0 for i in range(C.shape[1])]
    
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
