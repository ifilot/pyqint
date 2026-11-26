# -*- coding: utf-8 -*-

from __future__ import annotations

import time
import warnings
from functools import partial
from operator import index
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt
from scipy import constants
from scipy.sparse.linalg import ArpackNoConvergence

from .._memory import check_hf_gradient_memory
from ..pyqint_core import PyQInt
from ..structure.element import ELEMENTS_DATA
from ..structure.molecule import Molecule
from .initial_guess import LINDEP_THRESHOLD, sad_density
from .stability import lowest_rhf_mode, rotate_rhf_occupied
from .vibrations import displace_basis, harmonic_analysis

Vec = npt.NDArray[np.float64]
Mat = npt.NDArray[np.float64]

# nucleus entry: (position, charge)
Nucleus = Tuple[Vec, int]

AU_TO_DEBYE = constants.elementary_charge * constants.value("Bohr radius") / (1e-21 / constants.c)

# couple of hardcoded variables for the DIIS algorithm
SUBSPACE_LENGTH = 4
SUBSPACE_START = 1

STABILITY_THRESHOLD = -1e-5

class HF:
    """
    Perform Hartree-Fock calculations
    """
    _eri_storage = "dense"
    _eri_engine = "hellsing"

    def __init__(
        self,
        mol: Molecule,
        basis: Union[str, Sequence[Any]],
        *,
        eri_storage: str = "dense",
        eri_engine: str = "hellsing",
    ) -> None:
        """
        Initialize the Hartree-Fock solver.

        Parameters
        ----------
        mol
            Molecule object.
        basis
            Basis set name or a list of CGF objects.
        eri_storage
            Use ``"dense"`` to return the four-dimensional ERI tensor, or
            ``"packed"`` to contract unique ERIs directly in RHF and UHF.
            Packed results contain ``packed_tei`` and set ``tetensor`` to None.
        eri_engine
            PyQInt engine used for the SCF integrals, such as ``"auto"``.
        """
        if eri_storage not in ("dense", "packed"):
            raise ValueError(f"Unknown ERI storage: {eri_storage!r}")
        self._eri_storage = eri_storage
        self._eri_engine = eri_engine
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
        *,
        guess: str = "sad",
        stability: bool = True,
        stability_max_restarts: int = 3,
        calc_dipole: bool = False,
        dipole_origin: Sequence[float] = (0.0, 0.0, 0.0),
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
            Include the analytic nuclear energy gradient under ``forces``.
            The physical nuclear force is the negative of this gradient.
        itermax
            Maximum number of SCF iterations per attempt.
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
        guess
            Initial density used to build the first Fock matrix. ``"sad"``
            adds the densities of neutral atoms in their own basis functions.
            Each open atomic shell shares its electrons equally among its
            degenerate orbitals. This default accounts for electronic
            repulsion from the first Fock build and can reach lower-energy
            solutions. The first molecular diagonalization sets the actual
            occupations, including for ions. ``"core"`` starts with zero
            density, so the first Fock matrix is the core Hamiltonian T + V.
            Supplied ``orbc_init`` takes precedence. Multicenter or
            off-nucleus basis functions fall back to the core guess.
        stability
            Check the converged determinant for a downhill real RHF orbital
            rotation. Follow a downhill mode and restart if one is found.
            Use ``guess="core", stability=False`` to reproduce the original
            core-start SCF trajectory.
        stability_max_restarts
            Maximum number of stability restarts. Set to zero to diagnose
            an instability without following it.
        calc_dipole
            Include the molecular dipole vector in debye under ``dipole``.
        dipole_origin
            Cartesian origin of the dipole moment in bohr.

        Returns
        -------
        dict
            Includes ``stable`` (True or False after a completed check,
            otherwise None), ``stability_eigenvalue`` (lowest eigenvalue of
            A + B, or None if unavailable), and ``stability_restarts``.
            The iteration history includes all SCF attempts. With zero
            iterations the result contains the chosen initial density.

        Raises
        ------
        ValueError
            If the number of electrons is odd or exceeds the capacity of the
            molecular orbitals, the guess is unknown, or the restart limit
            is negative.
        """
        if self._nelec % 2 != 0:
            raise ValueError(
                f"RHF requires an even number of electrons, not {self._nelec}; "
                "use uhf() for open-shell systems."
            )
        if guess not in ("core", "sad"):
            raise ValueError(f"Unknown initial guess: {guess!r}")
        stability_max_restarts = index(stability_max_restarts)
        if stability_max_restarts < 0:
            raise ValueError("stability_max_restarts must be nonnegative.")
        if calc_forces:
            tolerance = self.__gradient_tolerance(itermax, tolerance)
        if tolerance_grad is None:
            tolerance_grad = np.sqrt(tolerance)

        # create empty dictionary for time tracking statistics
        time_stats = {}

        # build integrals
        integrator = PyQInt(eri_engine=self._eri_engine)
        start = time.time()
        packed = self._eri_storage == "packed"
        S, T, V, tei = integrator.build_integrals_openmp(self._cgfs, self._nuclei, packed=packed)
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

        if orbc_init is None:
            initial_density = sad_density(self._cgfs, self._nuclei) if guess == "sad" else None
        else:
            occ_init = [2 if i < self._nelec // 2 else 0 for i in range(orbc_init.shape[1])]
            initial_density = np.einsum('ik,jk,k->ij', orbc_init, orbc_init, occ_init)
        P = np.zeros_like(S) if initial_density is None else initial_density

        hcore = T + V
        coulomb_exchange = partial(self.__coulomb_exchange, integrator, tei)
        energies = []
        time_stats['iterations'] = []
        time_stats['self_convergence'] = 0.0
        stable = lowest_eigenvalue = None
        restarts = 0

        for attempt in range(1 + stability_max_restarts):
            start = time.time()
            P, F, J, K, energy, orbe, C, converged = self.__rhf_scf(
                hcore, S, X, nuc_rep, coulomb_exchange, P, occ,
                itermax, use_diis, verbose, tolerance, tolerance_grad,
                initial_density is not None or attempt > 0, energies, time_stats['iterations'],
            )
            time_stats['self_convergence'] += time.time() - start
            if not stability or not converged:
                stable = lowest_eigenvalue = None
                if restarts:
                    warnings.warn("RHF stability restart did not converge.", RuntimeWarning, stacklevel=2)
                break

            try:
                lowest_eigenvalue, direction = lowest_rhf_mode(C, F, self._nelec // 2, coulomb_exchange)
            except ArpackNoConvergence:
                stable = lowest_eigenvalue = None
                warnings.warn("RHF stability eigensolver did not converge.", RuntimeWarning, stacklevel=2)
                break

            stable = lowest_eigenvalue is None or lowest_eigenvalue >= STABILITY_THRESHOLD
            if stable or attempt == stability_max_restarts:
                break
            restarted_density = self.__follow_rhf_mode(C, direction, energy, hcore, nuc_rep, coulomb_exchange)
            if restarted_density is None:
                break
            P = restarted_density
            restarts += 1
            if verbose:
                print(f"Restarting RHF along an unstable mode (eigenvalue {lowest_eigenvalue:.2e} Ht).")

        if stable is False:
            warnings.warn(
                f"RHF solution is internally unstable (lowest eigenvalue {lowest_eigenvalue:.2e} Ht) "
                f"after {restarts} stability restart(s).",
                RuntimeWarning, stacklevel=2,
            )

        gradient = None
        if calc_forces:
            ew_density = np.einsum('ik,jk,k,k->ij', C, C, occ, orbe)
            gradient = self.__nuclear_gradient(P, ew_density)

        # build solution dictionary
        sol = {
            "energy": energy,
            "converged": converged,
            "stable": stable,
            "stability_eigenvalue": lowest_eigenvalue,
            "stability_restarts": restarts,
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
            'tetensor': None if packed else tei,
            "time_stats" : time_stats,
            "ecore": np.sum(P * (T + V)),
            "ekin": np.einsum('ij,ji', T, P),
            "enuc": np.einsum('ij,ji', V, P),
            "erep": 0.5 * np.einsum('ij,ij', P, J),
            "ex": -0.25 * np.einsum('ij,ij', P, K),
            "enucrep": nuc_rep,
            "nelec": self._nelec,
            "mol": self._mol,
            "forces": gradient,
        }

        if packed:
            sol["packed_tei"] = tei
        if calc_dipole:
            sol["dipole"] = self.__molecular_dipole(integrator, P, dipole_origin)
        return sol

    def __rhf_scf(self, hcore, S, X, nuc_rep, coulomb_exchange, P, occ,
                  itermax, use_diis, verbose, tolerance, tolerance_grad,
                  include_first_diis, energies, iteration_times):
        """Run one RHF attempt with a fresh DIIS subspace."""
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
            J, K = coulomb_exchange(P)
            G = J - 0.5 * K

            # build Fock matrix
            F = hcore + G

            # transform Fock matrix
            Fprime = X.transpose() @ F @ X

            # diagonalize F
            orbe, Cprime = np.linalg.eigh(Fprime)

            # back-transform
            C = X @ Cprime

            # calculate energy E
            energy = 0.5 * np.einsum('ij,ji', P, hcore + F)

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
            if niter > 0 or include_first_diis:
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
            iteration_times.append(iterend - iterstart)

            # print info for this iteration
            if verbose:
                print("Iteration: %2i | Energy: %12.6f Ht | Error: %8.2e | Time: %6.4f s" % (niter, energy, errmax, iteration_times[-1]))

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

        if itermax > 0 and not converged:
            warnings.warn(
                f"RHF did not converge within {itermax} iterations "
                f"(last energy change {ediff if niter > 4 else float('nan'):.2e} Ht, "
                f"error {errmax:.2e}).",
                RuntimeWarning, stacklevel=3,
            )

        # build the density matrix from the final orbitals and the Fock
        # matrix and energy that belong to it; without any iterations, the
        # initial guess density is retained
        if C is not None:
            P = np.einsum('ik,jk,k->ij', C, C, occ)
        J, K = coulomb_exchange(P)
        G = J - 0.5 * K
        F = hcore + G
        energy = 0.5 * np.einsum('ji,ij', P, hcore + F) + nuc_rep
        if energies:
            energies[-1] = energy

        return P, F, J, K, energy, orbe, C, converged

    @staticmethod
    def __follow_rhf_mode(orbitals, direction, energy, hcore, nuc_rep, coulomb_exchange):
        """Find a downhill step along the unstable orbital rotation."""
        for step in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125):
            best_density = None
            best_energy = energy
            for sign in (1, -1):
                occupied = rotate_rhf_occupied(orbitals, sign * step * direction)
                density = 2 * (occupied @ occupied.T)
                coulomb, exchange = coulomb_exchange(density)
                fock = hcore + coulomb - 0.5 * exchange
                trial_energy = 0.5 * np.einsum('ij,ji', density, hcore + fock) + nuc_rep
                if trial_energy < best_energy:
                    best_energy = trial_energy
                    best_density = density
            if best_density is not None:
                return best_density
        return None

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
            *,
            guess: str = "sad",
            calc_forces: bool = False,
            calc_dipole: bool = False,
            dipole_origin: Sequence[float] = (0.0, 0.0, 0.0),
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
            Number of electrons. Defaults to the molecular electron count
            after applying its charge.
        tolerance_grad
            Threshold for the commutators FPS - SPF; defaults to
            sqrt(tolerance).
        guess
            ``"sad"`` (default) starts from the neutral atomic densities
            described in :meth:`rhf`, split equally between alpha and beta.
            This accounts for electronic repulsion in the first Fock build.
            The first molecular diagonalization sets the two spin electron
            counts from the requested charge and multiplicity. ``"core"``
            starts both spin densities at zero. Supplied ``orbc_init`` takes
            precedence. Multicenter and off-nucleus basis functions fall back
            to the core guess.
        calc_forces
            Include the analytic nuclear energy gradient under ``forces``.
            The physical nuclear force is the negative of this gradient.
        calc_dipole
            Include the molecular dipole vector in debye under ``dipole``.
        dipole_origin
            Cartesian origin of the dipole moment in bohr.

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

        if guess not in ("core", "sad"):
            raise ValueError(f"Unknown initial guess: {guess!r}")

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

        if calc_forces:
            tolerance = self.__gradient_tolerance(itermax, tolerance)
        if tolerance_grad is None:
            tolerance_grad = np.sqrt(tolerance)

        time_stats = {}

        # build integrals
        integrator = PyQInt(eri_engine=self._eri_engine)
        start = time.time()
        packed = self._eri_storage == "packed"
        S, T, V, tei = integrator.build_integrals_openmp(self._cgfs, self._nuclei, packed=packed)
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
        include_first_diis = True
        if orbc_init is None:
            density = sad_density(self._cgfs, self._nuclei) if guess == "sad" else None
            include_first_diis = density is not None
            Pa = np.zeros_like(S) if density is None else 0.5 * density
            Pb = Pa.copy()
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
            J, Ka, Kb = self.__coulomb_exchange(integrator, tei, Pa, Pb)

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
            if niter > 0 or include_first_diis:
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
            J_new, Ka_new, Kb_new = self.__coulomb_exchange(integrator, tei, Pa_new, Pb_new)

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
        J, Ka, Kb = self.__coulomb_exchange(integrator, tei, Pa, Pb)
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

        gradient = None
        if calc_forces:
            # At self-consistency P F P is the energy-weighted density of a spin channel.
            ew_density = Pa @ Fa @ Pa + Pb @ Fb @ Pb
            gradient = self.__nuclear_gradient(Ptot, ew_density, Pa, Pb)

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
            "tetensor": None if packed else tei,
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
            "forces": gradient,
        }

        if packed:
            sol["packed_tei"] = tei
        if calc_dipole:
            sol["dipole"] = self.__molecular_dipole(integrator, Ptot, dipole_origin)
        return sol

    def vibrational_frequencies(
        self,
        *,
        multiplicity: Optional[int] = None,
        step: float = 1e-3,
        masses: Optional[Sequence[float]] = None,
        itermax: int = 100,
        use_diis: bool = True,
        tolerance: float = 1e-12,
        tolerance_grad: float = 1e-8,
        ortho: str = "canonical",
        orbc_init: Optional[Union[Mat, Dict[str, Mat]]] = None,
    ) -> Dict[str, Any]:
        """Calculate harmonic frequencies at the supplied geometry.

        Central differences of analytic energy gradients build the Cartesian
        Hessian. Use a stationary geometry to interpret the harmonic modes.

        Parameters
        ----------
        multiplicity
            Spin multiplicity for UHF. The default ``None`` selects RHF.
        step
            Cartesian displacement in bohr.
        masses
            Atomic masses in unified atomic mass units, in molecule order.
            Defaults to the element masses stored by PyQInt.
        itermax, use_diis, tolerance, tolerance_grad, ortho
            SCF options for the reference and all displaced geometries.
        orbc_init
            Optional initial orbitals for the reference calculation, using
            the same format as ``rhf`` or ``uhf`` for the chosen method.

        Returns
        -------
        dict
            ``frequencies`` in cm⁻¹, with negative values denoting imaginary
            frequencies, normalized Cartesian ``normal_modes`` with shape
            (modes, atoms, 3), the Cartesian ``hessian`` in hartree/bohr²,
            and the ``masses`` used for mass weighting.
        """
        if not np.isfinite(step) or step <= 0:
            raise ValueError("Vibrational displacement must be finite and positive.")
        if not self._nuclei:
            raise ValueError("Vibrational analysis requires at least one nucleus.")
        if itermax < 1:
            raise ValueError("Vibrational analysis requires at least one SCF iteration.")
        coordinates = np.asarray([position for position, _ in self._nuclei], dtype=float)
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("Nuclear coordinates must be finite.")
        if masses is None:
            masses = [ELEMENTS_DATA[charge]["atomic-mass"] for _, charge in self._nuclei]
        masses = np.array(masses, dtype=float)
        if masses.shape != (len(coordinates),) or not np.all(np.isfinite(masses) & (masses > 0)):
            raise ValueError("Supply one finite, positive mass per nucleus.")

        tolerance = self.__gradient_tolerance(itermax, tolerance)
        options = dict(itermax=itermax, use_diis=use_diis, tolerance=tolerance,
                       tolerance_grad=tolerance_grad, ortho=ortho)
        method = "rhf" if multiplicity is None else "uhf"
        if multiplicity is not None:
            options["multiplicity"] = multiplicity
        reference = getattr(self, method)(orbc_init=orbc_init, **options)
        if not reference["converged"]:
            raise RuntimeError(f"Reference {method.upper()} calculation did not converge for vibrational analysis.")
        if multiplicity is None:
            orbitals = reference["orbc"]
        else:
            orbitals = {"alpha": reference["orbc_alpha"], "beta": reference["orbc_beta"]}
        del reference

        def gradient(displacement):
            molecule = Molecule(self._mol.get_name())
            molecule.set_charge(self._mol.get_charge())
            for (symbol, _), position in zip(self._mol, coordinates + displacement):
                molecule.add_atom(symbol, *position)
            basis = self._basis
            if not isinstance(basis, str):
                basis = displace_basis(self._cgfs, coordinates, displacement)
            solver = HF(molecule, basis, eri_storage=self._eri_storage, eri_engine=self._eri_engine)
            result = getattr(solver, method)(calc_forces=True, orbc_init=orbitals, **options)
            if not result["converged"]:
                atom, axis = divmod(np.argmax(np.abs(displacement)), 3)
                delta = displacement[atom, axis]
                raise RuntimeError(
                    f"{method.upper()} did not converge for atom {atom + 1} displaced {delta:+g} bohr along {'xyz'[axis]}."
                )
            # The legacy result key "forces" contains the energy gradient.
            return result["forces"].ravel()

        hessian = np.empty((coordinates.size, coordinates.size))
        displacement = np.zeros_like(coordinates)
        for i in range(coordinates.size):
            displacement.flat[i] = step
            plus = gradient(displacement)
            displacement.flat[i] = -step
            minus = gradient(displacement)
            displacement.flat[i] = 0
            hessian[:, i] = (plus - minus) / (2 * step)
        hessian = (hessian + hessian.T) / 2
        frequencies, modes = harmonic_analysis(hessian, coordinates, masses)
        return {"frequencies": frequencies, "normal_modes": modes, "hessian": hessian, "masses": masses}

    def __molecular_dipole(self, integrator, density, origin):
        """Combine nuclear and electronic dipoles relative to the same origin."""
        origin = np.asarray(origin, dtype=np.float64)
        integrals = integrator.build_dipole_openmp(self._cgfs, origin=origin)
        nuclear = np.zeros(3)
        for position, charge in self._nuclei:
            nuclear += charge * (position - origin)
        electronic = np.einsum('ij,dij->d', density, integrals)
        return AU_TO_DEBYE * (nuclear - electronic)

    def __coulomb_exchange(self, integrator, tei, density_alpha, density_beta=None):
        """Contract the ERIs with the restricted or spin-resolved densities."""
        if self._eri_storage == "packed":
            return integrator.contract_packed_tei(tei, density_alpha, density_beta)

        density = density_alpha if density_beta is None else density_alpha + density_beta
        J = np.einsum('kl,ijlk->ij', density, tei)
        Ka = np.einsum('kl,iklj->ij', density_alpha, tei)
        if density_beta is None:
            return J, Ka
        Kb = np.einsum('kl,iklj->ij', density_beta, tei)
        return J, Ka, Kb

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

    def __gradient_tolerance(self, itermax: int, tolerance: float) -> float:
        """Check derivative workspace and tighten the SCF energy tolerance."""
        if itermax < 1:
            raise ValueError("Calculating forces requires at least one SCF iteration.")
        check_hf_gradient_memory(len(self._cgfs), len(self._nuclei), packed=self._eri_storage == "packed")
        return min(tolerance, 1e-12)

    def __nuclear_gradient(self, density: Mat, ew_density: Mat,
                           density_a: Optional[Mat] = None, density_b: Optional[Mat] = None) -> Mat:
        """Contract analytic derivative integrals into the Cartesian energy gradient."""
        integrator = PyQInt()
        overlap, kinetic, attraction, repulsion = integrator.build_geometric_derivatives_openmp(
            self._cgfs, self._nuclei,
        )
        gradient = np.einsum('adij,ij->ad', kinetic, density)
        gradient += np.einsum('adij,ij->ad', attraction, density)
        gradient -= np.einsum('adij,ij->ad', overlap, ew_density)

        # Derivative ERIs share the packed ordering used for the SCF integrals.
        for atom, derivatives in enumerate(repulsion):
            for axis, packed in enumerate(derivatives):
                if density_b is None:
                    J, K = integrator.contract_packed_tei(packed, density)
                    exchange = 0.5 * np.sum(density * K)
                else:
                    J, Ka, Kb = integrator.contract_packed_tei(packed, density_a, density_b)
                    exchange = np.sum(density_a * Ka) + np.sum(density_b * Kb)
                gradient[atom, axis] += 0.5 * (np.sum(density * J) - exchange)

        # Each nuclear pair contributes equal and opposite energy derivatives.
        for atom, (position, charge) in enumerate(self._nuclei):
            for other in range(atom + 1, len(self._nuclei)):
                other_position, other_charge = self._nuclei[other]
                delta = other_position - position
                distance = np.linalg.norm(delta)
                if distance > 1e-4:
                    derivative = charge * other_charge * delta / distance ** 3
                    gradient[atom] += derivative
                    gradient[other] -= derivative
        return gradient
