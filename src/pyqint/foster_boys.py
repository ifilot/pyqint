# -*- coding: utf-8 -*-

"""
Foster–Boys orbital localization.

This module implements the Foster–Boys procedure for constructing
localized molecular orbitals from canonical Hartree–Fock orbitals.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional

import numpy as np
import numpy.typing as npt
import scipy.optimize

from .pyqint_core import PyQInt


Vec = npt.NDArray[np.float64]
Mat = npt.NDArray[np.float64]


class FosterBoys:
    """
    Foster–Boys orbital localization procedure.

    This class is *stateful* and intended for one localization task.
    Users should interact only via `run()`.
    """

    def __init__(
        self,
        hf_result: Dict[str, Any],
        *,
        seed: Optional[int] = None,
        maxiter: int = 1000,
    ) -> None:
        """
        Parameters
        ----------
        hf_result
            Result dictionary returned by the Hartree–Fock procedure.
        seed
            Random seed for reproducibility.
        maxiter
            Maximum number of Foster–Boys iterations.
        """
        if 'orbe_alpha' in hf_result.keys():
            raise Exception('PopulationAnalysis is not yet supported for UHF')

        # Canonical HF quantities (read-only)
        self._orbc_canonical: Mat = hf_result["orbc"]
        self._orbe_canonical: Vec = hf_result["orbe"]
        self._mol = hf_result["mol"]
        self._nuclei = hf_result["nuclei"]
        self._nelec: int = hf_result["nelec"]
        self._H: Mat = hf_result["fock"]
        self._cgfs = hf_result["cgfs"]
        self._overlap = hf_result["overlap"]
        self._fock = hf_result["fock"]
        self.__density = hf_result["density"]

        # Algorithm parameters
        self._maxiter: int = maxiter
        self._rng = np.random.default_rng(seed)

        # Occupation mask (restricted closed-shell)
        nocc = self._nelec // 2
        self._occ: Vec = np.array(
            [1.0 if i < nocc else 0.0 for i in range(len(self._cgfs))]
        )

        # Precompute dipole tensor (dominant cost)
        self._dipole_tensor: Mat = self._build_dipole_tensor(self._cgfs)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, nr_runners: int = 1) -> Dict[str, Any]:
        """
        Run the Foster–Boys localization.

        Multiple random initializations can be used to reduce the
        probability of converging to a local minimum.

        Parameters
        ----------
        nr_runners
            Number of independent random initializations.

        Returns
        -------
        dict
            Localization result.
        """
        best_result: Optional[Dict[str, Any]] = None
        best_r2: float = -np.inf

        for _ in range(nr_runners):
            result = self._single_runner()
            if result["r2final"] > best_r2:
                best_r2 = result["r2final"]
                best_result = result

        assert best_result is not None
        return best_result

    # ------------------------------------------------------------------
    # Core algorithm
    # ------------------------------------------------------------------

    def _single_runner(self) -> Dict[str, Any]:
        """
        Execute one Foster–Boys optimization run.
        """
        C = self._random_orthogonal_initial_guess(self._orbc_canonical)

        r2_old = 0.0
        for niter in range(self._maxiter):
            C, r2_new = self._mix_orbitals(C)
            if abs(r2_new - r2_old) < 1e-7:
                break
            r2_old = r2_new
        else:
            raise RuntimeError("Foster–Boys localization did not converge.")

        orbe, orbc = self._compute_orbital_energies(C)

        return {
            "orbe": orbe,
            "orbc": orbc,
            "overlap": self._overlap,
            "fock": self._fock,
            "nriter": niter + 1,
            "mol": self._mol,
            "r2start": self._compute_r2(self._orbc_canonical),
            "r2final": self._compute_r2(orbc),
            "nelec": self._nelec,
            "cgfs": self._cgfs,
            "nuclei": self._nuclei,
            "density": self.__density,
        }

    # ------------------------------------------------------------------
    # Foster–Boys mechanics
    # ------------------------------------------------------------------

    def _mix_orbitals(self, C: Mat) -> tuple[Mat, float]:
        """
        Perform pairwise orbital rotations to maximize the Boys functional.
        """
        nocc = self._nelec // 2
        r2_start = self._compute_r2(C)
        r2_best = r2_start

        for i in range(nocc):
            for j in range(i + 1, nocc):
                res = scipy.optimize.minimize(
                    self._evaluate_rotation,
                    0.0,
                    args=(C, i, j),
                    bounds=[(-np.pi, np.pi)],
                    tol=1e-12,
                )

                alpha = res.x[0]
                C_new = self._rotate_pair(C, i, j, alpha)

                r2 = self._compute_r2(C_new)
                if r2 > r2_best:
                    C = C_new
                    r2_best = r2

        return C, r2_best

    def _evaluate_rotation(self, alpha: float, C: Mat, i: int, j: int) -> float:
        """
        Objective function for a 2×2 orbital rotation.
        """
        C_new = self._rotate_pair(C, i, j, alpha)
        return -self._compute_r2(C_new)

    def _compute_r2(self, C: Mat) -> float:
        """
        Compute the Foster–Boys localization functional.
        """
        dip = np.einsum("ji,ki,jkl->il", C, C, self._dipole_tensor)
        return float(np.einsum("ij,i->", dip**2, self._occ))

    # ------------------------------------------------------------------
    # Linear algebra helpers
    # ------------------------------------------------------------------

    def _rotate_pair(self, C: Mat, i: int, j: int, alpha: float) -> Mat:
        """
        Apply a 2×2 unitary rotation to orbitals i and j.
        """
        C_new = C.copy()
        C_new[:, i] = np.cos(alpha) * C[:, i] + np.sin(alpha) * C[:, j]
        C_new[:, j] = -np.sin(alpha) * C[:, i] + np.cos(alpha) * C[:, j]
        return C_new

    def _random_orthogonal_initial_guess(self, C: Mat, nops: int = 100) -> Mat:
        """
        Generate a randomized orthogonal transformation of occupied orbitals.
        """
        nocc = self._nelec // 2
        for _ in range(nops):
            i, j = self._rng.choice(nocc, size=2, replace=False)
            angle = self._rng.uniform(0.0, 2.0 * np.pi)
            C = self._rotate_pair(C, i, j, angle)
        return C

    # ------------------------------------------------------------------
    # Precomputation
    # ------------------------------------------------------------------

    def _build_dipole_tensor(self, cgfs: list) -> Mat:
        """
        Precompute the dipole integral tensor ⟨χ_i | r_k | χ_j⟩.
        """
        n = len(cgfs)
        tensor = np.zeros((n, n, 3))
        integrator = PyQInt()

        for i, c1 in enumerate(cgfs):
            for j, c2 in enumerate(cgfs):
                for k in range(3):
                    tensor[i, j, k] = integrator.dipole(c1, c2, k, 0.0)

        return tensor

    def _compute_orbital_energies(self, C: Mat) -> tuple[Vec, Mat]:
        """
        Compute MO energies in the localized basis.
        """
        energies = np.array([C[:, i] @ self._H @ C[:, i] for i in range(C.shape[1])])
        idx = np.argsort(energies)
        return energies[idx], C[:, idx]