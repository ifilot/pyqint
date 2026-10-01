# -*- coding: utf-8 -*-

"""
Foster-Boys orbital localization.

This module implements the Foster-Boys procedure for constructing
localized molecular orbitals from canonical Hartree-Fock orbitals.
Both restricted (RHF) and unrestricted (UHF) results are supported; for
the latter, the alpha and beta orbitals are localized independently.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional

import numpy as np
import numpy.typing as npt
import scipy.optimize

from .pyqint_core import PyQInt
from .spin import get_spin_channel, is_unrestricted


Vec = npt.NDArray[np.float64]
Mat = npt.NDArray[np.float64]


class FosterBoys:
    """
    Foster-Boys orbital localization procedure.

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
            Result dictionary returned by the Hartree-Fock procedure.
        seed
            Random seed for reproducibility.
        maxiter
            Maximum number of Foster-Boys iterations.
        """
        self._res = hf_result
        self._unrestricted: bool = is_unrestricted(hf_result)

        # Canonical HF quantities shared by all spin channels (read-only)
        self._mol = hf_result["mol"]
        self._nuclei = hf_result["nuclei"]
        self._nelec: int = hf_result["nelec"]
        self._cgfs = hf_result["cgfs"]
        self._overlap = hf_result["overlap"]

        # Algorithm parameters
        self._maxiter: int = maxiter
        self._rng = np.random.default_rng(seed)

        # Precompute dipole tensor (dominant cost)
        self._dipole_tensor: Mat = self._build_dipole_tensor(self._cgfs)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, nr_runners: int = 1) -> Dict[str, Any]:
        """
        Run the Foster-Boys localization.

        Multiple random initializations can be used to reduce the
        probability of converging to a local minimum.

        Parameters
        ----------
        nr_runners
            Number of independent random initializations.

        Returns
        -------
        dict
            Localization result. For RHF input, the keys 'orbc', 'orbe',
            'fock', 'nriter', 'r2start' and 'r2final' are provided. For UHF
            input, these keys carry an '_alpha' or '_beta' suffix, mirroring
            the output of the UHF procedure.
        """
        result: Dict[str, Any] = {
            "overlap": self._overlap,
            "mol": self._mol,
            "nelec": self._nelec,
            "cgfs": self._cgfs,
            "nuclei": self._nuclei,
            "density": self._res["density"],
        }

        if not self._unrestricted:
            result.update(self._localize_channel(None, nr_runners))
            return result

        for spin in ("alpha", "beta"):
            channel = self._localize_channel(spin, nr_runners)
            for key, value in channel.items():
                result[key + "_" + spin] = value

        for key in ("nalpha", "nbeta", "multiplicity"):
            result[key] = self._res[key]

        return result

    # ------------------------------------------------------------------
    # Core algorithm
    # ------------------------------------------------------------------

    def _localize_channel(self, spin: Optional[str], nr_runners: int) -> Dict[str, Any]:
        """
        Localize the occupied orbitals of a single spin channel, retaining
        the best result out of `nr_runners` random initializations.
        """
        channel = get_spin_channel(self._res, spin)
        C0: Mat = channel["orbc"]
        F: Mat = channel["fock"]
        nocc: int = channel["nocc"]

        best_result: Optional[Dict[str, Any]] = None
        best_r2: float = -np.inf

        for _ in range(nr_runners):
            result = self._single_runner(C0, F, nocc)
            if result["r2final"] > best_r2:
                best_r2 = result["r2final"]
                best_result = result

        assert best_result is not None
        best_result["density"] = channel["density"]
        return best_result

    def _single_runner(self, C0: Mat, F: Mat, nocc: int) -> Dict[str, Any]:
        """
        Execute one Foster-Boys optimization run.
        """
        # with fewer than two occupied orbitals, there is nothing to rotate
        if nocc < 2:
            C = C0
            niter = -1
        else:
            C = self._random_orthogonal_initial_guess(C0, nocc)

            r2_old = 0.0
            for niter in range(self._maxiter):
                C, r2_new = self._mix_orbitals(C, nocc)
                if abs(r2_new - r2_old) < 1e-7:
                    break
                r2_old = r2_new
            else:
                raise RuntimeError("Foster-Boys localization did not converge.")

        orbe, orbc = self._compute_orbital_energies(C, F)

        return {
            "orbe": orbe,
            "orbc": orbc,
            "fock": F,
            "nriter": niter + 1,
            "r2start": self._compute_r2(C0, nocc),
            "r2final": self._compute_r2(orbc, nocc),
        }

    # ------------------------------------------------------------------
    # Foster-Boys mechanics
    # ------------------------------------------------------------------

    def _mix_orbitals(self, C: Mat, nocc: int) -> tuple[Mat, float]:
        """
        Perform pairwise orbital rotations to maximize the Boys functional.
        """
        r2_start = self._compute_r2(C, nocc)
        r2_best = r2_start

        for i in range(nocc):
            for j in range(i + 1, nocc):
                res = scipy.optimize.minimize(
                    self._evaluate_rotation,
                    0.0,
                    args=(C, i, j, nocc),
                    bounds=[(-np.pi, np.pi)],
                    tol=1e-12,
                )

                alpha = res.x[0]
                C_new = self._rotate_pair(C, i, j, alpha)

                r2 = self._compute_r2(C_new, nocc)
                if r2 > r2_best:
                    C = C_new
                    r2_best = r2

        return C, r2_best

    def _evaluate_rotation(self, alpha: float, C: Mat, i: int, j: int,
                           nocc: int) -> float:
        """
        Objective function for a 2×2 orbital rotation.
        """
        C_new = self._rotate_pair(C, i, j, alpha)
        return -self._compute_r2(C_new, nocc)

    def _compute_r2(self, C: Mat, nocc: int) -> float:
        """
        Compute the Foster-Boys localization functional over the
        `nocc` occupied orbitals.
        """
        Cocc = C[:, :nocc]
        dip = np.einsum("ji,ki,jkl->il", Cocc, Cocc, self._dipole_tensor)
        return float(np.sum(dip**2))

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

    def _random_orthogonal_initial_guess(self, C: Mat, nocc: int,
                                         nops: int = 100) -> Mat:
        """
        Generate a randomized orthogonal transformation of occupied orbitals.
        """
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

    def _compute_orbital_energies(self, C: Mat, F: Mat) -> tuple[Vec, Mat]:
        """
        Compute MO energies in the localized basis.
        """
        energies = np.array([C[:, i] @ F @ C[:, i] for i in range(C.shape[1])])
        idx = np.argsort(energies)
        return energies[idx], C[:, idx]