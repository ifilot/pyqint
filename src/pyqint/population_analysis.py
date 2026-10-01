# -*- coding: utf-8 -*-

from __future__ import annotations

from typing import Dict, Any, List, Optional

import numpy as np
import numpy.typing as npt

from .spin import get_spin_channel, is_unrestricted

Vec = npt.NDArray[np.float64]
Mat = npt.NDArray[np.float64]

class PopulationAnalysis:
    """
    Population Analysis class

    This class implements:
      - Mulliken population analysis
      - Lowkin population analysis
      - MOOP: Molecular Orbital Overlap Population
      - MOHP: Molecular Orbital Hamilton Population
      - MOBI: Molecular Orbital Bond Index

    The analysis operates on the output of either a restricted (RHF) or an
    unrestricted (UHF) Hartree–Fock calculation. For UHF results, the
    orbital-resolved analyses (MOOP, MOHP, MOBI) require the spin channel
    to be specified via ``spin='alpha'`` or ``spin='beta'``.
    """

    def __init__(self, res: Dict[str, Any]) -> None:
        """
        Parameters
        ----------
        res
            Result dictionary returned by a Hartree–Fock calculation
            (either RHF or UHF).
        """
        self.res = res

        # Whether the result originates from an unrestricted calculation
        self.unrestricted: bool = is_unrestricted(res)

        # Overlap matrix
        self.S: Mat = res["overlap"]

        # Total density matrix
        self.P: Mat = res["density"]

        # Spin density matrix (alpha - beta); vanishes for RHF
        if self.unrestricted:
            self.Pspin: Mat = res["density_alpha"] - res["density_beta"]
        else:
            self.Pspin = np.zeros_like(self.P)

        # Number of electrons
        self.nelec: int = res["nelec"]

        # Nuclear positions: [(position, charge), ...]
        self.nuclei = res["nuclei"]

        # Contracted Gaussian basis functions
        self.cgfs = res["cgfs"]

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def mulliken(self, n: int) -> float:
        """
        Perform Mulliken Population Analysis

        n: index of nucleus

        Returns the atomic charge of nucleus n.
        """
        return self.nuclei[n][1] - self._mulliken_population(n, self.P)

    def lowdin(self, n: int) -> float:
        """
        Perform Löwdin Population Analysis

        n: index of nucleus

        Returns the atomic charge of nucleus n.
        """
        return self.nuclei[n][1] - self._lowdin_population(n, self.P)

    def mulliken_spin(self, n: int) -> float:
        """
        Compute the Mulliken spin population (N_alpha - N_beta) of a nucleus.

        n: index of nucleus

        For RHF results the spin population is identically zero.
        """
        return self._mulliken_population(n, self.Pspin)

    def lowdin_spin(self, n: int) -> float:
        """
        Compute the Löwdin spin population (N_alpha - N_beta) of a nucleus.

        n: index of nucleus

        For RHF results the spin population is identically zero.
        """
        return self._lowdin_population(n, self.Pspin)

    def moop(self, n1: int, n2: int, spin: Optional[str] = None) -> Vec:
        """
        Compute the Molecular Orbital Overlap Population (MOOP).

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.
        spin
            Spin channel ('alpha' or 'beta'); required for UHF results
            and must be omitted for RHF results.

        Returns
        -------
        ndarray
            MOOP values for each molecular orbital.
        """
        return self._population_analysis(n1, n2, 'overlap', spin)

    def mohp(self, n1: int, n2: int, spin: Optional[str] = None) -> Vec:
        """
        Compute the Molecular Orbital Hamilton Population (MOHP).

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.
        spin
            Spin channel ('alpha' or 'beta'); required for UHF results
            and must be omitted for RHF results.

        Returns
        -------
        ndarray
            MOHP values for each molecular orbital.
        """
        return self._population_analysis(n1, n2, 'fock', spin)

    def mobi(self, n1: int, n2: int, spin: Optional[str] = None) -> Vec:
        """
        Compute the Molecular Orbital Bond Index (MOBI).

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.
        spin
            Spin channel ('alpha' or 'beta'); required for UHF results
            and must be omitted for RHF results.

        Returns
        -------
        ndarray
            MOBI values for each molecular orbital.
        """
        return self._population_analysis(n1, n2, 'density', spin)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _population_analysis(self, n1: int, n2: int, kind: str,
                             spin: Optional[str]) -> Vec:
        """
        Shared implementation of MOHP/MOOP/MOBI.

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.
        kind
            Either 'overlap' (S), 'fock' (H) or 'density' (P).
        spin
            Spin channel for UHF results, None for RHF results.

        Returns
        -------
        ndarray
            Population coefficients per molecular orbital.

        Notes
        -----
        For UHF results each spin orbital holds a single electron and the
        prefactor is halved with respect to RHF. For MOBI, the spin density
        enters with a factor 2, consistent with the open-shell Mayer bond
        order. As such, for a closed-shell system the alpha and beta
        contributions sum to the RHF result.
        """
        if n1 == n2:
            raise ValueError(
                "Population analysis requires two distinct atoms."
            )

        channel = get_spin_channel(self.res, spin)
        orbc = channel["orbc"]

        if kind == 'overlap':
            matrix = self.S
        elif kind == 'fock':
            matrix = channel["fock"]
        else:
            # scale spin densities to the closed-shell convention
            matrix = channel["density"] * (2.0 / channel["occ_factor"])

        # Determine basis functions belonging to each nucleus
        idx1, idx2 = self._basis_indices_for_atoms(n1, n2)

        C1 = orbc[idx1, :]
        C2 = orbc[idx2, :]
        M12 = matrix[np.ix_(idx1, idx2)]
        coeff = channel["occ_factor"] * np.einsum('ik,ij,jk->k', C1, M12, C2,
                                                  optimize=True)

        return coeff

    def _basis_indices_for_atom(self, n: int) -> List[int]:
        """
        Determine which basis functions are centered on nucleus n.
        """
        nuc = self.nuclei[n][0]
        return [i for i, cgf in enumerate(self.cgfs)
                if np.linalg.norm(cgf.p - nuc) < 1e-3]

    def _mulliken_population(self, n: int, P: Mat) -> float:
        """
        Number of electrons assigned to nucleus n by Mulliken partitioning
        of the density matrix P.
        """
        # determine the overlap weighted density matrix
        overlap_density = P @ self.S

        # add up electrons in basis functions localized on nucleus
        return float(sum(overlap_density[i, i]
                         for i in self._basis_indices_for_atom(n)))

    def _lowdin_population(self, n: int, P: Mat) -> float:
        """
        Number of electrons assigned to nucleus n by Löwdin partitioning
        of the density matrix P.
        """
        # diagonalize S
        s, U = np.linalg.eigh(self.S)

        # construct transformation matrix X, using Löwdin orthogonalization
        X = U @ np.diag(np.sqrt(s)) @ U.transpose()

        # determine P in (local) orthonormalized basis
        P_prime = X @ P @ X

        # add up electrons in basis functions localized on nucleus
        return float(sum(P_prime[i, i]
                         for i in self._basis_indices_for_atom(n)))

    def _basis_indices_for_atoms(self, n1: int, n2: int) -> tuple[List[int], List[int]]:
        """
        Determine which basis functions belong to two nuclei.

        Basis functions are assigned to atoms by comparing their centers
        to nuclear positions within a small tolerance.

        Parameters
        ----------
        n1, n2
            Indices of the nuclei.

        Returns
        -------
        (list, list)
            Indices of basis functions centered on atom n1 and n2.
        """
        nuc1 = self.nuclei[n1][0]
        nuc2 = self.nuclei[n2][0]

        idx1: List[int] = []
        idx2: List[int] = []

        for i, cgf in enumerate(self.cgfs):
            if np.linalg.norm(cgf.p - nuc1) < 1e-3:
                idx1.append(i)
            if np.linalg.norm(cgf.p - nuc2) < 1e-3:
                idx2.append(i)

        return idx1, idx2
