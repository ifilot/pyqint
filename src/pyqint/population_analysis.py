# -*- coding: utf-8 -*-

from __future__ import annotations

from typing import Dict, Any, List

import numpy as np
import numpy.typing as npt

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

    The analysis operates on the output of a restricted Hartree–Fock
    calculation and assumes doubly occupied orbitals.
    """

    def __init__(self, res: Dict[str, Any]) -> None:
        """
        Parameters
        ----------
        res
            Result dictionary returned by a Hartree–Fock calculation.
        """
        # Molecular orbital coefficients (AO -> MO)
        self.orbc: Mat = res["orbc"]

        # Overlap matrix
        self.S: Mat = res["overlap"]

        # Orbital energies
        self.orbe: Vec = res["orbe"]

        # Density matrix
        self.P: Mat = res["density"]

        # Number of electrons (restricted, closed-shell assumed)
        self.nelec: int = res["nelec"]

        # Nuclear positions: [(position, charge), ...]
        self.nuclei = res["nuclei"]

        # Fock (Hamiltonian) matrix
        self.H: Mat = res["fock"]

        # Contracted Gaussian basis functions
        self.cgfs = res["cgfs"]

        # Occupation mask (2 electrons per occupied MO)
        nocc = self.nelec // 2
        self.occ: Vec = np.array(
            [1.0 if i < nocc else 0.0 for i in range(len(self.cgfs))]
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def mulliken(self, n: int):
        """
        Perform Mulliken Population Analysis

        n: index of nucleus
        """

        # figure out which basis functions belongs to the nucleus
        cgfs_n = []
        nuc = self.nuclei[n][0]
        for i,cgf in enumerate(self.cgfs):
            if np.linalg.norm(cgf.p - nuc) < 1e-3:
                cgfs_n.append(i)

        # detemrine the overlap weighted density matrix
        overlap_density = self.P @ self.S

        # add up electrons in basis functions localized on nucleus
        mulliken = 0
        for i in cgfs_n: # loop over atom cgfs
            mulliken += overlap_density[i,i]

        atomic_charge = self.nuclei[n][1] - mulliken

        return atomic_charge
    
    def lowdin(self, n: int):
        """
        Perform Löwdin Population Analysis

        n: index of nucleus
        """

        # diagonalize S
        s, U = np.linalg.eigh(self.S)

        # construct transformation matrix X, using Löwdin orthogonalization
        X = U @ np.diag(np.sqrt(s)) @ U.transpose()

        # figure out which basis functions belongs to the nucleus
        cgfs_n = []
        nuc = self.nuclei[n][0]
        for i,cgf in enumerate(self.cgfs):
            if np.linalg.norm(cgf.p - nuc) < 1e-3:
                cgfs_n.append(i)

        # determine P in (local) orthonormalized basis
        P_prime = X @ self.P @ X

        # add up electrons in basis functions localized on nucleus
        lowdin = 0
        for i in cgfs_n: # loop over atoms cgfs
            lowdin += P_prime[i,i]

        atomic_charge = self.nuclei[n][1] - lowdin

        return atomic_charge

    def moop(self, n1: int, n2: int) -> Vec:
        """
        Compute the Molecular Orbital Overlap Population (MOOP).

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.

        Returns
        -------
        ndarray
            MOOP values for each molecular orbital.
        """
        return self._population_analysis(n1, n2, matrix=self.S)

    def mohp(self, n1: int, n2: int) -> Vec:
        """
        Compute the Molecular Orbital Hamilton Population (MOHP).

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.

        Returns
        -------
        ndarray
            MOHP values for each molecular orbital.
        """
        return self._population_analysis(n1, n2, matrix=self.H)

    def mobi(self, n1: int, n2: int) -> Vec:
        """
        Compute the Molecular Orbital Bond Index (MOBI).

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.

        Returns
        -------
        ndarray
            MOBI values for each molecular orbital.
        """
        return self._population_analysis(n1, n2, matrix=self.P)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _population_analysis(self, n1: int, n2: int, matrix: Mat) -> Vec:
        """
        Shared implementation of MOHP/MOOP/MOBI.

        Parameters
        ----------
        n1, n2
            Indices of the two nuclei.
        matrix
            Either the overlap matrix (S), the Hamiltonian matrix (H), 
            or density matrix (P).

        Returns
        -------
        ndarray
            Population coefficients per molecular orbital.
        """
        if n1 == n2:
            raise ValueError(
                "Population analysis requires two distinct atoms."
            )

        # Determine basis functions belonging to each nucleus
        idx1, idx2 = self._basis_indices_for_atoms(n1, n2)

        C1 = self.orbc[idx1, :]
        C2 = self.orbc[idx2, :]
        M12 = matrix[np.ix_(idx1, idx2)]
        coeff = 2.0 * np.einsum('ik,ij,jk->k', C1, M12, C2, optimize=True)

        return coeff

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
