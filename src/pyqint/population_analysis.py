# -*- coding: utf-8 -*-

from . import PyQInt
import numpy as np

class PA:
    """
    Class for performing Mulliken and Löwdin Population Analysis
    """
    def __init__(self, res):
        # copy objects from Hartree-Fock result dictionary
        self.P = res['density']
        self.S = res['overlap']
        self.cgfs = res['cgfs']
        self.orbc = res['orbc']
        self.nuclei = res['nuclei']

    def mulliken(self, n):
        """
        Perform Mulliken Population Analysis

        n: index of nucleus
        """

        # figure out which basis functions belongs to the nucleus
        C = self.orbc
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
    
    def lowdin(self, n):
        """
        Perform Löwdin Population Analysis

        n: index of nucleus
        """

        # diagonalize S
        s, U = np.linalg.eigh(self.S)

        # construct transformation matrix X, using Löwdin orthogonalization
        X = U @ np.diag(np.sqrt(s)) @ U.transpose()

        # figure out which basis functions belongs to the nucleus
        C = self.orbc
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
