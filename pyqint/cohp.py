# -*- coding: utf-8 -*-

from .pyqint import PyQInt
import numpy as np

class COHP:
    """
    Class for performing Crystal Orbital Hamilton Population Analysis
    """
    def __init__(self, res):
        # copy objects from Hartree-Fock result dictionary
        self.orbc_canonical = res['orbc']
        self.orbe_canonical = res['orbe']
        self.nelec = res['nelec']
        self.nuclei = res['nuclei']
        self.H = res['fock']
        self.cgfs = res['cgfs']
        self.maxiter = 1000
        self.occ = [1 if i < self.nelec//2 else 0 for i in range(0, len(self.cgfs))]

    def run(self, C, n1, n2):
        """
        Run the COHP algorithm, take a coefficient matrix as input

        C: coefficient matrix
        n1: index of nucleus 1
        n2: index of nucleus 2
        """
        if n1 == n2:
            raise Exception('Cannot perform COHP for the same atoms')

        # figure out which basis functions below to which nucleus
        cgfs1 = []
        cgfs2 = []
        nuc1 = self.nuclei[n1][0]
        nuc2 = self.nuclei[n2][0]
        for i,cgf in enumerate(self.cgfs):
            if np.linalg.norm(cgf.p - nuc1) < 1e-3:
                cgfs1.append(i)
            if np.linalg.norm(cgf.p - nuc2) < 1e-3:
                cgfs2.append(i)

        N = len(C)
        cohp = np.zeros(N)
        for k in range(0, N): # loop over molecular orbitals
            for i in cgfs1: # loop over atom A
                for j in cgfs2: # loop over atom B
                    cohp[k] += 2.0 * self.H[i,j] * C[i,k] * C[j,k]

        return cohp
