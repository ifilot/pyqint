# -*- coding: utf-8 -*-

from .pyqint import PyQInt
import numpy as np

class MOPA:
    """
    Class for performing Molecular Orbital Population Analysis
    """
    def __init__(self, res):
        # copy objects from Hartree-Fock result dictionary
        self.orbc = res['orbc']
        self.S = res['overlap']
        self.orbe = res['orbe']
        self.nelec = res['nelec']
        self.nuclei = res['nuclei']
        self.H = res['fock']
        self.cgfs = res['cgfs']
        self.maxiter = 1000
        self.occ = [1 if i < self.nelec//2 else 0 for i in range(0, len(self.cgfs))]

    def moop(self, n1, n2):
        """
        Run the MOOP algorithm, take a coefficient matrix as input

        n1: index of nucleus 1
        n2: index of nucleus 2
        """
        if n1 == n2:
            raise Exception('Cannot perform MOOP for the same atoms')

        # figure out which basis functions below to which nucleus
        C = self.orbc
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
        moop_coeff = np.zeros(N)
        for k in range(0, N): # loop over molecular orbitals
            for i in cgfs1: # loop over atom A
                for j in cgfs2: # loop over atom B
                    moop_coeff[k] += 2.0 * self.S[i,j] * C[i,k] * C[j,k]

        return moop_coeff
    
    def mohp(self, n1, n2):
        """
        Run the MOHP algorithm, take a coefficient matrix as input

        C: coefficient matrix
        n1: index of nucleus 1
        n2: index of nucleus 2
        """
        if n1 == n2:
            raise Exception('Cannot perform MOHP for the same atoms')

        # figure out which basis functions below to which nucleus
        C = self.orbc
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
        mohp_coeff = np.zeros(N)
        for k in range(0, N): # loop over molecular orbitals
            for i in cgfs1: # loop over atom A
                for j in cgfs2: # loop over atom B
                    mohp_coeff[k] += 2.0 * self.H[i,j] * C[i,k] * C[j,k]

        return mohp_coeff
