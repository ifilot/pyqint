# -*- coding: utf-8 -*-

from . import HF
import numpy as np

class GeometryOptimization:
    """
    Class for performing Crystal Orbital Hamilton Population Analysis
    """
    def __init__(self, mol, basis, verbose=False):
        self.mol = mol
        self.basis = basis
        self.verbose = verbose

    def run(self):
        """
        Perform geometry optimization

        Currently, a very simple steepest descent algorithm is implemented, but
        this should of course be improved using something like conjugate
        gradient
        """
        # perform initial HF calculation
        res = HF().rhf(self.mol, self.basis, calc_forces=True, verbose=self.verbose)
        forces = res['forces']

        # steepest descent optimization parameter
        eta = 0.5

        # start geometry optimization procedure
        nriter = 0
        while np.max(np.linalg.norm(forces, axis=1) > 1e-4) and nriter < 100:
            for i,atom in enumerate(self.mol.atoms):
                atom[1] -= eta * forces[i,:]

            res = HF().rhf(self.mol, self.basis, calc_forces=True, verbose=self.verbose,
                           orbc_init=res['orbc'])
            forces = res['forces']
            nriter += 1

        return res
