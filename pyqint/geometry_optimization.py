# -*- coding: utf-8 -*-

from .pyqint import HF
import numpy as np

class GeometryOptimization:
    """
    Class for performing Crystal Orbital Hamilton Population Analysis
    """
    def __init__(self, mol, basis):
        self.mol = mol
        self.basis = basis

    def run(self):
        """
        Perform geometry optimization
        """
        # perform initial HF calculation
        HF().rhf(self.mol, self.basis, calc_forces=True)


