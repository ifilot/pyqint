import unittest
from pyqint import PyQInt, cgf, gto, Molecule
from copy import deepcopy
import numpy as np
import multiprocessing
import os

class TestCGF(unittest.TestCase):

    def testFunctionsCGF(self):
        """
        Test functions embedded in CGF class
        """

        # construct integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule()
        mol.add_atom('H', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')

        # test values at these coordinates
        coords = []
        for x in np.linspace(0, 10, 10):
            coords.append([x, x, x])
        coords = np.array(coords)
        amps = cgfs[0].get_amps(coords)

        ans = [6.2825e-01, 7.1229e-02, 6.8672e-03, 3.0000e-04, 3.7662e-06,
               1.3536e-08, 1.3927e-11, 4.1021e-15, 3.4591e-19, 8.3505e-24]

        np.testing.assert_almost_equal(amps, ans, 4)

if __name__ == '__main__':
    unittest.main()
