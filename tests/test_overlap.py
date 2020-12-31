import unittest
from pyqint.pyqint import PyQInt, cgf, gto
from copy import deepcopy
import numpy as np

class TestOverlap(unittest.TestCase):

    def test_gto_overlap(self):
        """
        Test kinetic integrals for GTOs

        Sij = <gto_i | gto_j>
        """

        # construct integrator object
        integrator = PyQInt()

        # test GTO
        gto1 = gto(0.154329, [0.0, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = gto(0.535328, [0.0, 0.0, 0.0], 0.623914, 0, 0, 0)
        gto3 = gto(0.444635, [0.0, 0.0, 0.0], 0.168855, 0, 0, 0)
        overlap = integrator.overlap_gto(gto1, gto1)
        result = 0.31055691838264465
        np.testing.assert_almost_equal(overlap, result, 8)

    def test_cgf_overlap(self):
        """
        Test kinetic integrals for contracted Gaussians

        Sij = <cgf_i | cgf_j>
        """

        # construct integrator object
        integrator = PyQInt()

        # build cgf for hydrogen separated by 1.4 a.u.
        cgf1 = cgf([0.0, 0.0, 0.0])

        cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
        cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
        cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

        cgf2 = deepcopy(cgf1)
        cgf2.p[2] = 1.4

        S = np.zeros((2,2))
        S[0,0] = integrator.overlap(cgf1, cgf1)
        S[0,1] = S[1,0] = integrator.overlap(cgf1, cgf2)
        S[1,1] = integrator.overlap(cgf2, cgf2)

        S11 = 1.0
        S12 = 0.65931845
        np.testing.assert_almost_equal(S[0,0], S11, 8)
        np.testing.assert_almost_equal(S[1,1], S11, 8)
        np.testing.assert_almost_equal(S[0,1], S12, 8)

if __name__ == '__main__':
    unittest.main()
