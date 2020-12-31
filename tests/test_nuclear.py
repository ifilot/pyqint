import unittest
from pyqint.pyqint import PyQInt, cgf, gto
from copy import deepcopy
import numpy as np

class TestNuclear(unittest.TestCase):

    def test_gto_nuclear(self):
        # construct integrator object
        integrator = PyQInt()

        # test GTO
        gto1 = gto(0.154329, [0.0, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = gto(0.535328, [0.0, 0.0, 0.0], 0.623914, 0, 0, 0)
        gto3 = gto(0.444635, [0.0, 0.0, 0.0], 0.168855, 0, 0, 0)
        nuclear = integrator.nuclear_gto(gto1, gto1, [0.0, 0.0, 1.0])
        result = -0.31049036979675293
        np.testing.assert_almost_equal(nuclear, result, 8)

    def test_cgf_nuclear(self):
        integrator = PyQInt()

        # build cgf for hydrogen seperated by 1.4 a.u.
        cgf1 = cgf([0.0, 0.0, 0.0])

        cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
        cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
        cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

        cgf2 = deepcopy(cgf1)
        cgf2.p[2] = 1.4

        V1 = np.zeros((2,2))
        V1[0,0] = integrator.nuclear(cgf1, cgf1, cgf1.p, 1)
        V1[0,1] = V1[1,0] = integrator.nuclear(cgf1, cgf2, cgf1.p, 1)
        V1[1,1] = integrator.nuclear(cgf2, cgf2, cgf1.p, 1)

        V2 = np.zeros((2,2))
        V2[0,0] = integrator.nuclear(cgf1, cgf1, cgf2.p, 1)
        V2[0,1] = V2[1,0] = integrator.nuclear(cgf1, cgf2, cgf2.p, 1)
        V2[1,1] = integrator.nuclear(cgf2, cgf2, cgf2.p, 1)

        V11 = -1.2266135215759277
        V12 = -0.5974172949790955
        V22 = -0.6538270711898804
        np.testing.assert_almost_equal(V1[0,0], V11, 8)
        np.testing.assert_almost_equal(V1[1,1], V22, 8)
        np.testing.assert_almost_equal(V1[0,1], V12, 8)
        np.testing.assert_almost_equal(V2[0,0], V22, 8)
        np.testing.assert_almost_equal(V2[1,1], V11, 8)
        np.testing.assert_almost_equal(V2[0,1], V12, 8)

if __name__ == '__main__':
    unittest.main()
