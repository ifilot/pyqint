import unittest
from pyqint.pyqint import PyQInt, cgf, gto
from copy import deepcopy
import numpy as np

class TestIntegrals(unittest.TestCase):

    def test_integrals(self):
        """
        Test automatic integral evaluation
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

        # build list of CGFs
        cgfs = []
        cgfs.append(cgf1)
        cgfs.append(cgf2)

        # position and charge of nuclei
        nuclei = [[cgf1.p, 1], [cgf2.p, 1]]

        # evaluate all integrals
        S, T, V, teint = integrator.build_integrals(cgfs, nuclei)

        # overlap integrals
        S_result = np.matrix([[1.0, 0.65931845],
                              [0.65931845, 1.0]])

        # kinetic integrals
        T_result = np.matrix([[0.7600315809249878, 0.2364544570446014],
                              [0.2364544570446014, 0.7600315809249878]])

        # nuclear attraction integrals
        V1_result = np.matrix([[-1.2266135215759277, -0.5974172949790955],
                               [-0.5974172949790955, -0.6538270711898804]])
        V2_result = np.matrix([[-0.6538270711898804, -0.5974172949790955],
                               [-0.5974172949790955, -1.2266135215759277]])
        V_result = V1_result + V2_result

        # two-electron integrals
        T1111 = 0.7746056914329529
        T1122 = 0.5696758031845093
        T1112 = 0.4441076219081878
        T1212 = 0.2970285713672638

        # perform tests
        np.testing.assert_almost_equal(S, S_result, 8)
        np.testing.assert_almost_equal(T, T_result, 8)
        np.testing.assert_almost_equal(V, V_result, 8)
        np.testing.assert_almost_equal(teint[integrator.teindex(0,0,0,0)], T1111, 8)
        np.testing.assert_almost_equal(teint[integrator.teindex(0,0,1,1)], T1122, 8)
        np.testing.assert_almost_equal(teint[integrator.teindex(0,0,0,1)], T1112, 8)
        np.testing.assert_almost_equal(teint[integrator.teindex(0,1,0,1)], T1212, 8)

if __name__ == '__main__':
    unittest.main()
