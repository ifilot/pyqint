import unittest
from pyqint.pyqint import PyQInt, cgf, gto
from copy import deepcopy
import numpy as np
import multiprocessing
import os

class TestIntegrals(unittest.TestCase):

    def test_integrals_h2(self):
        """
        Test automatic integral evaluation for hydrogen molecule
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
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = integrator.build_integrals(cgfs, nuclei, npar=ncpu, verbose=False)

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

    def test_integrals_h2o(self):
        """
        Test automatic integral evaluation for water molecule
        """

        # construct integrator object
        integrator = PyQInt()

        # define sto-3g coefficient for hydrogen
        coeff_hydrogen = {
            "1s": [[  0.154329, 3.425251],
                   [  0.535328, 0.623914],
                   [  0.444635, 0.168855]]
        }

        # define sto-3g coefficients for oxygen
        coeff_oxygen = {
            "1s": [[ 0.154329,130.709320],
                   [ 0.535328, 23.808861],
                   [ 0.444635,  6.443608]],
            "2s": [[-0.099967,  5.033151],
                   [ 0.399513,  1.169596],
                   [ 0.700115,  0.380389]],
            "2p": [[ 0.155916,  5.033151],
                   [ 0.607684,  1.169596],
                   [ 0.391957,  0.380389]],
        }

        # build cgfs
        cgfs = []
        h_pos = [[0.7570, 0.5860, 0.0000], [-0.757, 0.5860, 0.0000]] # positions of the two hydrogen atoms
        o_pos = [0,0,0] # position of the oxygen atom
        for pos in h_pos:
            for key, value in coeff_hydrogen.items():
                if key.endswith("s"):
                    cgfs.append(cgf(pos))
                    for coeff, alpha in value:
                        cgfs[-1].add_gto(coeff, alpha, 0, 0, 0)

        for key, value in coeff_oxygen.items():
            if key.endswith("s"):
                cgfs.append(cgf(o_pos))
                for coeff, alpha in value:
                    cgfs[-1].add_gto(coeff, alpha, 0, 0, 0)
            if key.endswith("p"):
                cgfs.append(cgf(o_pos)) # x
                for coeff, alpha in value:
                    cgfs[-1].add_gto(coeff, alpha, 1, 0, 0)
                cgfs.append(cgf(o_pos)) # y
                for coeff, alpha in value:
                    cgfs[-1].add_gto(coeff, alpha, 0, 1, 0)
                cgfs.append(cgf(o_pos)) # z
                for coeff, alpha in value:
                    cgfs[-1].add_gto(coeff, alpha, 0, 0, 1)

        # position and charge of nuclei
        nuclei = [[h_pos[0], 1], [h_pos[1], 1], [o_pos, 8]]

        # evaluate all integrals
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = integrator.build_integrals(cgfs, nuclei, npar=ncpu, verbose=False)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint.npy'))

        # assessment
        self.assertEqual(len(cgfs), 7)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(teint, teint_result, 4)

if __name__ == '__main__':
    unittest.main()
