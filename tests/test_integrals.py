import unittest
from pyqint import Evaluator, Molecule
import numpy as np
import multiprocessing

class TestIntegrals(unittest.TestCase):

    def test_integrals_h2(self):
        """
        Test automatic integral evaluation for hydrogen molecule
        """

        # construct integrator object
        evaluator = Evaluator()

        # build hydrogen molecule
        mol = Molecule()
        mol.add_atom('H', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        # evaluate all integrals
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = evaluator.build_integrals(cgfs, nuclei, npar=ncpu, verbose=False)

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
        evaluator = Evaluator()

        # build water molecule
        mol = Molecule()
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        mol.add_atom('O', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')

        # evaluate all integrals
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = evaluator.build_integrals(cgfs, nuclei, npar=ncpu, verbose=False)

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
