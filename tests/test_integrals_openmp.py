import unittest
from pyqint import PyQInt, Molecule
import numpy as np
import multiprocessing
import os
from nose.tools import nottest

class TestIntegralsOpenMP(unittest.TestCase):

    def test_integrals_h2o_openmp_sto3g(self):
        """
        Test automatic integral evaluation for water molecule
        """

        # construct integrator object
        integrator = PyQInt()

        # build water molecule
        mol = Molecule("H2O")
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        mol.add_atom('O', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')

        # evaluate all integrals
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_sto3g.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_sto3g.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_sto3g.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_sto3g.npy'))

        # assessment
        self.assertEqual(len(cgfs), 7)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(teint, teint_result, 4)

    def test_integrals_h2o_openmp_sto6g(self):
        """
        Test automatic integral evaluation for water molecule
        """

        # construct integrator object
        integrator = PyQInt()

        # build water molecule
        mol = Molecule("H2O")
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        mol.add_atom('O', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto6g')

        # evaluate all integrals
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_sto6g.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_sto6g.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_sto6g.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_sto6g.npy'))

        # assessment
        self.assertEqual(len(cgfs), 7)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(teint, teint_result, 4)

    def test_integrals_h2o_openmp_p321(self):
        """
        Test automatic integral evaluation for water molecule
        """

        # construct integrator object
        integrator = PyQInt()

        # build water molecule
        mol = Molecule("H2O")
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        mol.add_atom('O', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('p321')

        # evaluate all integrals
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_p321.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_p321.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_p321.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_p321.npy'))

        # assessment
        self.assertEqual(len(cgfs), 13)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(teint, teint_result, 4)

    def test_integrals_h2o_openmp_p631(self):
        """
        Test automatic integral evaluation for water molecule
        """

        # construct integrator object
        integrator = PyQInt()

        # build water molecule
        mol = Molecule("H2O")
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        mol.add_atom('O', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('p631')

        # evaluate all integrals
        ncpu = multiprocessing.cpu_count()
        S, T, V, teint = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_p631.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_p631.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_p631.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_p631.npy'))

        # assessment
        self.assertEqual(len(cgfs), 13)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(teint, teint_result, 4)

if __name__ == '__main__':
    unittest.main()
