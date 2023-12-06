import unittest
from pyqint import PyQInt, Molecule
import numpy as np
import multiprocessing
import os

class TestIntegralsOpenMP(unittest.TestCase):

    def test_integrals_h2o_openmp_sto3g(self):
        """
        Test OpenMP integral evaluation for H2O using sto-3g basis set
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
        S, T, V, tetensor = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_sto3g.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_sto3g.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_sto3g.npy'))

        # grab unique two-electron integrals
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_sto3g.npy'))

        # repackage unique two-electron integrals into a tensor object
        tetensor_result = self.__repackage_tensor(teint_result, len(cgfs))

        # assessment
        self.assertEqual(len(cgfs), 7)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(tetensor, tetensor_result, 4)

    def test_integrals_h2o_openmp_sto6g(self):
        """
        Test OpenMP integral evaluation for H2O using sto-6g basis set
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
        S, T, V, tetensor = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_sto6g.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_sto6g.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_sto6g.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_sto6g.npy'))

        # repackage unique two-electron integrals into a tensor object
        tetensor_result = self.__repackage_tensor(teint_result, len(cgfs))

        # assessment
        self.assertEqual(len(cgfs), 7)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(tetensor, tetensor_result, 4)

    def test_integrals_h2o_openmp_p321(self):
        """
        Test OpenMP integral evaluation for H2O using p-321 basis set
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
        S, T, V, tetensor = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_p321.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_p321.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_p321.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_p321.npy'))

        # repackage unique two-electron integrals into a tensor object
        tetensor_result = self.__repackage_tensor(teint_result, len(cgfs))

        # assessment
        self.assertEqual(len(cgfs), 13)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(tetensor, tetensor_result, 4)

    def test_integrals_h2o_openmp_p631(self):
        """
        Test OpenMP integral evaluation for H2O using p-631 basis set
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
        S, T, V, tetensor = integrator.build_integrals_openmp(cgfs, nuclei)

        # load results from npy files
        S_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_overlap_p631.npy'))
        T_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_kinetic_p631.npy'))
        V_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_nuclear_p631.npy'))
        teint_result = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_teint_p631.npy'))

        # repackage unique two-electron integrals into a tensor object
        tetensor_result = self.__repackage_tensor(teint_result, len(cgfs))

        # assessment
        self.assertEqual(len(cgfs), 13)
        self.assertEqual(S.shape[0], S.shape[1])
        self.assertEqual(S.shape[0], len(cgfs))
        np.testing.assert_almost_equal(S, S_result, 4)
        np.testing.assert_almost_equal(T, T_result, 4)
        np.testing.assert_almost_equal(V, V_result, 4)
        np.testing.assert_almost_equal(tetensor, tetensor_result, 4)

    def __repackage_tensor(self, teint_result, N):
        # construct integrator object
        integrator = PyQInt()

        # repackage unique two-electron integrals into a tensor object
        tetensor_result = np.zeros((N,N,N,N))
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        idx = integrator.teindex(i,j,k,l)
                        tetensor_result[i,j,k,l] = teint_result[idx]

        return tetensor_result

if __name__ == '__main__':
    unittest.main()
