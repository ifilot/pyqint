import unittest
from pyqint import PyQInt, GTO, Molecule
import numpy as np

class TestKinetic(unittest.TestCase):

    def test_gto_kinetic(self):
        """
        Test kinetic integrals for primitive GTOs

        Tij = <gto_I | -1/2 nabla^{2} | gto_j>
        """

        # construct integrator object
        integrator = PyQInt()

        # test GTO
        gto1 = GTO(0.154329, [0.0, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = GTO(0.535328, [0.0, 0.0, 0.0], 0.623914, 0, 0, 0)
        gto3 = GTO(0.444635, [0.0, 0.0, 0.0], 0.168855, 0, 0, 0)
        kinetic = integrator.kinetic_gto(gto1, gto1)
        result = 1.595603108406067
        np.testing.assert_almost_equal(kinetic, result, 4)

    def test_cgf_kinetic(self):
        """
        Test kinetic integrals for contracted Gaussians

        Tij = <cgf_i | -1/2 nabla^{2} | cgf_j>
        """

        # construct integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule("H2")
        mol.add_atom('H', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        T = np.zeros((2,2))
        T[0,0] = integrator.kinetic(cgfs[0], cgfs[0])
        T[0,1] = T[1,0] = integrator.kinetic(cgfs[0], cgfs[1])
        T[1,1] = integrator.kinetic(cgfs[1], cgfs[1])

        T11 = 0.7600315809249878
        T12 = 0.2364544570446014
        np.testing.assert_almost_equal(T[0,0], T11, 4)
        np.testing.assert_almost_equal(T[1,1], T11, 4)
        np.testing.assert_almost_equal(T[0,1], T12, 4)

if __name__ == '__main__':
    unittest.main()
