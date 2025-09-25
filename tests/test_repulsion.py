import unittest
from pyqint import PyQInt, GTO, Molecule
import numpy as np

class TestRepulsion(unittest.TestCase):

    def test_gto_repulsion(self):
        """
        Test two-electron integrals for primitive GTOs

        (ij|kl) = <gto_i gto_j | r_ij | gto_k gto_l>
        """

        # construct integrator object
        integrator = PyQInt()

        # test GTO
        gto1 = GTO(0.154329, [0.0, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = GTO(0.535328, [0.0, 0.0, 0.0], 0.623914, 0, 0, 0)
        gto3 = GTO(0.444635, [0.0, 0.0, 0.0], 0.168855, 0, 0, 0)
        repulsion = integrator.repulsion_gto(gto1, gto1, gto1, gto1)
        result = 0.20141123130697272
        np.testing.assert_almost_equal(repulsion, result, 4)

    def test_cgf_repulsion(self):
        """
        Test two-electron integrals for contracted Gaussians

        (ij|kl) = <cgf_i cgf_j | r_ij | cgf_k cgf_l>
        """

        # construct integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule("H2")
        mol.add_atom('H', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        T1111 = integrator.repulsion(cgfs[0], cgfs[0], cgfs[0], cgfs[0])
        T1122 = integrator.repulsion(cgfs[0], cgfs[0], cgfs[1], cgfs[1])
        T1112 = integrator.repulsion(cgfs[0], cgfs[0], cgfs[0], cgfs[1])
        T2121 = integrator.repulsion(cgfs[1], cgfs[0], cgfs[1], cgfs[0])
        T1222 = integrator.repulsion(cgfs[0], cgfs[1], cgfs[1], cgfs[1])
        T2211 = integrator.repulsion(cgfs[1], cgfs[1], cgfs[0], cgfs[0])

        np.testing.assert_almost_equal(T1111, 0.7746056914329529, 4)
        np.testing.assert_almost_equal(T1122, 0.5696758031845093, 4)
        np.testing.assert_almost_equal(T1112, 0.4441076656879812, 4)
        np.testing.assert_almost_equal(T2121, 0.2970285713672638, 4)

        # test similarity between two-electron integrals
        np.testing.assert_almost_equal(T1222, T1112, 4)
        np.testing.assert_almost_equal(T1122, T2211, 4)

    def test_two_electron_indices(self):
        """
        Test unique two-electron indices
        """
        integrator = PyQInt()

        np.testing.assert_almost_equal(integrator.teindex(1,1,2,1), integrator.teindex(1,1,1,2), 4)
        np.testing.assert_almost_equal(integrator.teindex(1,1,2,1), integrator.teindex(2,1,1,1), 4)
        np.testing.assert_almost_equal(integrator.teindex(1,2,1,1), integrator.teindex(2,1,1,1), 4)
        np.testing.assert_almost_equal(integrator.teindex(1,1,1,2), integrator.teindex(1,1,2,1), 4)
        self.assertNotEqual(integrator.teindex(1,1,1,1), integrator.teindex(1,1,2,1))
        self.assertNotEqual(integrator.teindex(1,1,2,1), integrator.teindex(1,1,2,2))

if __name__ == '__main__':
    unittest.main()
