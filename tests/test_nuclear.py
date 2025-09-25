import unittest
from pyqint import PyQInt, GTO, Molecule
import numpy as np

class TestNuclear(unittest.TestCase):

    def test_gto_nuclear(self):
        """
        Test nuclear attraction integral for GTOs

        V^{(c)}_ij = <gto_i | -1 / |r-Rc| | gto_j>
        """

        # construct integrator object
        integrator = PyQInt()

        # test GTO
        gto1 = GTO(0.154329, [0.0, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = GTO(0.535328, [0.0, 0.0, 0.0], 0.623914, 0, 0, 0)
        gto3 = GTO(0.444635, [0.0, 0.0, 0.0], 0.168855, 0, 0, 0)
        nuclear = integrator.nuclear_gto(gto1, gto1, [0.0, 0.0, 1.0])
        result = -0.31049036979675293
        np.testing.assert_almost_equal(nuclear, result, 4)

    def test_cgf_nuclear(self):
        """
        Test nuclear attraction integrals for contracted Gaussians

        V^{(c)}_ij = <cgf_i | -Zc / |r-Rc| | cgf_j>
        """

        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule("H2")
        mol.add_atom('H', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        V1 = np.zeros((2,2))
        V1[0,0] = integrator.nuclear(cgfs[0], cgfs[0], cgfs[0].p, 1)
        V1[0,1] = V1[1,0] = integrator.nuclear(cgfs[0], cgfs[1], cgfs[0].p, 1)
        V1[1,1] = integrator.nuclear(cgfs[1], cgfs[1], cgfs[0].p, 1)

        V2 = np.zeros((2,2))
        V2[0,0] = integrator.nuclear(cgfs[0], cgfs[0], cgfs[1].p, 1)
        V2[0,1] = V2[1,0] = integrator.nuclear(cgfs[0], cgfs[1], cgfs[1].p, 1)
        V2[1,1] = integrator.nuclear(cgfs[1], cgfs[1], cgfs[1].p, 1)

        V11 = -1.2266135215759277
        V12 = -0.5974172949790955
        V22 = -0.6538270711898804
        np.testing.assert_almost_equal(V1[0,0], V11, 4)
        np.testing.assert_almost_equal(V1[1,1], V22, 4)
        np.testing.assert_almost_equal(V1[0,1], V12, 4)
        np.testing.assert_almost_equal(V2[0,0], V22, 4)
        np.testing.assert_almost_equal(V2[1,1], V11, 4)
        np.testing.assert_almost_equal(V2[0,1], V12, 4)

if __name__ == '__main__':
    unittest.main()
