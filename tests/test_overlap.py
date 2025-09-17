import unittest
from pyqint import PyQInt, GTO, Molecule
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
        gto1 = GTO(0.154329, [0.0, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = GTO(0.535328, [0.0, 0.0, 0.0], 0.623914, 0, 0, 0)
        gto3 = GTO(0.444635, [0.0, 0.0, 0.0], 0.168855, 0, 0, 0)
        overlap = integrator.overlap_gto(gto1, gto1)
        result = 0.31055691838264465
        np.testing.assert_almost_equal(overlap, result, 4)

    def test_cgf_overlap(self):
        """
        Test kinetic integrals for contracted Gaussians

        Sij = <cgf_i | cgf_j>
        """

        # construct integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule("H2")
        mol.add_atom('H', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        S = np.zeros((2,2))
        S[0,0] = integrator.overlap(cgfs[0], cgfs[0])
        S[0,1] = S[1,0] = integrator.overlap(cgfs[0], cgfs[1])
        S[1,1] = integrator.overlap(cgfs[1], cgfs[1])

        S11 = 1.0
        S12 = 0.65931845
        np.testing.assert_almost_equal(S[0,0], S11, 4)
        np.testing.assert_almost_equal(S[1,1], S11, 4)
        np.testing.assert_almost_equal(S[0,1], S12, 4)

if __name__ == '__main__':
    unittest.main()
