import unittest
from pyqint import PyQInt, cgf, gto, Molecule
from copy import deepcopy
import numpy as np

class TestRepulsionDeriv(unittest.TestCase):

    # def testDerivH2O(self):
    #     """
    #     Test Derivatives of dihydrogen
    #     """

    #     # build integrator object
    #     integrator = PyQInt()

    #     # build hydrogen molecule
    #     mol = Molecule()
    #     mol.add_atom('O', 0.0, 0.0, 0.0)
    #     mol.add_atom('H', 0.7570, 0.5860, 0.0)
    #     mol.add_atom('H', -0.7570, 0.5860, 0.0)
    #     cgfs, nuclei = mol.build_basis('sto3g')

    #     # calculate derivative towards H1 in the x-direction
    #     fx1 = integrator.repulsion_deriv(cgfs[2], cgfs[2], nuclei[1][0], 0) # px
    #     fx2 = integrator.repulsion_deriv(cgfs[2], cgfs[3], nuclei[1][0], 0) # py

    #     ans1 = calculate_force_finite_difference(mol, 1, 2, 2, 0)
    #     ans2 = calculate_force_finite_difference(mol, 1, 3, 3, 0)

    #     # assert that the repulsion of two CGFs that spawn from
    #     # the same nucleus will not change in energy due to a
    #     # change of the nucleus coordinates
    #     np.testing.assert_almost_equal(fx1, ans1, 4)
    #     np.testing.assert_almost_equal(fx2, ans2, 4)

    #     # assert that the cross-terms will change
    #     fx3 = integrator.repulsion_deriv(cgfs[2], cgfs[5], nuclei[1][0], 0)
    #     fx4 = integrator.repulsion_deriv(cgfs[2], cgfs[5], nuclei[1][0], 0)

    #     ans3 = calculate_force_finite_difference(mol, 1, 2, 5, 0)
    #     ans4 = calculate_force_finite_difference(mol, 1, 2, 5, 0)

    #     np.testing.assert_almost_equal(fx3, ans3, 4)
    #     self.assertFalse(fx3 == 0.0)
    #     np.testing.assert_almost_equal(fx4, ans4, 4)
    #     self.assertFalse(fx4 == 0.0)

    def testDerivH2(self):
        """
        Test Derivatives of dihydrogen
        """

        # build integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule('H2')
        mol.add_atom('H', -0.5, 0.0, 0.0)
        mol.add_atom('H',  0.5, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')

        # calculate derivative towards H1 in the x-direction
        fx1 = integrator.repulsion_deriv(cgfs[0], cgfs[0], cgfs[0], cgfs[0], nuclei[0][0], 0)
        fx2 = integrator.repulsion_deriv(cgfs[1], cgfs[1], cgfs[1], cgfs[1], nuclei[0][0], 0)

        # mol | nuc_id | cgf_id1 | cgf_id2 | cgf_id3 | cgf_id4 | coord
        ans1 = calculate_force_finite_difference(mol, 0, 0, 0, 0, 0, 0)
        ans2 = calculate_force_finite_difference(mol, 0, 1, 1, 1, 1, 0)

        # assert that the repulsion of two CGFs that spawn from
        # the same nucleus will not change in energy due to a
        # change of the nucleus coordinates
        np.testing.assert_almost_equal(fx1, ans1, 4)
        np.testing.assert_almost_equal(fx1, 0.0, 4)
        np.testing.assert_almost_equal(fx2, ans2, 4)
        np.testing.assert_almost_equal(fx2, 0.0, 4)

        # assert that the cross-terms will change
        fx3 = integrator.repulsion_deriv(cgfs[0], cgfs[0], cgfs[1], cgfs[1], nuclei[0][0], 0)
        fx4 = integrator.repulsion_deriv(cgfs[0], cgfs[0], cgfs[1], cgfs[1], nuclei[1][0], 0)
        fx5 = integrator.repulsion_deriv(cgfs[1], cgfs[0], cgfs[1], cgfs[0], nuclei[0][0], 0)
        fx6 = integrator.repulsion_deriv(cgfs[1], cgfs[0], cgfs[1], cgfs[0], nuclei[1][0], 0)

        # mol | nuc_id | cgf_id1 | cgf_id2 | cgf_id3 | cgf_id4 | coord
        ans3 = calculate_force_finite_difference(mol, 0, 0, 0, 1, 1, 0)
        ans4 = calculate_force_finite_difference(mol, 1, 0, 0, 1, 1, 0)
        ans5 = calculate_force_finite_difference(mol, 0, 1, 0, 1, 0, 0)
        ans6 = calculate_force_finite_difference(mol, 1, 1, 0, 1, 0, 0)

        np.testing.assert_almost_equal(fx3, ans3, 4)
        np.testing.assert_almost_equal(fx4, ans4, 4)
        np.testing.assert_almost_equal(fx5, ans5, 4)
        np.testing.assert_almost_equal(fx6, ans6, 4)

    def testDerivRepulsionGTO(self):
        # construct integrator object
        integrator = PyQInt()

        # build gtos
        gto1 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto3 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto4 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 0, 0, 0)

        # perform integration
        rep1 = integrator.repulsion_gto_deriv(gto1, gto2, gto3, gto4, 0)
        rep2 = integrator.repulsion_gto_deriv(gto2, gto1, gto3, gto4, 0)

        # perform finite difference
        ans1 = calculate_deriv_gto(gto1, gto2, gto3, gto4, 0)

        np.testing.assert_almost_equal(rep1, ans1, 4)

def calculate_force_finite_difference(mol, nuc_id, cgf_id1, cgf_id2, cgf_id3, cgf_id4, coord):
    # build integrator object
    integrator = PyQInt()

    # distance
    diff = 0.00001

    mol1 = deepcopy(mol)
    mol1.atoms[nuc_id][1][coord] -= diff / 2
    mol2 = deepcopy(mol)
    mol2.atoms[nuc_id][1][coord] += diff / 2

    # build hydrogen molecule
    cgfs1, nuclei = mol1.build_basis('sto3g')
    left = integrator.repulsion(cgfs1[cgf_id1], cgfs1[cgf_id2], cgfs1[cgf_id3], cgfs1[cgf_id4])
    cgfs2, nuclei = mol2.build_basis('sto3g')
    right = integrator.repulsion(cgfs2[cgf_id1], cgfs2[cgf_id2], cgfs2[cgf_id3], cgfs2[cgf_id4])

    return (right - left) / diff

def calculate_deriv_gto(gto1, gto2, gto3, gto4, coord):
    # build integrator object
    integrator = PyQInt()

    # distance
    diff = 0.00001
    p = np.zeros(3)
    p[coord] = diff

    gto1_new1 = gto(gto1.c, gto1.p - 0.5 * p, gto1.alpha, gto1.l, gto1.m, gto1.n)
    gto1_new2 = gto(gto1.c, gto1.p + 0.5 * p, gto1.alpha, gto1.l, gto1.m, gto1.n)

    # build hydrogen molecule
    left = integrator.repulsion_gto(gto1_new1, gto2, gto3, gto4)
    right = integrator.repulsion_gto(gto1_new2, gto2, gto3, gto4)

    return (right - left) / diff

if __name__ == '__main__':
    unittest.main()
