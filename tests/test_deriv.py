import unittest
from pyqint import PyQInt, cgf, gto, Molecule
from copy import deepcopy
import numpy as np
import multiprocessing
import os

class TestDeriv(unittest.TestCase):

    def testDerivBF_s1(self):
        """
        Test Derivative of the s-type basis functions in x-direction
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = [-0.5, 0.0, 0.0]
        gto1 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 0, 0, 0)

        # derivative towards bf coordinates
        t1 = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 0)

        # test the first derivative of the second hydrogen atom in the
        # Cartesian direction
        t2 = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 0)

        # also calculate this integral using finite difference
        fd_01 = calculate_fx_bf_finite_difference(gto1, nucleus)
        fd_02 = calculate_fx_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(-2 * t1, fd_01, 4)
        np.testing.assert_almost_equal(-2 * t2, fd_02, 4)

    def testDerivBF_s2(self):
        """
        Test Derivative of the s-type basis functions in y-direction
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = [-0.5, 0.0, 0.0]
        gto1 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 0, 0, 0)
        gto2 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 0, 0, 0)

        # derivative towards bf coordinates
        t1 = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 1)

        # test the first derivative of the second hydrogen atom in the
        # Cartesian direction
        t2 = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 1)

        # also calculate this integral using finite difference
        fd_01 = calculate_fy_bf_finite_difference(gto1, nucleus)
        fd_02 = calculate_fy_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(-2 * t1, fd_01, 4)
        np.testing.assert_almost_equal(-2 * t2, fd_02, 4)

    def testDerivBF_px(self):
        """
        Test Derivative of the px-type basis functions in x-direction
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = [-0.5, 0.0, 0.0]
        gto1 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 1, 0, 0)
        gto2 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 1, 0, 0)

        # derivative towards bf coordinates
        t1 = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 0)

        # test the first derivative of the second hydrogen atom in the
        # Cartesian direction
        t2 = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 0)

        # also calculate this integral using finite difference
        fx_fd_01 = calculate_fx_bf_finite_difference(gto1, nucleus)
        fx_fd_02 = calculate_fx_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(-2 * t1, fx_fd_01, 4)
        np.testing.assert_almost_equal(-2 * t2, fx_fd_02, 4)

    def testDerivBF_pxy(self):
        """
        Test Derivative of the px-type basis functions in y-direction
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = [-0.5, 0.0, 0.0]
        gto1 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 1, 0, 0)
        gto2 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 1, 0, 0)

        # derivative towards bf coordinates
        t1 = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 1)

        # test the first derivative of the second hydrogen atom in the
        # Cartesian direction
        t2 = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 1)

        # also calculate this integral using finite difference
        fx_fd_01 = calculate_fy_bf_finite_difference(gto1, nucleus)
        fx_fd_02 = calculate_fy_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(-2 * t1, fx_fd_01, 4)
        np.testing.assert_almost_equal(-2 * t2, fx_fd_02, 4)

    def testDerivBF_pyy(self):
        """
        Test Derivative of the py-type basis functions in y-direction
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = [-0.5, 0.0, 0.0]
        gto1 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 0, 1, 0)
        gto2 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 0, 1, 0)

        # derivative towards bf coordinates
        t1 = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 1)

        # test the first derivative of the second hydrogen atom in the
        # Cartesian direction
        t2 = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 1)

        # also calculate this integral using finite difference
        fx_fd_01 = calculate_fy_bf_finite_difference(gto1, nucleus)
        fx_fd_02 = calculate_fy_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(-2 * t1, fx_fd_01, 4)
        np.testing.assert_almost_equal(-2 * t2, fx_fd_02, 4)

    def testDerivBF_d(self):
        """
        Test Derivative of the dx2-type basis functions in x-direction
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = [-0.5, 0.0, 0.0]
        gto1 = gto(0.154329, [-0.50, 0.0, 0.0], 3.425251, 2, 0, 0)
        gto2 = gto(0.154329, [ 0.50, 0.0, 0.0], 3.425251, 2, 0, 0)

        # derivative towards bf coordinates
        t1 = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 0)

        # test the first derivative of the second hydrogen atom in the
        # Cartesian direction
        t2 = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 0)

        # also calculate this integral using finite difference
        fx_fd_01 = calculate_fx_bf_finite_difference(gto1, nucleus)
        fx_fd_02 = calculate_fx_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(-2 * t1, fx_fd_01, 4)
        np.testing.assert_almost_equal(-2 * t2, fx_fd_02, 4)

    def testDerivOpt_s1(self):
        """
        Test Derivative of the nuclear attraction operator in x-direction
        for s-type orbitals
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = np.array([-0.5, 0.0, 0.0])
        gto1 = gto(0.154329,  nucleus, 3.425251, 0, 0, 0)
        gto2 = gto(0.154329, -nucleus, 3.425251, 0, 0, 0)

        t1a = integrator.nuclear_gto_deriv_op(gto1, gto1, nucleus, 0)
        t1b = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 0)

        t2a = integrator.nuclear_gto_deriv_op(gto2, gto2, nucleus, 0)
        t2b = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 0)

        # also calculate this integral using finite difference
        fd_01 = calculate_fx_op_finite_difference(gto1, nucleus)
        fd_02a = calculate_fx_op_finite_difference(gto2, nucleus)
        fd_02b = calculate_fx_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(t1a, 0.0, 4)
        np.testing.assert_almost_equal(t1b, 0.0, 4)

        np.testing.assert_almost_equal(-2.0 * t2b, fd_02b, 4)
        np.testing.assert_almost_equal(t2a, fd_02a, 4)

    def testDerivOpt_s2(self):
        """
        Test Derivative of the nuclear attraction operator in y-direction
        for s-type orbitals
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = np.array([-0.5, 0.0, 0.0])
        gto1 = gto(0.154329,  nucleus, 3.425251, 0, 0, 0)
        gto2 = gto(0.154329, -nucleus, 3.425251, 0, 0, 0)

        t1a = integrator.nuclear_gto_deriv_op(gto1, gto1, nucleus, 1)
        t1b = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 1)

        t2a = integrator.nuclear_gto_deriv_op(gto2, gto2, nucleus, 1)
        t2b = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 1)

        # also calculate this integral using finite difference
        fd_01 = calculate_fy_op_finite_difference(gto1, nucleus)
        fd_02a = calculate_fy_op_finite_difference(gto2, nucleus)
        fd_02b = calculate_fy_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(t1a, 0.0, 4)
        np.testing.assert_almost_equal(t1b, 0.0, 4)

        np.testing.assert_almost_equal(-2.0 * t2b, fd_02b, 4)
        np.testing.assert_almost_equal(t2a, fd_02a, 4)

    def testDerivOpt_p1(self):
        """
        Test Derivative of the nuclear attraction operator in x-direction
        for p-type orbitals
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = np.array([-0.5, 0.0, 0.0])
        gto1 = gto(0.154329,  nucleus, 3.425251, 1, 0, 0)
        gto2 = gto(0.154329, -nucleus, 3.425251, 1, 0, 0)

        t1a = integrator.nuclear_gto_deriv_op(gto1, gto1, nucleus, 0)
        t1b = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 0)

        t2a = integrator.nuclear_gto_deriv_op(gto2, gto2, nucleus, 0)
        t2b = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 0)

        # also calculate this integral using finite difference
        fd_01 = calculate_fx_op_finite_difference(gto1, nucleus)
        fd_02a = calculate_fx_op_finite_difference(gto2, nucleus)
        fd_02b = calculate_fx_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(t1a, 0.0, 4)
        np.testing.assert_almost_equal(t1b, 0.0, 4)

        np.testing.assert_almost_equal(-2.0 * t2b, fd_02b, 4)
        np.testing.assert_almost_equal(t2a, fd_02a, 4)

    def testDerivOpt_d1(self):
        """
        Test Derivative of the nuclear attraction operator in x-direction
        for d-type orbitals
        """

        # construct integrator object
        integrator = PyQInt()

        # build gtos
        nucleus = np.array([-0.5, 0.0, 0.0])
        gto1 = gto(0.154329,  nucleus, 3.425251, 2, 0, 0)
        gto2 = gto(0.154329, -nucleus, 3.425251, 2, 0, 0)

        t1a = integrator.nuclear_gto_deriv_op(gto1, gto1, nucleus, 0)
        t1b = integrator.nuclear_gto_deriv_bf(gto1, gto1, nucleus, 0)

        t2a = integrator.nuclear_gto_deriv_op(gto2, gto2, nucleus, 0)
        t2b = integrator.nuclear_gto_deriv_bf(gto2, gto2, nucleus, 0)

        # also calculate this integral using finite difference
        fd_01 = calculate_fx_op_finite_difference(gto1, nucleus)
        fd_02a = calculate_fx_op_finite_difference(gto2, nucleus)
        fd_02b = calculate_fx_bf_finite_difference(gto2, nucleus)

        # testing
        np.testing.assert_almost_equal(t1a, 0.0, 4)
        np.testing.assert_almost_equal(t1b, 0.0, 4)

        np.testing.assert_almost_equal(-2.0 * t2b, fd_02b, 4)
        np.testing.assert_almost_equal(t2a, fd_02a, 4)

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

        fx1 = integrator.nuclear_deriv(cgfs[0], cgfs[0], nuclei[0][0], nuclei[0][1], 0)
        fx2 = integrator.nuclear_deriv(cgfs[1], cgfs[1], nuclei[0][0], nuclei[0][1], 0)

        ans1 = calculate_force_finite_difference(cgfs[0], cgfs[0], nuclei[0][0], nuclei[0][1], 0)
        ans2 = calculate_force_finite_difference(cgfs[1], cgfs[1], nuclei[0][0], nuclei[0][1], 0)

        np.testing.assert_almost_equal(fx1, ans1, 4)
        np.testing.assert_almost_equal(fx2, ans2, 4)


def calculate_force_finite_difference(cgf1, cgf2, nucleus, charge, coord):
    # build integrator object
    integrator = PyQInt()

    # distance
    diff = 0.00001

    # build hydrogen molecule
    nucleus[coord] -= diff / 2.0
    left = integrator.nuclear(cgf1, cgf2, nucleus, charge)
    nucleus[coord] += diff
    right = integrator.nuclear(cgf1, cgf2, nucleus, charge)

    return (right - left) / diff

def calculate_fx_bf_finite_difference(_gto, nucleus):
    """
    Calculate the basis function derivative in the x-direction using
    a finite difference method
    """
    # build integrator object
    integrator = PyQInt()

    diff = 0.01
    gto1 = gto(0.154329, [ _gto.p[0] - diff / 2, _gto.p[1], _gto.p[2]], 3.425251, _gto.l, _gto.m, _gto.n)
    gto2 = gto(0.154329, [ _gto.p[0] + diff / 2, _gto.p[1], _gto.p[2]], 3.425251, _gto.l, _gto.m, _gto.n)

    left = integrator.nuclear_gto(gto1, gto1, nucleus)
    right = integrator.nuclear_gto(gto2, gto2, nucleus)

    return (right - left) / diff

def calculate_fy_bf_finite_difference(_gto, nucleus):
    """
    Calculate the basis function derivative in the y-direction using
    a finite difference method
    """
    # build integrator object
    integrator = PyQInt()

    diff = 0.01
    gto1 = gto(0.154329, [ _gto.p[0], _gto.p[1] - diff / 2, _gto.p[2]], 3.425251, _gto.l, _gto.m, _gto.n)
    gto2 = gto(0.154329, [ _gto.p[0], _gto.p[1] + diff / 2, _gto.p[2]], 3.425251, _gto.l, _gto.m, _gto.n)

    left = integrator.nuclear_gto(gto1, gto1, nucleus)
    right = integrator.nuclear_gto(gto2, gto2, nucleus)

    return (right - left) / diff

def calculate_fx_op_finite_difference(_gto, nucleus):
    """
    Calculate the operator derivative in the x-direction using
    a finite difference method
    """
    # build integrator object
    integrator = PyQInt()

    diff = 0.000001
    nucleus[0] -= diff / 2.0
    left = integrator.nuclear_gto(_gto, _gto, nucleus)
    nucleus[0] += diff
    right = integrator.nuclear_gto(_gto, _gto, nucleus)

    return (right - left) / diff

def calculate_fy_op_finite_difference(_gto, nucleus):
    """
    Calculate the operator derivative in the y-direction using
    a finite difference method
    """
    # build integrator object
    integrator = PyQInt()

    diff = 0.000001
    nucleus[1] -= diff / 2.0
    left = integrator.nuclear_gto(_gto, _gto, nucleus)
    nucleus[1] += diff
    right = integrator.nuclear_gto(_gto, _gto, nucleus)

    return (right - left) / diff

if __name__ == '__main__':
    unittest.main()
