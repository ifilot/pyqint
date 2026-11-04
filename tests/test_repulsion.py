import unittest
from pyqint import PyQInt, GTO, Molecule, MoleculeBuilder
import numpy as np

def cartesian(cgf):
    """Cartesian orders (l, m, n) of a CGF"""
    gto = cgf.gtos[0]
    return (gto.l, gto.m, gto.n)

def d_indices(cgfs):
    """Indices of the first dxx, dxy and dxz functions"""
    i = [cartesian(c) for c in cgfs].index((2, 0, 0))
    assert [cartesian(c) for c in cgfs[i:i+3]] == [(2, 0, 0), (1, 1, 0), (1, 0, 1)]
    return i, i + 1, i + 2

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
        repulsion = integrator.repulsion_gto(gto1, gto1, gto1, gto1)
        result = 0.20141123130697272
        np.testing.assert_almost_equal(repulsion, result, 8)

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
        cgfs, _ = mol.build_basis('sto3g')

        T1111 = integrator.repulsion(cgfs[0], cgfs[0], cgfs[0], cgfs[0])
        T1122 = integrator.repulsion(cgfs[0], cgfs[0], cgfs[1], cgfs[1])
        T1112 = integrator.repulsion(cgfs[0], cgfs[0], cgfs[0], cgfs[1])
        T2121 = integrator.repulsion(cgfs[1], cgfs[0], cgfs[1], cgfs[0])
        T1222 = integrator.repulsion(cgfs[0], cgfs[1], cgfs[1], cgfs[1])
        T2211 = integrator.repulsion(cgfs[1], cgfs[1], cgfs[0], cgfs[0])

        np.testing.assert_almost_equal(T1111, 0.7746059442114874, 8)
        np.testing.assert_almost_equal(T1122, 0.5696759264718783, 8)
        np.testing.assert_almost_equal(T1112, 0.444107658891185, 8)
        np.testing.assert_almost_equal(T2121, 0.29702854118104927, 8)

        # test similarity between two-electron integrals
        np.testing.assert_almost_equal(T1222, T1112, 8)
        np.testing.assert_almost_equal(T1122, T2211, 8)

    def test_cgf_repulsion_p(self):
        # construct integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = MoleculeBuilder.from_name('H2O')
        cgfs, _ = mol.build_basis('sto3g')

        px_px_py_py = integrator.repulsion(cgfs[2], cgfs[2], cgfs[3], cgfs[3])
        px_px_px_py = integrator.repulsion(cgfs[2], cgfs[2], cgfs[2], cgfs[3])
        px_py_px_py = integrator.repulsion(cgfs[2], cgfs[3], cgfs[2], cgfs[3])
        px_py_pz_s =  integrator.repulsion(cgfs[2], cgfs[3], cgfs[4], cgfs[0])

        np.testing.assert_almost_equal(px_px_py_py, 0.7852702009215765, 8)
        np.testing.assert_almost_equal(px_px_px_py, 0.00000000000000000, 8)
        np.testing.assert_almost_equal(px_py_px_py, 0.04744444436276907, 8)
        np.testing.assert_almost_equal(px_py_pz_s,  0.00000000000000000, 8)

    def test_cgf_repulsion_d(self):
        # construct integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule("Fe")
        mol.add_atom('Fe', 0.0, 0.0, 0.0)
        mol.add_atom('C', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        # Fe 3d (xx, xy, xz) functions and the Fe 2px, 2py and 3pz functions
        dxx, dxy, dxz = d_indices(cgfs)
        px2, py2, pz3 = 2, 3, 8
        self.assertEqual(cartesian(cgfs[pz3]), (0, 0, 1))

        T1 = integrator.repulsion(cgfs[dxx], cgfs[dxx], cgfs[dxx], cgfs[dxx])
        T2 = integrator.repulsion(cgfs[dxx], cgfs[dxx], cgfs[pz3], cgfs[pz3])
        T3 = integrator.repulsion(cgfs[dxy], cgfs[dxy], cgfs[px2], cgfs[px2])
        T4 = integrator.repulsion(cgfs[dxz], cgfs[dxz], cgfs[py2], cgfs[py2])

        # reference values verified against PySCF 2.14
        np.testing.assert_almost_equal(T1, 1.1397183799612702, 8)
        np.testing.assert_almost_equal(T2, 0.9716112348728463, 8)
        np.testing.assert_almost_equal(T3, 1.2524693412052625, 8)
        np.testing.assert_almost_equal(T4, 1.2197180065968953, 8)

    def test_cgf_repulsion_d_fallback(self):
        # construct integrator object without any caching
        integrator = PyQInt(0,0)

        # build hydrogen molecule
        mol = Molecule("Fe")
        mol.add_atom('Fe', 0.0, 0.0, 0.0)
        mol.add_atom('C', 0.0, 0.0, 1.4)
        cgfs, nuclei = mol.build_basis('sto3g')

        # Fe 3d (xx, xy, xz) functions and the Fe 2px, 2py and 3pz functions
        dxx, dxy, dxz = d_indices(cgfs)
        px2, py2, pz3 = 2, 3, 8
        self.assertEqual(cartesian(cgfs[pz3]), (0, 0, 1))

        T1 = integrator.repulsion(cgfs[dxx], cgfs[dxx], cgfs[dxx], cgfs[dxx])
        T2 = integrator.repulsion(cgfs[dxx], cgfs[dxx], cgfs[pz3], cgfs[pz3])
        T3 = integrator.repulsion(cgfs[dxy], cgfs[dxy], cgfs[px2], cgfs[px2])
        T4 = integrator.repulsion(cgfs[dxz], cgfs[dxz], cgfs[py2], cgfs[py2])

        # reference values verified against PySCF 2.14
        np.testing.assert_almost_equal(T1, 1.1397183799612702, 8)
        np.testing.assert_almost_equal(T2, 0.9716112348728463, 8)
        np.testing.assert_almost_equal(T3, 1.2524693412052625, 8)
        np.testing.assert_almost_equal(T4, 1.2197180065968953, 8)

    def test_two_electron_indices(self):
        """
        Test unique two-electron indices
        """
        integrator = PyQInt()

        np.testing.assert_almost_equal(integrator.teindex(1,1,2,1), integrator.teindex(1,1,1,2), 7)
        np.testing.assert_almost_equal(integrator.teindex(1,1,2,1), integrator.teindex(2,1,1,1), 7)
        np.testing.assert_almost_equal(integrator.teindex(1,2,1,1), integrator.teindex(2,1,1,1), 7)
        np.testing.assert_almost_equal(integrator.teindex(1,1,1,2), integrator.teindex(1,1,2,1), 7)
        self.assertNotEqual(integrator.teindex(1,1,1,1), integrator.teindex(1,1,2,1))
        self.assertNotEqual(integrator.teindex(1,1,2,1), integrator.teindex(1,1,2,2))

if __name__ == '__main__':
    unittest.main()
