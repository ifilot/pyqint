import pickle
import unittest
from pyqint import CGF, PyQInt, GTO, Molecule, MoleculeBuilder
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

    def test_multicenter_repulsion(self):
        """Standalone ERIs and derivatives preserve explicit primitive centers."""
        def atom_centered(x):
            cgf = CGF([x, 0.0, 0.0])
            cgf.add_gto(1.0, 1.0, 0, 0, 0)
            return cgf

        a = atom_centered(0.0)
        b = CGF([2.1, 0.0, 0.0])
        b.add_gto_with_position(1.0, [1.4, 0.0, 0.0], 0.7, 0, 0, 0)
        b.add_gto_with_position(1.0, [2.8, 0.0, 0.0], 0.7, 0, 0, 0)

        step = 1e-5
        for engine in ("hellsing", "rys", "auto"):
            with self.subTest(engine=engine):
                integrator = PyQInt(eri_engine=engine)
                tensor = integrator.build_tei_openmp([a, b])
                self.assertAlmostEqual(integrator.repulsion(a, a, b, b), tensor[0, 0, 1, 1], places=13)

                plus = integrator.build_tei_openmp([atom_centered(step), b])[0, 0, 1, 1]
                minus = integrator.build_tei_openmp([atom_centered(-step), b])[0, 0, 1, 1]
                derivative = integrator.repulsion_deriv(a, a, b, b, a.p, 0)
                self.assertAlmostEqual(derivative, (plus - minus) / (2.0 * step), places=9)

    def test_select_eri_engine(self):
        """All engines agree on the same water/STO-3G ERI tensor."""
        cgfs, nuclei = MoleculeBuilder.from_name("H2O").build_basis("sto3g")
        default = PyQInt()
        reference = default.build_tei_openmp(cgfs)
        for name in ("hellsing_direct", "huzinaga", "huzinaga_scalar", "rys"):
            with self.subTest(engine=name):
                integrator = PyQInt(eri_engine=name)
                result = integrator.build_tei_openmp(cgfs)
                np.testing.assert_allclose(result, reference, rtol=1e-12, atol=1e-14)

        rys = PyQInt(eri_engine="rys")
        _, _, _, result = rys.build_integrals_openmp(cgfs, nuclei)
        np.testing.assert_allclose(result, reference, rtol=1e-12, atol=1e-14)
        _, _, _, reference_deriv = default.build_geometric_derivatives_openmp(cgfs, nuclei)
        _, _, _, result_deriv = rys.build_geometric_derivatives_openmp(cgfs, nuclei)
        np.testing.assert_allclose(result_deriv, reference_deriv, rtol=1e-12, atol=1e-14)

    def test_auto_matches_hellsing_for_packaged_basis(self):
        """Automatic dispatch matches cached Hellsing with a Pople basis."""
        cgfs, nuclei = MoleculeBuilder.from_name("H2O").build_basis("p321")
        auto = PyQInt(lmax=0, eri_engine="auto")
        hellsing = PyQInt(eri_engine="hellsing")
        expected = hellsing.build_tei_openmp(cgfs)
        actual, profile = auto.build_tei_openmp(cgfs, profile=True)
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)
        np.testing.assert_array_equal(actual, auto.build_tei_openmp(cgfs))
        self.assertEqual(actual.dtype, np.float64)
        self.assertTrue(actual.flags.c_contiguous)
        for stage in ("pair_setup", "quartet_planning", "integral_evaluation", "tensor_assembly", "numpy_conversion"):
            self.assertGreaterEqual(profile[stage], 0.0)
        self.assertGreater(profile["prepared_hellsing_quartets"], 0)
        self.assertGreater(profile["rys_quartets"], 0)
        self.assertGreater(profile["rys_shell_batches"], 0)
        pair_count = len(cgfs) * (len(cgfs) + 1) // 2
        planned = profile["parity_zero_quartets"] + profile["prepared_hellsing_quartets"] + profile["rys_quartets"]
        self.assertEqual(planned, pair_count * (pair_count + 1) // 2)
        np.testing.assert_allclose(auto.build_integrals_openmp(cgfs, nuclei)[3], expected, rtol=1e-12, atol=1e-14)

        # Interleave Cartesian components to exercise canonical shell grouping.
        order = list(range(0, len(cgfs), 2)) + list(range(1, len(cgfs), 2))
        permuted = auto.build_tei_openmp([cgfs[i] for i in order])
        np.testing.assert_allclose(permuted, expected[np.ix_(order, order, order, order)], rtol=1e-12, atol=1e-14)

        deriv = auto.build_geometric_derivatives_openmp(cgfs, nuclei)
        reference = hellsing.build_geometric_derivatives_openmp(cgfs, nuclei)
        for result, target in zip(deriv, reference):
            np.testing.assert_array_equal(result, target)

        quartet = [cgfs[i].gtos[0] for i in (0, 1, 2, 3)]
        self.assertEqual(auto.repulsion_gto(*quartet), hellsing.repulsion_gto(*quartet))
        self.assertEqual(auto.repulsion_gto_deriv(*quartet, 0), hellsing.repulsion_gto_deriv(*quartet, 0))

    def test_auto_reflection_zeros_and_derivatives(self):
        """Reflection zeros apply to the tensor at an exact common plane."""
        a = CGF([0.0, 0.0, 0.4])
        a.add_gto(1.0, 0.8, 0, 0, 0)

        def pz(z):
            cgf = CGF([1.1, 0.3, z])
            cgf.add_gto(1.0, 0.7, 0, 0, 1)
            return cgf

        b = pz(0.4)
        auto = PyQInt(eri_engine="auto")
        tensor, profile = auto.build_tei_openmp([a, b], profile=True)
        self.assertGreater(profile["parity_zero_quartets"], 0)
        self.assertEqual(tensor[0, 0, 0, 1], 0.0)
        step = 1e-5
        plus = auto.build_tei_openmp([a, pz(0.4 + step)])[0, 0, 0, 1]
        minus = auto.build_tei_openmp([a, pz(0.4 - step)])[0, 0, 0, 1]
        reference = (plus - minus) / (2 * step)
        self.assertGreater(abs(reference), 1e-3)
        packed = auto.build_geometric_derivatives_openmp([a, b], [(a.p, 1), (b.p, 1)])[3]
        self.assertAlmostEqual(packed[1, 2, auto.teindex(0, 0, 0, 1)], reference, places=9)

    def test_auto_rys_shell_batches(self):
        """The Rys shell path preserves Cartesian integrals on distinct centers."""
        cgfs = []
        for position in ([0.7, -0.2, 0.1], [-0.4, 0.6, -0.3]):
            for powers in ((2, 0, 0), (0, 2, 0), (0, 0, 2)):
                cgf = CGF(position)
                cgf.add_gto(0.8, 0.6, *powers)
                cgf.add_gto(0.2, 1.3, *powers)
                cgfs.append(cgf)

        # Shell membership also requires matching contraction coefficients.
        different_radial_contraction = CGF(cgfs[0].p)
        different_radial_contraction.add_gto(0.7, 0.6, 2, 0, 0)
        different_radial_contraction.add_gto(0.3, 1.3, 2, 0, 0)
        cgfs.append(different_radial_contraction)

        auto = PyQInt(eri_engine="auto")
        actual, profile = auto.build_tei_openmp(cgfs, profile=True)
        self.assertGreater(profile["rys_shell_batches"], 0)
        self.assertGreater(profile["rys_quartets"], profile["rys_shell_batches"])
        for engine in ("hellsing", "rys"):
            with self.subTest(engine=engine):
                expected = PyQInt(eri_engine=engine).build_tei_openmp(cgfs)
                np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)

    def test_rys_primitive_quadrature_regimes(self):
        """Exercise one to four roots at coincident and separated product centers."""
        hellsing = PyQInt()
        rys = PyQInt(eri_engine="rys")
        for distance in (0.0, 1e-6, 2.0, 20.0):
            a, b = [0.0, 0.0, 0.0], [distance, 0.0, 0.0]
            s_a = GTO(1.0, a, 0.9, 0, 0, 0)
            s_b = GTO(1.0, b, 0.7, 0, 0, 0)
            for power_a, power_b in ((1, 0), (1, 1), (2, 2), (3, 3)):
                left = GTO(1.0, a, 0.9, power_a, 0, 0)
                right = GTO(1.0, b, 0.7, power_b, 0, 0)
                with self.subTest(distance=distance, powers=(power_a, power_b)):
                    quartet = (left, s_a, right, s_b)
                    expected = hellsing.repulsion_gto(*quartet)
                    actual = rys.repulsion_gto(*quartet)
                    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)

    def test_packed_tei_matches_dense(self):
        mol = MoleculeBuilder.from_name("h2o")
        cgfs, nuclei = mol.build_basis("p631")
        integrator = PyQInt(eri_engine="auto")
        dense = integrator.build_integrals_openmp(cgfs, nuclei)
        packed = integrator.build_integrals_openmp(iter(cgfs), nuclei, packed=True)
        tei, profile = integrator.build_tei_openmp(cgfs, packed=True, profile=True)

        for expected, actual in zip(dense[:3], packed[:3]):
            np.testing.assert_array_equal(actual, expected)
        np.testing.assert_array_equal(tei, packed[3])
        indices = [integrator.teindex(*quartet) for quartet in np.ndindex(dense[3].shape)]
        np.testing.assert_array_equal(tei[indices], dense[3].ravel())
        self.assertEqual(profile["tensor_assembly"], 0.0)
        self.assertEqual(profile["numpy_conversion"], 0.0)
        self.assertGreater(profile["parity_zero_quartets"], 0)

        n = 512
        pairs = n * (n + 1) // 2
        self.assertEqual(integrator.teindex(n - 1, n - 1, n - 1, n - 1), pairs * (pairs + 1) // 2 - 1)

    def test_eri_engine_pickle(self):
        """Preserve the chosen engine on round trips and read the old state format."""
        for engine in ("rys", "auto"):
            with self.subTest(engine=engine):
                integrator = PyQInt(2, 8, eri_engine=engine)
                restored = pickle.loads(pickle.dumps(integrator))
                self.assertEqual(restored.__getstate__(), (2, 8, engine))
        restored = PyQInt()
        restored.__setstate__((2, 8))
        self.assertEqual(restored.__getstate__(), PyQInt(2, 8).__getstate__())
        with self.assertRaisesRegex(ValueError, "Unknown ERI engine"):
            PyQInt(eri_engine="unknown")

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
