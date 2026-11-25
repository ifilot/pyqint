import unittest

import numpy as np

from pyqint import CGF, HF, MoleculeBuilder, PyQInt, SALCBasis, build_salcs


class TestSALCs(unittest.TestCase):

    def test_h2_integrals_and_hartree_fock(self):
        mol = MoleculeBuilder.from_name("h2")
        aos, nuclei = mol.build_basis("sto3g")
        weights = np.array([[1., 1.], [1., -1.]])
        salcs = build_salcs(iter(aos), weights, labels=["sigma_g", "sigma_u"])

        self.assertIsInstance(salcs, SALCBasis)
        self.assertEqual(salcs.labels, ("sigma_g", "sigma_u"))
        np.testing.assert_array_equal(weights, [[1, 1], [1, -1]])
        integrator = PyQInt()
        source = integrator.build_integrals_openmp(aos, nuclei)
        transformed = integrator.build_integrals_openmp(salcs.cgfs, nuclei)
        b = salcs.transform
        for actual, original in zip(transformed[:3], source[:3]):
            np.testing.assert_allclose(actual, b @ original @ b.T, atol=1e-13)
        np.testing.assert_allclose(transformed[0], np.eye(2), atol=1e-13)
        eri = np.einsum("ai,bj,ck,dl,ijkl->abcd", b, b, b, b, source[3])
        np.testing.assert_allclose(transformed[3], eri, atol=1e-13)

        ao_hf = HF(mol, aos).rhf(calc_forces=True)
        salc_hf = HF(mol, salcs.cgfs, eri_storage="packed", eri_engine="auto").rhf(calc_forces=True)
        self.assertTrue(ao_hf["converged"] and salc_hf["converged"])
        self.assertAlmostEqual(salc_hf["energy"], ao_hf["energy"], places=12)
        np.testing.assert_allclose(salc_hf["forces"], ao_hf["forces"], atol=1e-11)

    def test_rectangular_transform_uses_normalized_source_cgfs(self):
        left = CGF([0., 0., -1.])
        left.add_gto(6., 1.2, 0, 0, 0)
        right = CGF([0., 0., 1.])
        right.add_gto(0.125, 0.7, 0, 0, 0)
        aos = [left, right]
        weights = np.array([[1., -5e-5]])
        result = build_salcs(aos, weights)

        self.assertIsNone(result.labels)
        self.assertEqual(result.transform.shape, (1, 2))
        integrator = PyQInt()
        overlap = np.array([[integrator.overlap(a, b) for b in aos] for a in aos])
        expected = weights / np.sqrt((weights @ overlap @ weights.T)[0, 0])
        np.testing.assert_allclose(result.transform, expected, atol=1e-15)
        point = [0.2, -0.3, 0.5]
        values = [cgf.get_amp(point) for cgf in aos]
        np.testing.assert_allclose(result.cgfs[0].get_amp(point), (expected @ values)[0], atol=1e-14)
        np.testing.assert_array_equal([gto.p for gto in result.cgfs[0].gtos], [left.p, right.p])

    def test_rejects_invalid_combinations(self):
        mol = MoleculeBuilder.from_name("h2")
        aos, _ = mol.build_basis("sto3g")
        with self.assertRaisesRegex(ValueError, "2 columns"):
            build_salcs(aos, [[1]])
        with self.assertRaisesRegex(ValueError, "real and finite"):
            build_salcs(aos, [[1, np.nan]])
        with self.assertRaisesRegex(ValueError, "real and finite"):
            build_salcs(aos, [[1, 1j]])
        with self.assertRaisesRegex(ValueError, "one string label"):
            build_salcs(aos, [[1, 1]], labels=["a", "b"])
        with self.assertRaisesRegex(ValueError, "row 0 .*all-zero"):
            build_salcs(aos, [[0, 0]], labels=["empty"])

        scaled = CGF(aos[0].p)
        for gto in aos[0].gtos:
            scaled.add_gto_with_position(2 * gto.c, gto.p, gto.alpha, gto.l, gto.m, gto.n)
        with self.assertRaisesRegex(ValueError, "row 0 .*cancelled.*cancels to zero"):
            build_salcs([aos[0], scaled], [[1, -1]], labels=["cancelled"])
