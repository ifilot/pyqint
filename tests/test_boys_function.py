import unittest
from pyqint import PyQInt, Molecule, HF, CGF, GTO
import numpy as np
from math import erf, pi, sqrt

def cartesian_powers(l):
    """
    Cartesian powers (i, j, k) with i + j + k = l, in the order xx...x first
    """
    return [(i, j, l - i - j) for i in range(l, -1, -1) for j in range(l - i, -1, -1)]

class TestBoysFunction(unittest.TestCase):

    def test_boys_function_order_zero(self):
        """
        Test F_0(T) through the nuclear attraction integral of an
        unnormalized s-type GTO with itself,

        <gto | 1 / |r-Rc| | gto> = (2 pi / p) F_0(p |Rc|^2), p = 2 alpha,

        over the range of T in which the Boys function is tabulated and beyond
        """
        integrator = PyQInt()
        alpha = 0.5
        p = 2.0 * alpha
        gto = GTO(1.0, [0.0, 0.0, 0.0], alpha, 0, 0, 0)

        for T in np.concatenate([[1e-8, 1e-6, 1e-3], np.linspace(0.01, 60.0, 337)]):
            rc = sqrt(T / p)
            exact = -2.0 * pi / p * 0.5 * sqrt(pi / T) * erf(sqrt(T))
            nuclear = integrator.nuclear_gto(gto, gto, [0.0, 0.0, rc])
            np.testing.assert_allclose(nuclear, exact, rtol=1e-13, atol=0)

    def test_repulsion_high_order(self):
        """
        Test two-electron integrals between g- and f-type functions, which
        require the Boys function up to order 14

        The reference values were obtained with PySCF 2.14.
        """
        # one primitive g shell (alpha = 0.15) and one primitive f shell
        # (alpha = 0.25), 0.8 bohr apart
        cgfs = []
        for p, l, alpha in (([0.0, 0.0, 0.0], 4, 0.15), ([0.0, 0.0, 0.8], 3, 0.25)):
            for (i, j, k) in cartesian_powers(l):
                cgf = CGF(p)
                cgf.add_gto(1.0, alpha, i, j, k)
                cgfs.append(cgf)

        # indices refer to g_xxxx ... g_zzzz (0-14) and f_xxx ... f_zzz (15-24)
        ans = {
            (14, 14, 24, 24): 0.30843614751021214,  # (g_zzzz g_zzzz | f_zzz f_zzz)
            (0, 0, 15, 15):   0.3042020411357217,   # (g_xxxx g_xxxx | f_xxx f_xxx)
            (14, 24, 14, 24): 0.11588532992276215,  # (g_zzzz f_zzz | g_zzzz f_zzz)
            (5, 18, 5, 18):   0.0076187458587338676,# (g_xxzz f_xyy | g_xxzz f_xyy)
            (2, 21, 2, 21):   0.0006440697028811094,# (g_xxxz f_yyy | g_xxxz f_yyy)
        }

        integrator = PyQInt()
        for idx, ref in ans.items():
            repulsion = integrator.repulsion(*[cgfs[i] for i in idx])
            np.testing.assert_allclose(repulsion, ref, rtol=1e-12, atol=1e-16)

    def test_diffuse_basis_bh3(self):
        """
        Test that the two-electron integrals of BH3 in the diffuse
        aug-cc-pVDZ basis set form a positive semidefinite matrix and that the
        Hartree-Fock energy matches PySCF (-26.07101412913098 Ht, using the
        same basis set and geometry)

        Inaccurate values of the Boys function make the ERI matrix
        indefinite, upon which the SCF procedure can collapse to an
        unphysical solution.
        """
        mol = Molecule()
        mol.add_atom('B', 0.0, 0.0, 0.0)
        for k in range(3):
            angle = 2.0 * np.pi * k / 3.0
            mol.add_atom('H', 1.60 * np.cos(angle), 1.60 * np.sin(angle), 0.0)

        res = HF(mol, 'aug-cc-pVDZ').rhf()
        N = len(res['cgfs'])
        eigenvalues = np.linalg.eigvalsh(res['tetensor'].reshape(N * N, N * N))

        self.assertGreater(eigenvalues.min(), -1e-10)
        self.assertTrue(res['converged'])
        np.testing.assert_almost_equal(res['energy'], -26.07101412913098, 8)

if __name__ == '__main__':
    unittest.main()
