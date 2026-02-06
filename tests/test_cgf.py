import unittest
from pyqint import PyQInt, Molecule, CGF, MoleculeBuilder
import numpy as np
import os

class TestCGF(unittest.TestCase):
    def test_functions_cgf(self):
        """
        Test getting amplitude from CGF
        """

        # construct integrator object
        integrator = PyQInt()

        # get compile info
        compile_info = integrator.get_compile_info()

        # build hydrogen molecule
        mol = Molecule()
        mol.add_atom('H', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')

        # test values at these coordinates
        coords = []
        for x in np.linspace(0, 10, 10):
            coords.append([x, x, x])

        # results
        ans = np.array(
            [
                6.282468937e-01,
                7.122864446e-02,
                6.867199515e-03,
                3.000012405e-04,
                3.766188787e-06,
                1.353554025e-08,
                1.392653953e-11,
                4.102087793e-15,
                3.459081532e-19,
                8.350466908e-24,
            ]
        )

        # test for each coord
        amps = []
        for i,coord in enumerate(coords):
            amp = cgfs[0].get_amp(coord)
            amps.append(amp)

        np.testing.assert_almost_equal(amps, ans, 6)

    def test_plot_grid(self):
        """
        Test plotting of 1b2 molecular orbital of H2O
        """
        # coefficients
        coeff = [8.37612e-17, -2.73592e-16,  -0.713011, -1.8627e-17, 9.53496e-17, -0.379323,  0.379323]

        # construct integrator object
        integrator = PyQInt()

        # build water molecule
        mol = Molecule('H2O')
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')

        # build grid
        x = np.linspace(-2, 2, 50)
        y = np.linspace(-2, 2, 50)
        xx, yy = np.meshgrid(x,y)
        zz = np.zeros(len(x) * len(y))
        grid = np.vstack([xx.flatten(), yy.flatten(), zz]).reshape(3,-1).T
        res = integrator.plot_wavefunction(grid, coeff, cgfs).reshape((len(y), len(x)))

        ans = np.load(os.path.join(os.path.dirname(__file__), 'results', 'h2o_orb_1b2.npy'))
        np.testing.assert_almost_equal(res, ans, 6)

    def test_spherical_harmonic_on_site(self):
        """
        Check if the overlap matrix of all spherical harmonics up to l=6 is the identity matrix
        """
        basis_functions = []
        p0 = [0.0, 0.0, 0.0]
        for l in range(7):
            for m in range(-l, l+1):
                orb = CGF(p0)
                orb.add_spherical_gto(1.0, 1.0, l, m)
                basis_functions.append(orb)
        # build overlap matrix
        nbf = len(basis_functions)
        S = np.zeros((nbf, nbf))
        integrator = PyQInt()
        for i, orb1 in enumerate(basis_functions):
            for j, orb2 in enumerate(basis_functions):
                S[i,j] = integrator.overlap(orb1, orb2)
        np.testing.assert_almost_equal(S, np.eye(S.shape[0]), 6)

    def test_single_primitive_normalization(self):
        """
        Test that a single primitive GTO has self-overlap of 1.0
        """
        integrator = PyQInt()
        p0 = [0.0, 0.0, 0.0]

        # Test various angular momenta with different exponents
        test_cases = [
            # (alpha, l, m, n)
            (1.0, 0, 0, 0),   # s orbital
            (0.5, 1, 0, 0),   # px orbital
            (2.0, 0, 1, 0),   # py orbital
            (1.5, 0, 0, 1),   # pz orbital
            (1.0, 2, 0, 0),   # dxx orbital
            (1.0, 1, 1, 0),   # dxy orbital
            (0.8, 0, 0, 2),   # dzz orbital
        ]

        for alpha, l, m, n in test_cases:
            cgf = CGF(p0)
            cgf.add_gto(1.0, alpha, l, m, n)
            S = integrator.overlap(cgf, cgf)
            self.assertAlmostEqual(S, 1.0, places=10,
                msg=f"Single primitive ({l},{m},{n}) with alpha={alpha} failed: S={S}")

    def test_contracted_basis_normalization(self):
        """
        Test that contracted basis functions (STO-3G) have self-overlap of 1.0
        This verifies the contraction normalization is applied correctly.
        """
        integrator = PyQInt()

        # Build H2 molecule with STO-3G basis
        mol = Molecule()
        mol.add_atom('H', 0.0, 0.0, 0.0)
        mol.add_atom('H', 1.4, 0.0, 0.0)  # ~0.74 Angstrom bond length
        cgfs, _ = mol.build_basis('sto3g')

        # All diagonal elements of overlap matrix should be 1.0
        for i, cgf in enumerate(cgfs):
            S_ii = integrator.overlap(cgf, cgf)
            self.assertAlmostEqual(S_ii, 1.0, places=8,
                msg=f"Contracted basis function {i} has self-overlap {S_ii}")

    def test_overlap_matrix(self):
        """
        Test overlap matrix properties for a multi-atom system:
        1. Diagonal elements should be 1.0 (normalisation)
        2. Matrix should be symmetric
        3. Off-diagonal elements should be < 1.0
        """
        integrator = PyQInt()

        # Build water molecule with STO-3G basis
        mol = Molecule('H2O')
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        cgfs, _ = mol.build_basis('sto3g')

        # Build full overlap matrix
        nbf = len(cgfs)
        S = np.zeros((nbf, nbf))
        for i in range(nbf):
            for j in range(nbf):
                S[i, j] = integrator.overlap(cgfs[i], cgfs[j])

        # Test 1: Diagonal elements should be 1.0
        np.testing.assert_almost_equal(np.diag(S), np.ones(nbf), decimal=8,
            err_msg="Diagonal elements of overlap matrix are not 1.0")

        # Test 2: Matrix should be symmetric
        np.testing.assert_almost_equal(S, S.T, decimal=8,
            err_msg="Overlap matrix is not symmetric")

        # Test 3: Off-diagonal elements should have magnitude < 1.0
        off_diag = S - np.diag(np.diag(S))
        self.assertTrue(np.all(np.abs(off_diag) < 1.0),
            msg="Off-diagonal overlap elements should have magnitude < 1.0")
    
    def test_nh3_symmetrized_cgf(self):
        """
        Test whether symmetrized CGF is normalized
        """
        integrator = PyQInt()

        # Build water molecule with STO-3G basis
        mol = MoleculeBuilder.from_name('NH3')
        cgfs, _ = mol.build_basis('sto3g')
        self.assertEqual(len(cgfs), 8)

        # build symmetrized CGF
        cgf = CGF()
        for j,c in enumerate(np.array([0,0,0,0,0,1,1,-2], dtype=int)):
            if c != 0:
                for g in cgfs[j].gtos:
                    cgf.add_gto_with_position(
                        g.c * c,
                        g.p,
                        g.alpha,
                        g.l,
                        g.m,
                        g.n
                    )
        overlap = integrator.overlap(cgf, cgf)
        self.assertEqual(overlap, 1.0)

if __name__ == '__main__':
    unittest.main()
