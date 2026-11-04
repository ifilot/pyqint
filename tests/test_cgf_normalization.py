import unittest
import pickle
from pyqint import PyQInt, CGF, GTO
import numpy as np

class TestCGFNormalization(unittest.TestCase):

    def setUp(self):
        """
        Build a deliberately unnormalized contraction of two s-type GTOs
        """
        self.a1, self.a2 = 1.0, 0.25
        self.cgf = CGF([0.0, 0.0, 0.0])
        self.cgf.add_gto(1.0, self.a1, 0, 0, 0)
        self.cgf.add_gto(1.0, self.a2, 0, 0, 0)

        # normalization constants of the primitives and their overlap
        self.n1 = (2.0 * self.a1 / np.pi)**0.75
        self.n2 = (2.0 * self.a2 / np.pi)**0.75
        s12 = (2.0 * np.sqrt(self.a1 * self.a2) / (self.a1 + self.a2))**1.5
        self.norm = 1.0 / np.sqrt(2.0 + 2.0 * s12)

    def test_overlap(self):
        """
        Test that the integrals treat the contraction as normalized
        """
        integrator = PyQInt()
        np.testing.assert_almost_equal(integrator.overlap(self.cgf, self.cgf), 1.0, 14)

    def test_amplitude(self):
        """
        Test that the amplitude includes the normalization of the contraction
        """
        for r in ([0.0, 0.0, 0.0], [0.3, -0.2, 0.5], [1.0, 2.0, -1.5]):
            r2 = np.dot(r, r)
            ans = self.norm * (self.n1 * np.exp(-self.a1 * r2) + self.n2 * np.exp(-self.a2 * r2))
            np.testing.assert_allclose(self.cgf.get_amp(r), ans, rtol=1e-14)
            np.testing.assert_allclose(self.cgf.get_amp_f(*r), ans, rtol=1e-14)

    def test_gradient(self):
        """
        Test that the gradient includes the normalization of the contraction
        """
        for r in ([0.3, -0.2, 0.5], [1.0, 2.0, -1.5]):
            r2 = np.dot(r, r)
            ans = -2.0 * self.norm * np.array(r) * (self.a1 * self.n1 * np.exp(-self.a1 * r2) +
                                                     self.a2 * self.n2 * np.exp(-self.a2 * r2))
            np.testing.assert_allclose(self.cgf.get_grad(r), ans, rtol=1e-14)

    def test_pickle(self):
        """
        Test that a pickled and restored CGF keeps its normalization
        """
        cgf = pickle.loads(pickle.dumps(self.cgf))
        r = [0.3, -0.2, 0.5]
        np.testing.assert_allclose(cgf.get_amp(r), self.cgf.get_amp(r), rtol=1e-15)

    def test_pickle_multicenter(self):
        """
        Test that a pickled and restored CGF keeps the positions of GTOs that
        are not located at the CGF position (e.g. symmetry-adapted functions)
        """
        cgf = CGF([0.0, 0.0, 0.0])
        cgf.add_gto_with_position(1.0, [0.0, 0.0, 0.7], 0.8, 0, 0, 0)
        cgf.add_gto_with_position(-1.0, [0.0, 0.0, -0.7], 0.8, 0, 0, 0)
        restored = pickle.loads(pickle.dumps(cgf))

        for r in ([0.3, -0.2, 0.5], [0.0, 0.1, -0.9]):
            np.testing.assert_allclose(restored.get_amp(r), cgf.get_amp(r), rtol=1e-14)

        # restored function retains its overlap with the original one
        np.testing.assert_almost_equal(PyQInt().overlap(cgf, restored), 1.0, 12)

    def test_scale_invariance(self):
        """
        Test that scaling all coefficients by the same factor yields the
        same, normalized basis function for integrals and amplitudes alike
        """
        scaled = CGF([0.0, 0.0, 0.0])
        scaled.add_gto(3.7, self.a1, 0, 0, 0)
        scaled.add_gto(3.7, self.a2, 0, 0, 0)

        integrator = PyQInt()
        np.testing.assert_almost_equal(integrator.overlap(scaled, scaled), 1.0, 14)
        np.testing.assert_almost_equal(integrator.overlap(scaled, self.cgf), 1.0, 14)

        r = [0.3, -0.2, 0.5]
        np.testing.assert_allclose(scaled.get_amp(r), self.cgf.get_amp(r), rtol=1e-14)
        np.testing.assert_allclose(scaled.get_grad(r), self.cgf.get_grad(r), rtol=1e-14)

        grid = np.array([r, [1.0, 2.0, -1.5]])
        np.testing.assert_allclose(integrator.plot_wavefunction(grid, np.array([1.0]), [scaled]),
                                   integrator.plot_wavefunction(grid, np.array([1.0]), [self.cgf]),
                                   rtol=1e-14)

    def test_gto_raw_convention(self):
        """
        Test the documented convention for individual GTOs: the GTO integrals
        are evaluated over the unnormalized primitive, the amplitude over the
        normalized one, and neither includes the expansion coefficient
        """
        integrator = PyQInt()
        p = [0.1, 0.2, 0.3]
        gto1 = GTO(1.0, p, self.a1, 0, 0, 0)
        gto2 = GTO(2.5, p, self.a1, 0, 0, 0)
        n = gto1.get_norm()
        np.testing.assert_allclose(n, self.n1, rtol=1e-14)

        # raw integrals; multiplying by the normalization constants yields
        # the integrals over the (normalized) single-primitive CGF
        cgf = CGF(p)
        cgf.add_gto(1.0, self.a1, 0, 0, 0)
        np.testing.assert_allclose(integrator.overlap_gto(gto1, gto1),
                                   (np.pi / (2.0 * self.a1))**1.5, rtol=1e-14)
        np.testing.assert_allclose(n**2 * integrator.overlap_gto(gto1, gto1),
                                   integrator.overlap(cgf, cgf), rtol=1e-14)
        np.testing.assert_allclose(n**2 * integrator.kinetic_gto(gto1, gto1),
                                   integrator.kinetic(cgf, cgf), rtol=1e-14)
        np.testing.assert_allclose(n**4 * integrator.repulsion_gto(gto1, gto1, gto1, gto1),
                                   integrator.repulsion(cgf, cgf, cgf, cgf), rtol=1e-12)

        # the expansion coefficient is ignored by integrals and amplitude
        np.testing.assert_allclose(integrator.overlap_gto(gto2, gto2),
                                   integrator.overlap_gto(gto1, gto1), rtol=1e-14)
        np.testing.assert_allclose(gto2.get_amp(*p), gto1.get_amp(*p), rtol=1e-14)
        np.testing.assert_allclose(gto1.get_amp(*p), n, rtol=1e-14)

if __name__ == '__main__':
    unittest.main()
