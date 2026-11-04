import unittest
import pickle
from pyqint import PyQInt, CGF
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

if __name__ == '__main__':
    unittest.main()
