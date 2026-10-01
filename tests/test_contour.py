import os
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')

from pyqint import Molecule, HF, ContourPlotter

class TestContour(unittest.TestCase):

    def setUp(self):
        self.mol = Molecule()
        self.mol.add_atom('H', 0.0, 0.0, -0.7)
        self.mol.add_atom('H', 0.0, 0.0,  0.7)

    def test_contour_rhf(self):
        """
        Test building a contour plot for an RHF result
        """
        res = HF(self.mol, 'sto3g').rhf()
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'h2.png')
            ContourPlotter.build_contourplot(res, filename, 'xz', 3.0, 21, 1, 2)
            self.assertTrue(os.path.exists(filename))

        with self.assertRaises(ValueError):
            ContourPlotter.build_contourplot(res, 'h2.png', 'xz', 3.0, 21,
                                             1, 2, spin='alpha')

    def test_contour_uhf(self):
        """
        Test building contour plots for both spin channels of a UHF result
        """
        res = HF(self.mol, 'sto3g').uhf(multiplicity=3)
        with tempfile.TemporaryDirectory() as tmpdir:
            for spin in ('alpha', 'beta'):
                filename = os.path.join(tmpdir, 'h2_%s.png' % spin)
                ContourPlotter.build_contourplot(res, filename, 'xz', 3.0, 21,
                                                 1, 2, spin=spin)
                self.assertTrue(os.path.exists(filename))

        with self.assertRaises(ValueError):
            ContourPlotter.build_contourplot(res, 'h2.png', 'xz', 3.0, 21, 1, 2)

if __name__ == '__main__':
    unittest.main()
