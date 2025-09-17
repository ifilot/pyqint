import unittest
from pyqint import GTO
import numpy as np

class testGTO(unittest.TestCase):

    def testNormalizationGTO(self):
        """
        Test getting amplitude from CGF
        """

        gto_h = GTO(1.0, (0,0,0), 0.4166, 0, 0, 0)
        norm = gto_h.get_norm()
        np.testing.assert_almost_equal(norm, 0.36957240951430304, 4)

if __name__ == '__main__':
    unittest.main()
