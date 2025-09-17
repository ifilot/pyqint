import unittest
from pyqint import CGF
import numpy as np
import itertools

class TestGrad(unittest.TestCase):
    """
    Test the calculation of analytical gradients by comparing these
    with the finite difference results
    """

    def test_cgf_grad(self):
        """
        Test gradient of basis functions
        """
        h = 1e-4 # set step size for finite difference

        pp = list(itertools.product(range(4), repeat=3))

        for exp in pp:
            for p in pp:
                l,m,n = exp
                cgft = CGF()
                cgft.add_gto(0.154329, 3.425251, l, m, n)
                cgft.add_gto(0.535328, 0.623914, l, m, n)
                cgft.add_gto(0.444635, 0.168855, l, m, n)

                grad = cgft.get_grad(p)
                np.testing.assert_almost_equal(grad, calculate_derivs_finite_diff(p, h, l,m,n), 4)
    
    def test_grad_density(self):
        """
        Test the gradient of the density from the gradient of the basis functions
        """
        h = 1e-4 # set step size for finite difference

        pp = list(itertools.product(range(4), repeat=3))
        for exp in pp:
            for p in pp:
                l,m,n = exp
                cgft = CGF()
                cgft.add_gto(0.154329, 3.425251, l, m, n)
                cgft.add_gto(0.535328, 0.623914, l, m, n)
                cgft.add_gto(0.444635, 0.168855, l, m, n)

                # calculate gradient of the density using chain rule
                grad = 2.0 * cgft.get_amp(p) * np.array(cgft.get_grad(p))
                np.testing.assert_almost_equal(grad, calculate_derivs_density_finite_diff(p, h, l,m,n), 4)

def calculate_derivs_finite_diff(p, h, l=0, m=0, n=0):
    """
    Determine the gradient using finite differences
    """
    grad = [0,0,0]
    for i in range(0,3):
        r = np.zeros(3)
        
        r[i] = -h
        cgft = CGF()
        cgft.add_gto(0.154329, 3.425251, l, m, n)
        cgft.add_gto(0.535328, 0.623914, l, m, n)
        cgft.add_gto(0.444635, 0.168855, l, m, n)
        vl = cgft.get_amp(p + r)
        r[i] = h
        vr = cgft.get_amp(p + r)
        grad[i] = (vr - vl) / (2. * h)
        
    return grad

def calculate_derivs_density_finite_diff(p, h, l=0, m=0, n=0):
    """
    Determine the gradient using finite differences
    """
    grad = [0,0,0]
    for i in range(0,3):
        r = np.zeros(3)
        
        r[i] = -h
        cgft = CGF()
        cgft.add_gto(0.154329, 3.425251, l, m, n)
        cgft.add_gto(0.535328, 0.623914, l, m, n)
        cgft.add_gto(0.444635, 0.168855, l, m, n)
        vl = cgft.get_amp(p + r)**2 # take square to get density
        r[i] = h
        vr = cgft.get_amp(p + r)**2 # take square to get density
        grad[i] = (vr - vl) / (2. * h)
        
    return grad

if __name__ == '__main__':
    unittest.main()