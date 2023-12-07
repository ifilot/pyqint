import unittest
from pyqint import PyQInt, Molecule
import numpy as np
import os

class TestDipole(unittest.TestCase):

    def test_cgf_dipole(self):
        """
        Test dipole integrals for contracted Gaussians

        Dijk = <cgf_i | d_k | cgf_j>
        """

        # construct integrator object
        integrator = PyQInt()

        # build water molecule
        mol = Molecule("H2O")
        mol.add_atom('O',  0.00000, -0.07579, 0.0000, unit='angstrom')
        mol.add_atom('H',  0.86681,  0.60144, 0.0000, unit='angstrom')
        mol.add_atom('H', -0.86681,  0.60144, 0.0000, unit='angstrom')
        cgfs, nuclei = mol.build_basis('sto3g')

        N = len(cgfs)
        D = np.zeros((N,N,3))
        for i in range(N):
            for j in range(i,N):
                for k in range(0,3):
                    D[i,j,k] = integrator.dipole(cgfs[i], cgfs[j], k)

        # test dipole integrals
        for i,cart in enumerate(['x','y','z']):
            exact = np.load(os.path.join(os.path.dirname(__file__),
                                         'results',
                                         'h2o_dipole_%s.npy' % cart))
            np.testing.assert_almost_equal(D[:,:,i], exact, decimal=4)

if __name__ == '__main__':
    unittest.main()
