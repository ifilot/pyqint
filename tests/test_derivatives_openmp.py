import unittest
from pyqint import PyQInt, Molecule
import numpy as np

class TestDerivativesOpenMP(unittest.TestCase):

    def test_hartree_fock_forces_h2o(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.1, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        # calculate forces using analytical derivatives
        cgfs, nuclei = mol.build_basis('sto3g')
        N = len(cgfs)

        # build containers
        S = np.zeros((3,3,N,N))
        T = np.zeros_like(S)
        V = np.zeros_like(S)

        # build integrator object
        integrator = PyQInt()

        ntei = integrator.teindex(N-1,N-1,N-1,N-1)+1 # calculate number of two-electron integrals
        teints = np.zeros((3,3,ntei))

        for n,deriv_nucleus in enumerate(nuclei):
            for d in range(3):
                for i,cgf1 in enumerate(cgfs):
                    for j,cgf2 in enumerate(cgfs):

                        # derivative of overlap matrix
                        S[n,d,i,j] += integrator.overlap_deriv(cgf1, cgf2, deriv_nucleus[0], d)

                        # derivative of kinetic matrix
                        T[n,d,i,j] += integrator.kinetic_deriv(cgf1, cgf2, deriv_nucleus[0], d)

                        # derivative nuclear electron attraction
                        for nucleus in nuclei:
                            V[n,d,i,j] += integrator.nuclear_deriv(cgf1, cgf2, nucleus[0], nucleus[1], deriv_nucleus[0], d)

                        ij = i*(i+1)//2+j;

                        for k,cgf3 in enumerate(cgfs):
                            for l,cgf4 in enumerate(cgfs):
                                kl = k * (k+1)//2+l;
                                if ij <= kl:
                                    idx = integrator.teindex(i,j,k,l)
                                    teints[n,d,idx] = integrator.repulsion_deriv(cgf1, cgf2, cgf3, cgf4, deriv_nucleus[0], d)


        S2, T2, V2, teints2 = integrator.build_geometric_derivatives_openmp(cgfs, nuclei)

        np.testing.assert_almost_equal(S, S2)
        np.testing.assert_almost_equal(T, T2)
        np.testing.assert_almost_equal(V, V2)
        np.testing.assert_almost_equal(teints, teints2)
