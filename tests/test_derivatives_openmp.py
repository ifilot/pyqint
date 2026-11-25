import unittest
from itertools import product
from pyqint import CGF, PyQInt, Molecule
import numpy as np

class TestDerivativesOpenMP(unittest.TestCase):

    def test_multicenter_salc_derivatives(self):
        """Differentiate all integrals while displacing one SALC center at a time."""
        positions = np.array([[-0.7, -0.3, 0.15], [0.9, 0.4, -0.2]])

        def build_basis(centers):
            cgfs = []
            for phase in (1.0, -1.0):
                salc = CGF()
                for center, weight in zip(centers, (1.0, phase)):
                    salc.add_gto_with_position(weight, center, 0.8, 0, 0, 0)
                    salc.add_gto_with_position(0.15 * weight, center, 0.4, 1, 0, 0)
                cgfs.append(salc)
            return cgfs, [(center, 1) for center in centers]

        step = 1e-5
        cgfs, nuclei = build_basis(positions)
        auto = PyQInt(eri_engine="auto")
        hellsing = PyQInt(eri_engine="hellsing")
        tensor, profile = auto.build_tei_openmp(cgfs, profile=True)
        np.testing.assert_allclose(tensor, hellsing.build_tei_openmp(cgfs), rtol=1e-12, atol=1e-14)
        self.assertEqual(profile["prepared_hellsing_quartets"], 0)
        self.assertEqual(profile["rys_quartets"], 0)
        reference = hellsing.build_geometric_derivatives_openmp(cgfs, nuclei)
        for result, target in zip(auto.build_geometric_derivatives_openmp(cgfs, nuclei), reference):
            np.testing.assert_array_equal(result, target)

        for engine in ("hellsing", "rys"):
            with self.subTest(engine=engine):
                integrator = PyQInt(eri_engine=engine)
                derivatives = integrator.build_geometric_derivatives_openmp(cgfs, nuclei)
                indices = np.array([integrator.teindex(*q) for q in product(range(2), repeat=4)])
                dense_eri = derivatives[3][..., indices.reshape((2, 2, 2, 2))]

                for nucleus, coord in product(range(2), range(3)):
                    with self.subTest(nucleus=nucleus, coord=coord):
                        shift = np.zeros_like(positions)
                        shift[nucleus, coord] = step
                        plus = integrator.build_integrals_openmp(*build_basis(positions + shift))
                        minus = integrator.build_integrals_openmp(*build_basis(positions - shift))
                        for deriv, pos, neg in zip((*derivatives[:3], dense_eri), plus, minus):
                            expected = (pos - neg) / (2.0 * step)
                            np.testing.assert_allclose(deriv[nucleus, coord], expected, atol=2e-9, rtol=1e-7)

                a, b = cgfs
                center = positions[0]
                standalone = (
                    integrator.overlap_deriv(a, b, center, 0),
                    integrator.kinetic_deriv(a, b, center, 0),
                    sum(integrator.nuclear_deriv(a, b, rc, charge, center, 0) for rc, charge in nuclei),
                    integrator.repulsion_deriv(a, a, b, b, center, 0),
                )
                expected = [d[0, 0, 0, 1] for d in derivatives[:3]]
                expected.append(dense_eri[0, 0, 0, 0, 1, 1])
                np.testing.assert_allclose(standalone, expected, atol=1e-12)

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
