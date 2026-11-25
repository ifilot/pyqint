import unittest
from pyqint import HF, Molecule, MoleculeBuilder, PyQInt
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
            np.testing.assert_almost_equal(D[:,:,i], exact, decimal=6)

        batch = integrator.build_dipole_openmp(iter(cgfs))
        np.testing.assert_allclose(np.triu(batch), D.transpose(2, 0, 1), atol=1e-14)
        np.testing.assert_allclose(batch, batch.transpose(0, 2, 1), atol=0)

    def test_neutral_molecular_dipole(self):
        solver = HF(MoleculeBuilder.from_name("h2o"), "sto3g", eri_storage="packed")
        reference = solver.rhf()
        origin = [0.25, -0.5, 1.5]
        result = solver.rhf(calc_dipole=True, dipole_origin=origin)
        unrestricted = solver.uhf(multiplicity=1, calc_dipole=True)

        self.assertNotIn("dipole", reference)
        self.assertEqual(result["energy"], reference["energy"])
        self.assertGreater(np.linalg.norm(result["dipole"]), 0.1)
        np.testing.assert_allclose(result["dipole"], unrestricted["dipole"], atol=2e-6)

    def test_charged_dipole_origin(self):
        mol = MoleculeBuilder.from_name("h2")
        mol.set_charge(1)
        solver = HF(mol, "sto3g")
        origin = np.array([0.5, 0.25, -0.125])
        centered = solver.uhf(multiplicity=2, calc_dipole=True)
        shifted = solver.uhf(multiplicity=2, calc_dipole=True, dipole_origin=origin)

        np.testing.assert_allclose(centered["dipole"], 0, atol=1e-12)
        np.testing.assert_allclose(shifted["dipole"], -2.54174647 * origin, atol=1e-8)

if __name__ == '__main__':
    unittest.main()
