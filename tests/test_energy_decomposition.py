import unittest
from pyqint import MoleculeBuilder, HF, Molecule
import numpy as np

class TestEnergyDecomposition(unittest.TestCase):

    def test_hartree_fock_h2o(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = MoleculeBuilder.from_name("h2o")

        res = HF(mol, 'sto3g').rhf()
        P = res['density']
        T = res['kinetic']
        V = res['nuclear']
        H = res['fock']
        enucrep = res['enucrep']
        energy = res['energy']
        
        np.testing.assert_almost_equal(0.5 * np.einsum('ji,ij', P, T+V+H) + enucrep, energy, decimal=16)

    def test_unrestricted_hartree_fock_ch3(self):
        """
        Test energy decomposition of an unrestricted Hartree-Fock calculation
        on the methyl radical
        """
        R = 2.039
        sqrt3 = np.sqrt(3.0)
        mol = Molecule()
        mol.add_atom('C', 0.0, 0.0, 0.0)
        mol.add_atom('H',  R, 0.0, 0.0)
        mol.add_atom('H', -0.5 * R,  0.5 * sqrt3 * R, 0.0)
        mol.add_atom('H', -0.5 * R, -0.5 * sqrt3 * R, 0.0)

        res = HF(mol, 'sto3g').uhf(multiplicity=2)
        Pa = res['density_alpha']
        Pb = res['density_beta']
        P = res['density']
        Hcore = res['hcore']
        Fa = res['fock_alpha']
        Fb = res['fock_beta']
        tetensor = res['tetensor']

        # the components sum to the total energy
        np.testing.assert_almost_equal(res['ecore'] + res['ej'] + res['ex_alpha'] +
                                       res['ex_beta'] + res['enucrep'],
                                       res['energy'], decimal=12)

        # total energy from the spin-resolved Fock matrices
        np.testing.assert_almost_equal(
            0.5 * np.einsum('ji,ij', Pa, Hcore + Fa) +
            0.5 * np.einsum('ji,ij', Pb, Hcore + Fb) + res['enucrep'],
            res['energy'], decimal=12)

        # individual components evaluated from the densities
        J = np.einsum('kl,ijlk->ij', P, tetensor)
        Ka = np.einsum('kl,iklj->ij', Pa, tetensor)
        Kb = np.einsum('kl,iklj->ij', Pb, tetensor)
        np.testing.assert_almost_equal(res['ecore'], np.einsum('ij,ji', P, Hcore), decimal=12)
        np.testing.assert_almost_equal(res['ej'], 0.5 * np.einsum('ij,ji', P, J), decimal=12)
        np.testing.assert_almost_equal(res['ex_alpha'], -0.5 * np.einsum('ij,ji', Pa, Ka), decimal=12)
        np.testing.assert_almost_equal(res['ex_beta'], -0.5 * np.einsum('ij,ji', Pb, Kb), decimal=12)

        # exchange stabilizes; the alpha channel holds the unpaired electron
        self.assertLess(res['ex_alpha'], res['ex_beta'])
        self.assertLess(res['ex_beta'], 0.0)

    def test_unrestricted_hartree_fock_closed_shell(self):
        """
        For a closed-shell system, the UHF energy components should match
        the RHF ones
        """
        mol = MoleculeBuilder.from_name("h2o")

        res_rhf = HF(mol, 'sto3g').rhf(tolerance=1e-12)
        res_uhf = HF(mol, 'sto3g').uhf(multiplicity=1, tolerance=1e-12)

        # the SCF convergence criterion acts on the total energy, which is
        # variational; individual components therefore converge less tightly
        np.testing.assert_almost_equal(res_uhf['ecore'], res_rhf['ecore'], decimal=5)
        np.testing.assert_almost_equal(res_uhf['ej'], res_rhf['erep'], decimal=5)
        np.testing.assert_almost_equal(res_uhf['ex_alpha'] + res_uhf['ex_beta'],
                                       res_rhf['ex'], decimal=5)
        np.testing.assert_almost_equal(res_uhf['ex_alpha'], res_uhf['ex_beta'], decimal=5)
        np.testing.assert_almost_equal(res_uhf['enucrep'], res_rhf['enucrep'], decimal=12)

if __name__ == '__main__':
    unittest.main()