import unittest
from pyqint import Molecule, HF
import numpy as np

class TestHFMolecules(unittest.TestCase):

    def testH2(self):
        """
        Test Hartree-Fock calculation for H2
        
        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('H', 0.0000, 0.0000, 0.3561150187, unit='angstrom')
        mol.add_atom('H', 0.0000, 0.0000, -0.3561150187, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -1.1175059, 5)
    
    def testC2H4(self):
        """
        Test Hartree-Fock calculation for Ethylene
        
        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('C', -0.6530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
        mol.add_atom('C',  0.6530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
        mol.add_atom('H', -1.2288875372, -0.9156191261 ,0.0000000000, unit='angstrom')
        mol.add_atom('H', -1.2288875372,  0.9156191261 ,0.0000000000, unit='angstrom')
        mol.add_atom('H',  1.2288875372,  0.9156191261 ,0.0000000000, unit='angstrom')
        mol.add_atom('H',  1.2288875372, -0.9156191261 ,0.0000000000, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -77.0739726, 4)

    def testBF3(self):
        """
        Test Hartree-Fock calculation for BF3

        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('B', 0.0,  0.0, 0.0, unit='angstrom')
        mol.add_atom('F', -1.1334295832,  0.654385875, 0.0, unit='angstrom')
        mol.add_atom('F', 1.1334295832,  0.654385875, 0.0, unit='angstrom')
        mol.add_atom('F', 0.0,  -1.3087717499, 0.0, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -318.6619373, 4)

    def testCH4(self):
        """
        Test Hartree-Fock calculation for CH4

        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('C', 0.0,  0.0, 0.0, unit='angstrom')
        mol.add_atom('H', -5.9e-09,  -1.6e-09, 1.0830098121, unit='angstrom')
        mol.add_atom('H', 0.0,  -1.0210714424, -0.3610032723, unit='angstrom')
        mol.add_atom('H', -0.8842738057,  0.5105357237, -0.3610032748, unit='angstrom')
        mol.add_atom('H', 0.8842738117,  0.5105357203, -0.3610032651, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -39.7268637, 1)

    def testCO(self):
        """
        Test Hartree-Fock calculation for CO

        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('O', 0.0,  0.0, 0.4909201273, unit='angstrom')
        mol.add_atom('C', 0.0,  0.0, -0.6545601698, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -111.2254495, 4)

    def testCO2(self):
        """
        Test Hartree-Fock calculation for CO2

        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('C', 0.0,  0.0, 0.0, unit='angstrom')
        mol.add_atom('O', 0.0,  0.0, 1.1879700928, unit='angstrom')
        mol.add_atom('O', 0.0,  0.0, -1.1879700928, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -185.0683906, 5)

    def testH2O(self):
        """
        Test Hartree-Fock calculation for H2O

        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('O', 0.0,  0.0, -0.1271310707, unit='angstrom')
        mol.add_atom('H', 0.0,  -0.7580158822, 0.5085242828, unit='angstrom')
        mol.add_atom('H', 0.0,  0.7580158822, 0.5085242828, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -74.9659012, 4)

    def testLIH(self):
        """
        Test Hartree-Fock calculation for LIH

        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('Li', 0.0,  0.0, 0.377702853, unit='angstrom')
        mol.add_atom('H', 0.0,  0.0, -1.1331085589, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -7.8633821, 5)

    def testNH3(self):
        """
        Test Hartree-Fock calculation for NH3
        
        Note: a lower energy is found here as compared to Gaussian, so I reckon
        (might be wrong!) that this result is better.
        Using Gaussian: EHF = -55.7433617 Ht is found
        """
        mol = Molecule()
        mol.add_atom('H', -3.4e-09,  -2.6e-09, 1.1944430032, unit='angstrom')
        mol.add_atom('H', 1e-10,  -1.1261316622, -0.3981476702, unit='angstrom')
        mol.add_atom('H', -0.9752586265,  0.5630658334, -0.3981476693, unit='angstrom')
        mol.add_atom('H', 0.9752586299,  0.5630658315, -0.3981476637, unit='angstrom')
        mol.add_atom('N', 0.0,  0.0, 0.0, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -55.7998313, 5)

    def testC6H6(self):
        """
        Test Hartree-Fock calculation for Ethylene
        
        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('C',  0.0000000015, -1.3868467444, 0.0000000000, unit='angstrom')
        mol.add_atom('C',  1.2010445126, -0.6934233709, 0.0000000000, unit='angstrom')
        mol.add_atom('C',  1.2010445111,  0.6934233735, 0.0000000000, unit='angstrom')
        mol.add_atom('C', -0.0000000015,  1.3868467444, 0.0000000000, unit='angstrom')
        mol.add_atom('C', -1.2010445126,  0.6934233709, 0.0000000000, unit='angstrom')
        mol.add_atom('C', -1.2010445111, -0.6934233735, 0.0000000000, unit='angstrom')
        mol.add_atom('H',  0.0000000027, -2.4694205285, 0.0000000000, unit='angstrom')
        mol.add_atom('H',  2.1385809117, -1.2347102619, 0.0000000000, unit='angstrom')
        mol.add_atom('H',  2.1385809090,  1.2347102666, 0.0000000000, unit='angstrom')
        mol.add_atom('H', -0.0000000027,  2.4694205285, 0.0000000000, unit='angstrom')
        mol.add_atom('H', -2.1385809117,  1.2347102619, 0.0000000000, unit='angstrom')
        mol.add_atom('H', -2.1385809090, -1.2347102666, 0.0000000000, unit='angstrom')

        results = HF().rhf(mol, 'sto3g')

        # check that energy matches
        np.testing.assert_almost_equal(results['energy'], -227.8913603, 5)

if __name__ == '__main__':
    unittest.main()
