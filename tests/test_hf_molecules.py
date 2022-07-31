import unittest
from pyqint import PyQInt, cgf, gto, Molecule, HF
from copy import deepcopy
import numpy as np
import multiprocessing
import os
from nose.tools import nottest

class TestHFMolecules(unittest.TestCase):

    def testH2(self):
        """
        Test Hartree-Fock calculation for H2
        
        Data is compared to results obtained from Gaussian
        """
        mol = Molecule()
        mol.add_atom('H', 0.0000, 0.0000, 0.3561150187, unit='angstrom')
        mol.add_atom('H', 0.0000, 0.0000, -0.3561150187, unit='angstrom')

        results = perform_hf(mol)

        # check that energy matches
        np.testing.assert_almost_equal(results, -1.1175059, 5)
    
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

        results = perform_hf(mol)

        # check that energy matches
        np.testing.assert_almost_equal(results, -77.0739544, 5)    
        
    @nottest
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

        results = perform_hf(mol)

        # check that energy matches
        np.testing.assert_almost_equal(results, -227.8913603, 5)

def perform_hf(mol):
    results = HF().rhf(mol, 'sto3g')
    return results['energy']

if __name__ == '__main__':
    unittest.main()
