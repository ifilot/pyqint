import unittest
from pyqint import PyQInt, cgf, gto, Molecule, HF
from copy import deepcopy
import numpy as np
import multiprocessing
import os
from nose.tools import nottest

class TestBasisSet(unittest.TestCase):
    def testOrder(self):
        results1 = build_mol_order1()
        results2 = build_mol_order2()
        
        N = len(results1['overlap'])
        P = np.zeros((N,N))
        pp = [0,1,2,3,4,10,11,5,6,7,8,9,12,13]
        for i,p in enumerate(pp):
            P[i,p] = 1
        print(P)
        
        for l in ['overlap', 'kinetic', 'nuclear', 'transform']:
            np.testing.assert_almost_equal(results1[l], 
                                           P.transpose().dot(results2[l].dot(P)), 4)
   
def build_mol_order1():
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

    results = HF().rhf(mol, basis='sto3g')
    
    return results
    
def build_mol_order2():
    """
    Test Hartree-Fock calculation for Ethylene
    
    Data is compared to results obtained from Gaussian
    """
    mol = Molecule()
    mol.add_atom('C', -0.6530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
    mol.add_atom('H', -1.2288875372, -0.9156191261 ,0.0000000000, unit='angstrom')
    mol.add_atom('H', -1.2288875372,  0.9156191261 ,0.0000000000, unit='angstrom')
    mol.add_atom('C',  0.6530176758,  0.0000000000 ,0.0000000000, unit='angstrom')
    mol.add_atom('H',  1.2288875372,  0.9156191261 ,0.0000000000, unit='angstrom')
    mol.add_atom('H',  1.2288875372, -0.9156191261 ,0.0000000000, unit='angstrom')

    results = HF().rhf(mol, basis='sto3g')
    
    return results

if __name__ == '__main__':
    unittest.main()
