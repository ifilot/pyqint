import unittest
from pyqint import MoleculeBuilder
import numpy as np
import os

class TestMoleculeBuilder(unittest.TestCase):
    
    def test_construction(self):
        """
        Build a molecule from point group
        """
        mol = MoleculeBuilder().build_complex_td(1.0, 'C', 'H')

    def test_load_molecule_from_name(self):
        """
        Build a molecule from filename
        """
        mol = MoleculeBuilder().from_name('ch4')

        np.testing.assert_almost_equal(mol.get_atoms()[0][1], 
                                       np.array([0.0,0.0,0.0], dtype=np.float64))
        np.testing.assert_almost_equal(mol.get_atoms()[1][1], 
                                       np.array([0.6327670,0.6327670,0.6327670]) * 1.8897259886)
        np.testing.assert_equal(mol.get_atoms()[0][0], 'C')
        
        mol.build_basis('sto3g')
        np.testing.assert_equal(mol.get_nelec(), 10)
        
    def test_load_molecule_from_file(self):
        """
        Build a molecule from file
        """
        fname = os.path.join(os.path.dirname(__file__), 'results', 'ch4.xyz')
        mol = MoleculeBuilder().from_file(fname)

        np.testing.assert_almost_equal(mol.get_atoms()[0][1], 
                                       np.array([0.0,0.0,0.0], dtype=np.float64))
        np.testing.assert_almost_equal(mol.get_atoms()[1][1], 
                                       np.array([0.6327670,0.6327670,0.6327670]) * 1.8897259886)
        np.testing.assert_equal(mol.get_atoms()[0][0], 'C')
        
        mol.build_basis('sto3g')
        np.testing.assert_equal(mol.get_nelec(), 10)

if __name__ == '__main__':
    unittest.main()
