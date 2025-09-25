import unittest
from pyqint import Molecule, CGF, HF, MoleculeBuilder
import numpy as np

class TestCustomBasisSet(unittest.TestCase):
    
    def test_custom_basis_set_h2(self):
        mol = Molecule()
        mol.add_atom('H', 0.0000, 0.0000, 0.3561150187, unit='angstrom')
        mol.add_atom('H', 0.0000, 0.0000, -0.3561150187, unit='angstrom')        
        nuclei = mol.get_nuclei()

        cgfs = []
        for n in nuclei:
            _cgf = CGF(n[0])

            _cgf.add_gto(0.154329, 3.425251, 0, 0, 0)
            _cgf.add_gto(0.535328, 0.623914, 0, 0, 0)
            _cgf.add_gto(0.444635, 0.168855, 0, 0, 0)

            cgfs.append(_cgf)

        res = HF().rhf(mol, basis=cgfs)
        np.testing.assert_almost_equal(res['energy'], -1.1175059, 5)

    def test_custom_basis_set_co(self):
        mol = MoleculeBuilder().from_name('CO')
        cgfs, nuclei = mol.build_basis('sto3g')

        res = HF().rhf(mol, basis=cgfs)
        np.testing.assert_almost_equal(res['energy'], -111.2192571, 4)

if __name__ == '__main__':
    unittest.main()
