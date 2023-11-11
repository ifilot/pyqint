import os
from .molecule import Molecule 

class MoleculeBuilder:
    """
    Class that builds molecules from templates or from point group descriptions
    """
    def __init__(self):
        pass

    def from_name(self, molname):
        """
        Build molecule from molname
        """
        fname = os.path.join(os.path.dirname(__file__), 'molecules', molname.lower() + '.xyz')
        return self.from_file(fname)
       
    def from_file(self, path, molname=None):
        """
        Build molecule from file and return it
        """
        with open(path, 'r') as f:
            lines = f.readlines()
            
            nratoms = int(lines[0].strip())

            mol = Molecule(molname)
            for line in lines[2:2+nratoms]:
                pieces = line.split()
                mol.add_atom(pieces[0], float(pieces[1]), float(pieces[2]), float(pieces[3]), unit='angstrom')

            return mol

    def build_complex_td(self, r, at1, at2, unit='bohr'):
        mol = Molecule()
        mol.add_atom(at1,  0,  0,  0, unit=unit)
        mol.add_atom(at2,  1,  1,  1, unit=unit)
        mol.add_atom(at2,  1, -1, -1, unit=unit)
        mol.add_atom(at2, -1, -1,  1, unit=unit)
        mol.add_atom(at2, -1,  1,  1, unit=unit)

        return mol
