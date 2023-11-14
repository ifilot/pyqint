from pyqint import Molecule, PyQInt, FosterBoys, GeometryOptimization, HF
from pyqint import BlenderRender
import pyqint
import numpy as np
import matplotlib.pyplot as plt
import os

#
# Plot the isosurfaces for the CO molecule
#

print(pyqint.__version__)
outpath = os.path.dirname(__file__)

def main():   
    mol = Molecule('CO')
    mol.add_atom('C', 0, 0, -1.08232106)
    mol.add_atom('O', 0, 0, 1.08232106)
    res = HF().rhf(mol, 'sto3g')
    resfb = FosterBoys(res).run()

    br = BlenderRender()
    br.render_molecular_orbitals(res['mol'], res['cgfs'], resfb['orbc'], outpath)

if __name__ == '__main__':
    main()
