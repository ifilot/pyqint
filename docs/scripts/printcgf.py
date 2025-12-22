from pyqint import PyQInt, Molecule
import numpy as np

# construct integrator object
integrator = PyQInt()

# build hydrogen molecule
mol = Molecule('H2')
mol.add_atom('H', 0.0, 0.0, 0.0)
mol.add_atom('H', 0.0, 0.0, 1.4)
cgfs, nuclei = mol.build_basis('sto3g')

for cgf in cgfs:
    print(cgfs)