from pyqint import GeometryOptimization, Molecule

mol = Molecule()
dist = 1.0
mol.add_atom('C', 0.1, 0.0, 0.1, unit='angstrom')
mol.add_atom('H', dist, dist, dist, unit='angstrom')
mol.add_atom('H', -dist, -dist, dist, unit='angstrom')
mol.add_atom('H', -dist, dist, -dist, unit='angstrom')
mol.add_atom('H', dist, -dist, -dist, unit='angstrom')

res = GeometryOptimization(mol, 'sto3g', verbose=True).run()