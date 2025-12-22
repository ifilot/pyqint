from pyqint import MoleculeBuilder

mol = MoleculeBuilder.from_name('ch4')
mol.name = 'CH4'

print(mol)