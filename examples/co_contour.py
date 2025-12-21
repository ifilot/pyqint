from pyqint import MoleculeBuilder, HF, ContourPlotter

mol = MoleculeBuilder.from_name('CO')
res = HF(mol, 'sto3g').rhf(verbose=True)

ContourPlotter().build_contourplot(res, 'co.png', 'yz', 3, 101, 2, 5)