from pyqint import MoleculeBuilder, HF, ContourPlotter

mol = MoleculeBuilder.from_name('CH4')
res = HF(mol, 'sto3g').rhf(verbose=True)

up = [0,0,1]
ContourPlotter.build_contourplot(res, 'ch4.png', (0,1,2,up), 3, 101, 3, 3)