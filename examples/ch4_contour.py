from pyqint import MoleculeBuilder, HF, ContourPlotter

mol = MoleculeBuilder().from_name('CH4')
res = HF().rhf(mol, 'sto3g', verbose=True)

up = [0,0,1]
ContourPlotter().build_contourplot(res, 'ch4.png', [0,1,2,up], 3, 101, 3, 3)