from pyqint import MoleculeBuilder, HF, ContourPlotter

mol = MoleculeBuilder.from_name('borane')
res = HF(mol, 'sto3g').rhf(verbose=True)

planes = ['xy' for i in range(8)]
planes[4] = 'xz'
ContourPlotter.build_contourplot(res, 'borane.png', planes, 3, 101, 2, 4)