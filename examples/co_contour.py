from pyqint import MoleculeBuilder, HF, ContourPlotter

mol = MoleculeBuilder.from_name('CO')
res = HF(mol, 'sto3g').rhf(verbose=True)

ContourPlotter.build_contourplot(
    res,
    'co.png',
    plane='yz',
    sz=3.0,
    npts=101,
    nrows=2,
    ncols=5,
    levels=15
)