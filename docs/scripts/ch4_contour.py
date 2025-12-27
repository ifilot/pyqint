from pyqint import MoleculeBuilder, HF, ContourPlotter

mol = MoleculeBuilder.from_name('CH4')
res = HF(mol, 'sto3g').rhf(verbose=True)

up = [0, 0, 1]
ContourPlotter.build_contourplot(
    res,
    'ch4_contour.png',
    plane=[0, 1, 2, up],
    sz=3.0,
    npts=101,
    nrows=3,
    ncols=3
)