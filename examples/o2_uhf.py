from pyqint import MoleculeBuilder, HF

mol = MoleculeBuilder.from_name('O2')
res = HF(mol, 'sto3g').uhf(verbose=True, multiplicity=3)

print(res['orbe_alpha'])
print(res['orbe_beta'])