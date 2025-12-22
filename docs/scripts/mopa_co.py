from pyqint import MoleculeBuilder, HF, MOPA

mol = MoleculeBuilder.from_name('CO')
res = HF().rhf(mol, 'sto3g')
mopa = MOPA(res)

mohp = mopa.mohp(0, 1)
moop = mopa.moop(0, 1)

print('MOHP and MOOP values of canonical Hartree–Fock orbitals')
for i, (e, h, o) in enumerate(zip(res['orbe'], mohp, moop)):
    print('%3i %12.4f %12.4f %12.4f' % (i+1, e, h, o))