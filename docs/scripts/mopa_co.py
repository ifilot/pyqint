from pyqint import MoleculeBuilder, HF, PopulationAnalysis

mol = MoleculeBuilder.from_name('CO')
res = HF(mol, 'sto3g').rhf()
pa = PopulationAnalysis(res)

mohp = pa.mohp(0, 1)
moop = pa.moop(0, 1)

print('MOHP and MOOP values of canonical Hartree–Fock orbitals')
for i, (e, h, o) in enumerate(zip(res['orbe'], mohp, moop)):
    print('%3i %12.4f %12.4f %12.4f' % (i+1, e, h, o))
