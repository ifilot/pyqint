from pyqint import MoleculeBuilder, HF, FosterBoys, MOPA
import numpy as np

mol = MoleculeBuilder.from_name('CO')
res = HF().rhf(mol, 'sto3g')
mopa = MOPA(res)
mohp = mopa.mohp(0, 1)
moop = mopa.moop(0, 1)

res_fb = FosterBoys(res).run()
mopa_fb = MOPA(res_fb)
mohp_fb = mopa_fb.mohp(0, 1)
moop_fb = mopa_fb.moop(0, 1)

print('Sum of MOHP (canonical orbitals):',
        np.sum(mohp[:7]))
print('Sum of MOHP (Foster-Boys orbitals):',
        np.sum(mohp_fb[:7]))

print('Sum of MOOP (canonical orbitals):',
        np.sum(moop[:7]))
print('Sum of MOOP (Foster-Boys orbitals):',
        np.sum(moop_fb[:7]))