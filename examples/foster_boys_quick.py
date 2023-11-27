from pyqint import Molecule, HF, COHP, FosterBoys
import numpy as np

d = 1.145414
mol = Molecule()
mol.add_atom('C', 0.0, 0.0, -d/2, unit='angstrom')
mol.add_atom('O', 0.0, 0.0,  d/2, unit='angstrom')

res = HF().rhf(mol, 'sto3g')
cohp = COHP(res).run(res['orbc'], 0, 1)

resfb = FosterBoys(res).run(nr_runners=10)
cohp_fb = COHP(res).run(resfb['orbc'], 0, 1)

print('COHP values of canonical Hartree-Fock orbitals')
for i,(e,chi) in enumerate(zip(res['orbe'], cohp)):
    print('%3i %12.4f %12.4f' % (i+1,e,chi))
print()

print('COHP values after Foster-Boys localization')
for i,(e,chi) in enumerate(zip(resfb['orbe'], cohp_fb)):
    print('%3i %12.4f %12.4f' % (i+1,e,chi))
print()

print('Sum of COHP coefficient canonical orbitals: ', np.sum(cohp[:7]))
print('Sum of COHP coefficient Foster-Boys orbitals: ', np.sum(cohp_fb[:7]))
print('Result FB: ', resfb['r2start'], resfb['r2final'])
