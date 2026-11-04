from pyqint import Molecule, HF, PopulationAnalysis, FosterBoys
import numpy as np

#
# Molecular Orbital Hamilton Population (MOHP) analysis of the C-O bond in
# CO, using the canonical and the Foster-Boys localized molecular orbitals
#

d = 1.145414
mol = Molecule()
mol.add_atom('C', 0.0, 0.0, -d/2, unit='angstrom')
mol.add_atom('O', 0.0, 0.0,  d/2, unit='angstrom')

res = HF(mol, 'sto3g').rhf()
mohp = PopulationAnalysis(res).mohp(0, 1)

resfb = FosterBoys(res).run(nr_runners=10)
mohp_fb = PopulationAnalysis(resfb).mohp(0, 1)

print('MOHP values of canonical Hartree-Fock orbitals')
for i,(e,chi) in enumerate(zip(res['orbe'], mohp)):
    print('%3i %12.4f %12.4f' % (i+1,e,chi))
print()

print('MOHP values after Foster-Boys localization')
for i,(e,chi) in enumerate(zip(resfb['orbe'], mohp_fb)):
    print('%3i %12.4f %12.4f' % (i+1,e,chi))
print()

# the sum over the occupied orbitals is invariant under orbital localization
nocc = res['nelec'] // 2
print('Sum of MOHP coefficients canonical orbitals:   ', np.sum(mohp[:nocc]))
print('Sum of MOHP coefficients Foster-Boys orbitals: ', np.sum(mohp_fb[:nocc]))
print('Result FB: ', resfb['r2start'], resfb['r2final'])
