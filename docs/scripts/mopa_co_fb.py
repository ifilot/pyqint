from pyqint import MoleculeBuilder, HF, FosterBoys, PopulationAnalysis
import numpy as np

mol = MoleculeBuilder.from_name('CO')
res = HF(mol, 'sto3g').rhf()
pa = PopulationAnalysis(res)
mohp = pa.mohp(0, 1)
moop = pa.moop(0, 1)

res_fb = FosterBoys(res).run()
pa_fb = PopulationAnalysis(res_fb)
mohp_fb = pa_fb.mohp(0, 1)
moop_fb = pa_fb.moop(0, 1)

# the sums over the occupied orbitals are invariant under orbital localization
nocc = res['nelec'] // 2

print('Sum of MOHP (canonical orbitals):',
        np.sum(mohp[:nocc]))
print('Sum of MOHP (Foster-Boys orbitals):',
        np.sum(mohp_fb[:nocc]))

print('Sum of MOOP (canonical orbitals):',
        np.sum(moop[:nocc]))
print('Sum of MOOP (Foster-Boys orbitals):',
        np.sum(moop_fb[:nocc]))
