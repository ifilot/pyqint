from pyqint import Molecule, HF, FosterBoys, PopulationAnalysis, ContourPlotter
import numpy as np

# planar methyl radical (doublet)
R = 2.039
sqrt3 = np.sqrt(3.0)
mol = Molecule()
mol.add_atom('C', 0.0, 0.0, 0.0)
mol.add_atom('H',  R, 0.0, 0.0)
mol.add_atom('H', -0.5 * R,  0.5 * sqrt3 * R, 0.0)
mol.add_atom('H', -0.5 * R, -0.5 * sqrt3 * R, 0.0)

res = HF(mol, 'sto3g').uhf(multiplicity=2)

# atomic charges and spin populations
pa = PopulationAnalysis(res)
print('Atom   Charge (M)   Spin (M)   Charge (L)   Spin (L)')
for n, atom in enumerate(mol):
    print('%2s  %12.6f %10.6f %12.6f %10.6f' % (atom[0],
          pa.mulliken(n), pa.mulliken_spin(n),
          pa.lowdin(n), pa.lowdin_spin(n)))

# spin-resolved MOHP for the first C-H bond
for spin in ('alpha', 'beta'):
    nocc = res['n' + spin]
    mohp = pa.mohp(0, 1, spin=spin)
    print('Sum of MOHP (%5s):' % spin, np.sum(mohp[:nocc]))

# localize the alpha and beta orbitals separately
res_fb = FosterBoys(res, seed=0).run(nr_runners=3)
pa_fb = PopulationAnalysis(res_fb)
for spin in ('alpha', 'beta'):
    nocc = res['n' + spin]
    mohp_fb = pa_fb.mohp(0, 1, spin=spin)
    print('Sum of MOHP (%5s, Foster-Boys):' % spin, np.sum(mohp_fb[:nocc]))

# contour plots of the localized alpha and beta orbitals
for spin in ('alpha', 'beta'):
    ContourPlotter.build_contourplot(
        res_fb,
        'ch3_fb_contour_%s.png' % spin,
        plane='xy',
        sz=4.0,
        npts=101,
        nrows=1,
        ncols=5,
        spin=spin,
    )
