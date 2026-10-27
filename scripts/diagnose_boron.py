"""Reproduce the boron basis report with the installed PyQInt package.

Example (the full BH3/aug-cc-pVTZ calculation can take several minutes)::

    OMP_NUM_THREADS=2 python scripts/diagnose_boron.py --pyscf --output results.json

PySCF is optional. Atomic comparisons use UHF doublets; BH3 uses RHF.
Cartesian basis functions are used in PySCF to match PyQInt.
"""

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
import pyqint
from pyqint import HF, Molecule, PyQInt


def geometry(system):
    if system in ('H', 'B'):
        return [(system, (0.0, 0.0, 0.0))]
    return [('B', (0.0, 0.0, 0.0))] + [
        ('H', (2.19 * np.cos(t), 2.19 * np.sin(t), 0.0))
        for t in (0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0)
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--systems', nargs='+', choices=['H', 'B', 'BH3'],
                        default=['H', 'B', 'BH3'])
    parser.add_argument('--basis', nargs='+',
                        default=['sto3g', 'sto6g', 'aug-cc-pVDZ', 'aug-cc-pVTZ'])
    parser.add_argument('--pyscf', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'pyqint_version': pyqint.__version__,
              'pyqint_path': pyqint.__file__, 'python': platform.python_version(),
              'numpy': np.__version__, 'results': []}
    if args.pyscf:
        import pyscf
        from pyscf import gto, scf
        report['pyscf_version'] = pyscf.__version__
    print(json.dumps({k: v for k, v in report.items() if k != 'results'}), flush=True)

    for system in args.systems:
        atoms = geometry(system)
        for basis in args.basis:
            mol = Molecule(system)
            for atom, pos in atoms:
                mol.add_atom(atom, *pos, unit='bohr')
            started = time.monotonic()
            solver = HF(mol, basis)
            if system == 'BH3':
                result = solver.rhf(tolerance=1e-9)
            else:
                result = solver.uhf(multiplicity=2, tolerance=1e-9)
            # Evaluate overlap separately to also support older result dictionaries.
            cgfs, _ = mol.build_basis(basis)
            q = PyQInt()
            overlap = np.array([[q.overlap(a, b) for b in cgfs] for a in cgfs])
            eigenvalues = np.linalg.eigvalsh(overlap)
            row = {'system': system, 'basis': basis,
                   'energy': float(result['energy']), 'n_basis': len(cgfs),
                   'overlap_min': float(eigenvalues[0]),
                   'overlap_max': float(eigenvalues[-1]),
                   'elapsed_seconds': time.monotonic() - started}
            if 'energies' in result:
                row['iterations'] = len(result['energies'])
                row['energy_history'] = [float(e) for e in result['energies']]
            if args.pyscf:
                ref = gto.M(atom=atoms, unit='Bohr', cart=True, verbose=0,
                            spin=0 if system == 'BH3' else 1,
                            basis={'sto3g': 'sto-3g', 'sto6g': 'sto-6g',
                                   'p631': '6-31g'}.get(basis, basis))
                mf = scf.RHF(ref) if system == 'BH3' else scf.UHF(ref)
                mf.conv_tol = 1e-11
                row['pyscf_energy'] = float(mf.kernel())
                row['pyscf_converged'] = bool(mf.converged)
                row['energy_difference'] = row['energy'] - row['pyscf_energy']
            report['results'].append(row)
            print(json.dumps(row), flush=True)
            if args.output:
                args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
