# Boron basis investigation

The reported BH3 energy of approximately -2434 hartree was **not reproduced**.
A separate error in the boron STO-6G basis definition was reproduced and fixed.

## Environment and geometry

Tested the pip-installed PyQInt 1.4.3 Linux CPython 3.12 wheel and this repository
(1.5.0), using Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, PySCF 2.14.0, and
Basis Set Exchange 0.12. Calculations used two OpenMP threads and one OpenBLAS
thread. These results do not establish behavior on other platforms or wheels.

BH3 has B at the origin and three H atoms at radius 2.19 bohr, separated by
120 degrees in the xy-plane. Molecular calculations used RHF with default
settings. Atomic calculations used UHF doublets, not RHF with an odd electron
count. PySCF used Cartesian functions to match PyQInt. Basis filenames are
case-sensitive on Linux: use `aug-cc-pVTZ`, not `aug-cc-pvtz`, and `p631` for
the packaged 6-31G basis.

## Released-wheel results

Energies are in hartree. PySCF references use its standard basis data.

| Basis | H UHF, PyQInt | B UHF, PyQInt | B UHF, PySCF | BH3 RHF, PyQInt |
| --- | ---: | ---: | ---: | ---: |
| STO-3G | -0.466581918 | -24.148987970 | -24.148988599 | -26.070690348 |
| STO-6G (before fix) | -0.471039172 | -14.505516492 | -24.394294559 | -25.138679703 |
| aug-cc-pVDZ | -0.499334315 | -24.530725753 | -24.530725753 | -26.389197409 |
| aug-cc-pVTZ | -0.499833840 | -24.532257535 | -24.532257535 | -26.399179272 |

PySCF's BH3/aug-cc-pVTZ Cartesian RHF reference is -26.399190296 hartree.
The current 1.5.0 checkout returned the same -26.399179272 hartree as the
released 1.4.3 wheel (13 SCF iterations, approximately 132 seconds).
The remaining difference is about 0.000011 hartree; it is not evidence of the
reported catastrophic energy. Sampled molecular electron-repulsion integrals
showed differences up to about 4.4e-7 from PySCF. No change to the integral or
SCF code was made in this investigation.

## Isolation of the faults

The H and B aug-cc-pVDZ and aug-cc-pVTZ JSON entries match Basis Set Exchange
exactly, after expanding its general contractions into individual shells.
For each isolated atom, overlap, kinetic, nuclear-attraction and full
electron-repulsion tensors were compared with PySCF using PyQInt's own basis
data. AO ordering and Cartesian normalization were matched explicitly. Maximum
absolute integral discrepancies were below 6.2e-14. PySCF also reproduced the
bad B/STO-6G energy when supplied with the erroneous PyQInt data, isolating
that problem to the basis data rather than its integral implementation.

The smallest overlap eigenvalues were positive:

| System/basis | Minimum overlap eigenvalue |
| --- | ---: |
| H/aug-cc-pVTZ | 0.002141747 |
| B/aug-cc-pVTZ | 0.000821647 |
| BH3/aug-cc-pVTZ | 0.000020887 |

Positive overlap eigenvalues alone do not establish that a basis definition
is correct: even the erroneous B/STO-6G overlap matrix is positive definite.
Also, STO-3G and STO-6G are not nested variational spaces, so energy ordering
alone is not a general proof of an error. Here, the explicit data comparison
and independent calculation establish it.

## Correction and validation

All 18 exponents in B's three STO-6G contractions were incorrect. For example,
the core exponent 1355.584234 should be 506.0118369; the valence exponent
5.364202 should be 0.2384089591. Replaced B's exponents using
[Basis Set Exchange STO-6G data](https://www.basissetexchange.org/), retaining
the existing rounded contraction coefficients and all other element entries.

After correction:

| System/STO-6G | PyQInt | PySCF standard basis |
| --- | ---: | ---: |
| B UHF | -24.394294632 | -24.394294559 |
| BH3 RHF | -26.331115541 | -26.331129768 |

Added atomic and molecular regression tests in `tests/test_boron_basis.py`.
The targeted run passed 32 tests and failed one existing N2-anion/6-31G test
(`test_unrestricted_hartree_fock_n2_anion`). That same test also fails against
the unmodified pip-installed 1.4.3 wheel, with energy -108.746646987 instead
of the asserted -108.750462418. It does not involve boron or STO-6G.

## Reproduction

`scripts/diagnose_boron.py` works with the installed package, including 1.4.3.
It records the imported package path/version, energy histories and overlap
eigenvalue bounds, with optional Cartesian PySCF references and JSON output:

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python scripts/diagnose_boron.py \
    --pyscf --output boron-results.json
```

For a quick check, add `--systems B BH3 --basis sto6g`. A full
BH3/aug-cc-pVTZ run has 130 Cartesian basis functions and took approximately
136 seconds with the released wheel in this environment.

The original -2434 result still requires the original script, platform/wheel
details and iteration output to explain. The available handoff is insufficient
to identify its cause, and the confirmed STO-6G data defect does not explain
an aug-cc-pVTZ result.
