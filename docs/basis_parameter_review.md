# Review of all packaged basis parameters

Scope: every exponent, contraction coefficient, angular-momentum assignment,
element label and nuclear charge in the seven packaged basis JSON files.
This is a parameter/provenance audit, not an exhaustive validation of all
integral algorithms, SCF states, or molecular geometries.

## Conclusions

The boron/STO-6G defect identified in the preceding investigation is the only
entry in the released 1.4.3 wheel that fails to agree with **either** BSE v0
or BSE v1 within rounding. After that correction, all 276 element/basis
entries match at least one of these references.

That does **not** certify the legacy data as correct. Several inherited
definitions differ substantially from BSE's curated v1 data, and some cause
large energy differences. An independent code using the same legacy data can
reproduce the same problem. In particular, default PySCF tables are not an
independent modern parameter reference for these cases.

Recommendation: migrate the legacy STO and Pople tables to explicitly pinned
BSE v1 data in a separately reviewed change, documenting intentional variants
and numerical changes. Preserve the legacy tables under explicit versioned
names if compatibility is required. Updating every table silently under its
existing name would change reproducibility of previous calculations.

## Inventory and results

Audited 3,274 contracted shells and 17,079 primitive records, including
repeated primitives belonging to separate contractions. Reference library:
`basis_set_exchange==0.12`, comparing basis versions 0 and 1 separately.

| Packaged file | Elements | Shells | Primitive records | Agree with v0 | Agree with v1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| aug-cc-pVDZ | 34 | 444 | 3,193 | 30 | 34 |
| aug-cc-pVTZ | 34 | 622 | 4,095 | 30 | 34 |
| aug-cc-pVQZ | 34 | 834 | 5,292 | 30 | 34 |
| p321 (3-21G) | 55 | 545 | 1,248 | 55 | 54 |
| p631 (6-31G) | 30 | 228 | 824 | 30 | 28 |
| sto3g | 53 | 393 | 1,179 | 53 | 42 |
| sto6g, with B correction | 36 | 208 | 1,248 | 36 | 19 |

All three augmented correlation-consistent tables match BSE v1 **exactly**,
including all numerical values, after splitting general contractions and
removing zero padding. Their Li, Be, Na and Mg differences against v0 reflect
the choice of a newer reference, rather than corruption of these local files.

No invalid atomic numbers, duplicate JSON keys, nonfinite numbers,
nonpositive exponents, empty contractions or nonpositive contraction norms
were found. Aluminum/aluminium and cesium/caesium are accepted spelling aliases.

## Findings requiring scientific attention

### High: inherited STO-6G data, especially magnesium

Mg's core exponents are identical to Na's. Its other shells also differ from
the curated definition. For example, the first core exponent is 2600.756771
instead of 3103.386324. This mistake is present in the legacy reference and
in PySCF's default table, so agreement between those tables is insufficient.

For Sc through Kr, the third S/P shell retains the older coefficients rather
than the reoptimized 3sp coefficients in v1. The first coefficient of that S
shell is -0.007943 instead of -0.009737395526. BSE's
[STO family notes](https://www.basissetexchange.org/family_notes/sto/) explain
the reoptimization and scaling conventions. Co additionally has documented
scaling differences. These distinctions should be recorded, not collectively
described as local transcription errors.

The full list differing from v1 is: Mg, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu,
Zn, Ga, Ge, As, Se, Br, Kr. Boron agrees with v1 after the preceding fix.

### High: inherited STO-3G data for Y through Cd

Y has several substantial exponent differences; the first exponent in its
inner D shell is 5.660043 instead of 28.96238417. Zr's outer D and S/P shells
have exchanged exponent sets relative to v1. For Y through Cd, a coefficient
in the fourth S contraction is 1.155787 instead of 1.115578745. These are
not six-decimal rounding effects. Co also differs through documented scaling
choices.

The full list differing from v1 is: Co, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd.

### Medium: copper/3-21G

Cu's outermost S/P exponent is 0.000751 in the legacy table versus
0.04054498690 in v1, over a factor of 50 apart. Its preceding S/P contraction
also differs. Both definitions exist in reference databases; this is an
inherited discrepancy, not a failure of JSON loading. Cu is the only 3-21G
element with a material normalized-function difference from v1.

### Version choice: potassium and calcium/6-31G

BSE v1 contains two additional valence D shells for K and Ca; the packaged
tables and PySCF's defaults use the older definitions without them. BSE's
[Pople notes](https://www.basissetexchange.org/family_notes/pople/) attribute
this difference to the 2001 revision. This requires an explicit convention
choice, rather than an automatic missing-shell correction.

### Precision and equivalent representations

The STO/Pople files generally use six decimal places. For very diffuse
exponents, fixed absolute rounding discards meaningful relative precision.
For example, K/3-21G uses 0.016063 instead of 0.01606254630. Preserve the
source precision in future imports.

Some apparent coefficient discrepancies are harmless: BSE has normalized
some contractions, and others differ by an overall sign. PyQInt normalizes
contractions internally. The audit matches same-angular-momentum normalized
functions, allowing a whole-contraction sign, rather than treating every
raw coefficient difference as an error. It also reports raw maximum
parameter differences and exact numerical equality separately.

## Independent numerical checks

The following PySCF 2.14 calculations use explicit custom tables from each
source, Cartesian functions, bohr units, energy tolerance 1e-10, and the same
spin in both comparisons. These calculations converged. Open-shell energies
are converged UHF solutions, not proofs of the globally lowest atomic state.

| Atom/basis | Method | Packaged parameters | BSE v1 parameters | Difference, packaged minus v1 |
| --- | --- | ---: | ---: | ---: |
| Mg/STO-6G | RHF | -197.025789410 | -198.660064861 | +1.634275451 |
| Y/STO-3G | UHF doublet | -3252.994781730 | -3299.027651261 | +46.032869531 |
| Cu/3-21G | UHF doublet | -1630.655242731 | -1630.715764309 | +0.060521577 |

PyQInt reproduces the packaged-data Mg and Y energies to within 3e-9 hartree,
isolating these large differences to the selected data. Cu's default PyQInt
run did not agree with the converged PySCF result; it is excluded as a
successful engine comparison. Exploratory Zn RHF runs did not converge in
PySCF either and are excluded from this table. Parameters alone do not
establish SCF convergence, and difficult atomic calculations need explicit
iteration/state checks.

## Importer defect corrected during this review

`scripts/bse_to_pyqint.py` previously assigned every coefficient column of a
combined SP or SPD shell to its first angular momentum. An SP shell became
two S contractions, and an SPD shell became three S contractions.

Corrected the angular-momentum mapping, added checks for mismatched column
lengths, and retained all nonzero coefficients instead of dropping values
below 1e-12. Added four regression tests for combined shells, general
contractions, tiny nonzero coefficients and malformed data. The supplied
augmented input files use single-angular-momentum shells; regenerating all
three produces exactly their existing output after this fix.

This converter is still scoped to the supported element/name map (through
Kr) and all-electron basis data. It is not a general ECP importer.

The review adds no further numerical table replacements beyond the boron
correction from the preceding investigation.

## Reproduce the full parameter audit

Install the optional audit dependency `basis_set_exchange==0.12`, then run:

```sh
python scripts/audit_basis_parameters.py --output basis-audit.json
```

The script uses the installed BSE dataset without network access. It records
every element's structural checks, v0/v1 comparisons, exact equality,
normalized-function differences, and raw maximum parameter differences. The
default normalized-function screening tolerance is 1e-4 in L2 norm; this is
a rounding/provenance screen, not a bound on molecular energy errors.

To audit an installed 1.4.3 wheel, pass its `pyqint/basissets` directory via
`--basis-dir`. That audit exits nonzero and identifies B/STO-6G as unmatched;
the corrected checkout has no unmatched entries. Legacy matches remain
visible alongside v1 differences and must not be read as scientific approval.

Validation: six new tests pass. The broader targeted run passed 36 tests and
failed the same pre-existing N2-anion/6-31G test documented in the boron
investigation. JSON regeneration checks and `git diff --check` pass.
