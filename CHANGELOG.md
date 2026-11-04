# Changelog

All notable changes to PyQInt are documented in this file. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and PyQInt adheres
to [Semantic Versioning](https://semver.org/).

Entries for versions before 1.5.0 were reconstructed from the release notes and
pull requests on GitHub. Versions before 1.0.0 are not covered; see the
[git history](https://github.com/ifilot/pyqint/commits/master) for those.

## [1.5.1] - 2026-11-04

### Added

- `PyQInt.plot_basis_functions(grid, cgfs)` evaluates all basis functions on a
  grid in a single call and returns an array of shape
  `(n_basis_functions, n_points)`.
- `PyQInt.plot_basis_function_gradients(grid, cgfs)` evaluates the gradients
  of all basis functions on a grid in a single call and returns an array of
  shape `(n_basis_functions, n_points, 3)`.

### Performance

- The two new routines evaluate all basis functions within a single OpenMP
  parallel region. Evaluating basis functions on a grid one at a time via
  `plot_basis_function` or `plot_wavefunction` starts a parallel region per
  basis function, whose overhead dominates for small grids and many threads.

## [1.5.0] - 2026-11-04

### Added

- Unrestricted Hartree-Fock (UHF) results are supported by `PopulationAnalysis`
  (Mulliken and Löwdin spin populations via `mulliken_spin` and `lowdin_spin`;
  per-spin MOOP, MOHP and MOBI via `spin='alpha'` or `spin='beta'`),
  `FosterBoys` (each spin channel is localized independently) and
  `ContourPlotter` (`spin` argument).
- UHF results report the spin contamination `<S^2>` (`s2`, together with the
  exact value `s2_exact`) and the electron counts `nalpha` and `nbeta`.
- SCF results report whether the calculation `converged` and the number of
  iterations `niter`; a `RuntimeWarning` is issued when the SCF does not
  converge. A new `tolerance_grad` argument sets the threshold for the orbital
  gradient (the commutator FPS - SPF).
- `scripts/build_basis_sets.py` regenerates all packaged basis sets from the
  Basis Set Exchange and verifies them (`--check`). `scripts/bse_to_pyqint.py`
  can read basis sets directly from the `basis_set_exchange` package
  (`--basis`, `--version`, `--elements`) and supports elements up to Rn.

### Changed

- **Package layout.** The modules are grouped into subpackages: `basis`,
  `structure`, `methods`, `analysis` and `visualization`. Imports from the top
  level (`from pyqint import HF, Molecule, ...`) are unchanged, but imports of
  modules have moved, e.g. `pyqint.hf` is now `pyqint.methods.hf` and
  `pyqint.molecule` is now `pyqint.structure.molecule`.
- **Basis sets.** All packaged basis sets are now verbatim conversions of
  version 1 of the corresponding Basis Set Exchange basis sets. This corrects
  erroneous parameters (see *Fixed*), restores the full precision of exponents
  and coefficients (previously rounded to six decimals) and adds the 2001
  revision d-shells of K and Ca in 6-31G. Calculated energies therefore change
  slightly for all elements. For Sc and heavier elements in STO-3G, STO-6G and
  3-21G, the basis functions now follow the Basis Set Exchange shell ordering
  (d-shells after the s- and p-shells), which changes the order of the atomic
  orbitals.
- **Boys function.** The Boys function is evaluated with a tabulated Taylor
  expansion to near machine precision for all orders and arguments. Values of
  two-electron and nuclear attraction integrals change by up to about 4e-6.
- **Basis function normalization.** Contracted Gaussian functions are always
  normalized as a whole, and amplitudes and gradients (`get_amp`, `get_grad`,
  `plot_wavefunction`, `plot_gradient`) now include this normalization,
  consistent with all integrals. Packaged basis sets are effectively
  unaffected; plots based on custom, unnormalized contractions change (and are
  now correct).
- **SCF convergence.** Convergence requires both the energy change and the
  orbital gradient to be below their thresholds, and the DIIS error vectors
  are computed from the density matrix that produced each Fock matrix. Canonical
  orthogonalization removes near-linear dependencies in the basis set
  (overlap eigenvalues below 1e-7), such that the number of molecular orbitals
  can be smaller than the number of basis functions.
- RHF raises a `ValueError` for an odd number of electrons; UHF raises a
  `ValueError` for a multiplicity that is incompatible with the number of
  electrons or the basis set.
- Package data (basis sets, molecules) is accessed via `importlib.resources`.

### Fixed

- Nonphysical, far too negative energies (e.g. -7501 Ht instead of -25.02 Ht
  for a compressed BH3 molecule in aug-cc-pVDZ) caused by the inaccurate
  interpolated Boys function, which rendered the two-electron integrals
  indefinite for diffuse basis sets.
- Incorrect basis set parameters: all exponents of B in STO-6G and the
  exponents of Mg in STO-6G (which were copied from Na). Outdated or deviating
  parameters were replaced for Sc-Kr in STO-6G, Co and Y-Cd in STO-3G and Cu
  in 3-21G.
- `scripts/bse_to_pyqint.py` assigned all columns of combined SP and SPD
  shells to their first angular momentum.
- Analytical forces failed when the basis set was supplied as a list of CGFs.
- `HF.rhf(itermax=0)` raised an `UnboundLocalError`.
- A failed DIIS extrapolation consumed an SCF iteration.
- Unpickled CGFs lost the positions of GTOs that are not located at the CGF
  center (e.g. symmetry-adapted basis functions).
- `MoleculeBuilder.from_file` raised a `NameError` for missing files.
- Compilation on macOS (missing `<cstddef>` include).
- Examples and documentation scripts that used the removed `COHP` and `MOPA`
  classes, the `purepytessel` package or an outdated `HF` API.

### Removed

- The undocumented `pyqint.isosurface` module, including its `.abo` writer.
  Use [PyTessel](https://github.com/ifilot/pytessel) for isosurfaces and
  [ABO Builder](https://github.com/ifilot/abobuilder) for `.abo` files.

### Documentation

- The integrals over individual GTOs (`overlap_gto`, `kinetic_gto`, ...) are
  documented to be evaluated over unnormalized primitives.
- New sections on UHF post-processing, basis function normalization and the
  provenance of the packaged basis sets.

## [1.4.3] - 2026-02-06

### Fixed

- `ContourPlotter` did not work for a single row of plots
  ([#76](https://github.com/ifilot/pyqint/issues/76)).
- Renormalization of CGFs did not work for symmetry-adapted CGFs
  ([#78](https://github.com/ifilot/pyqint/issues/78)).

## [1.4.2] - 2026-02-01

### Added

- Color schemes for the `BlenderRender` class (experimental).
- Compilation details via `integrator.get_compile_info()`.

### Fixed

- OpenMP parallelization did not work under Windows.
- Normalization of contracted Gaussian functions
  ([#72](https://github.com/ifilot/pyqint/issues/72)).

### Documentation

- Section on using PyQInt with [ABO Builder](https://github.com/ifilot/abobuilder)
  and [Managlyph](https://github.com/ifilot/managlyph).
- Fixed math rendering, labels and gallery issues
  ([#73](https://github.com/ifilot/pyqint/pull/73)) and added ethane and bond
  index examples ([#74](https://github.com/ifilot/pyqint/pull/74)).

## [1.4.1] - 2026-01-06

### Added

- Compilation settings are exposed via the integrator API.

### Fixed

- OpenMP parallelization did not work under Windows.

## [1.4.0] - 2026-01-04

### Changed

- Two-electron integrals are evaluated with a cached Hellsing-based scheme,
  replacing the Huzinaga kernel, with an average speed-up of 3.2x.
- Geometric derivatives of the nuclear attraction integrals are evaluated
  analytically instead of by finite differences.

## [1.3.1] - 2026-01-02

### Added

- Option to highlight elements in matrices in `MatrixPlotter`.
- Benchmark for the evaluation of two-electron integrals.

## [1.3.0] - 2025-12-28

### Added

- Unrestricted Hartree-Fock method.
- New molecules and expanded population analysis methods.

### Changed

- All two-electron integral routines use a cached Boys function, giving a
  15-20% speed-up.

## [1.2.0] - 2025-12-22

### Added

- `PopulationAnalysis` class with Molecular Orbital Hamilton and Overlap
  Populations and Mulliken and Löwdin charges, replacing `COHP`.
- `MatrixPlotter` class for visualizing matrices.
- `ContourPlotter` class for contour plots of all molecular orbitals of a
  calculation.
- Examples on symmetry-adapted basis sets and on molecular orbital diagrams
  using PyMoDia.
- More molecules in `MoleculeBuilder`.

### Changed

- The build system changed from `setuptools` to `meson`.
- Post-processing classes share a unified API (initialization from a
  Hartree-Fock result) and several class methods became static.
- Changed output of the `GeometryOptimization` class.
- Electronic steps of the Hartree-Fock routine are 4-100x faster.

### Removed

- `COHP` class (superseded by `PopulationAnalysis`).

### Fixed

- Warnings on missing PyTessel and tqdm packages when not using Blender
  rendering.

## [1.1.1] - 2025-11-23

### Removed

- Conda packages; PyQInt is distributed via PyPI.

## [1.1.0] - 2025-09-25

### Changed

- The `gto` and `cgf` classes are renamed to `GTO` and `CGF`.

## [1.0.1] - 2025-09-03

### Documentation

- Example Hartree-Fock calculation now saves and shows the orbitals
  ([#52](https://github.com/ifilot/pyqint/issues/52)).
- Fixed links to PyTessel ([#53](https://github.com/ifilot/pyqint/issues/53))
  and revised the installation instructions
  ([#54](https://github.com/ifilot/pyqint/issues/54),
  [#55](https://github.com/ifilot/pyqint/issues/55)).
- Expanded README and added contribution guidelines.

## [1.0.0] - 2025-04-24

First stable release.

[1.5.1]: https://github.com/ifilot/pyqint/compare/v1.5.0...v1.5.1
[1.5.0]: https://github.com/ifilot/pyqint/compare/v1.4.3...v1.5.0
[1.4.3]: https://github.com/ifilot/pyqint/compare/v1.4.2...v1.4.3
[1.4.2]: https://github.com/ifilot/pyqint/compare/v1.4.1...v1.4.2
[1.4.1]: https://github.com/ifilot/pyqint/compare/v1.4.0...v1.4.1
[1.4.0]: https://github.com/ifilot/pyqint/compare/v1.3.1...v1.4.0
[1.3.1]: https://github.com/ifilot/pyqint/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/ifilot/pyqint/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/ifilot/pyqint/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/ifilot/pyqint/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/ifilot/pyqint/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/ifilot/pyqint/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/ifilot/pyqint/releases/tag/v1.0.0
