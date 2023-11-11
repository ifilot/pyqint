# PyQInt

[![build](https://github.com/ifilot/pyqint/actions/workflows/build.yml/badge.svg)](https://github.com/ifilot/pyqint/actions/workflows/build.yml)
[![docs](https://github.com/ifilot/pyqint/actions/workflows/docs.yml/badge.svg)](https://github.com/ifilot/pyqint/actions/workflows/docs.yml)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/version.svg)](https://anaconda.org/ifilot/pyqint)
[![PyPI](https://img.shields.io/pypi/v/pyqint?style=flat-square)](https://pypi.org/project/pyqint/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Purpose

PyQInt is a Python package for calculating one- and two-electron integrals as
encountered in electronic structure calculations. Since integral evaluation can
be quite computationally intensive, the evaluation is programmed in C++ and
connected to Python using Cython.

PyQInt mainly serves as an educational package to teach students how to perform
(simple) electronic structure calculations wherein the most difficult task,
i.e. the integral evaluation, is already encapsulated in a handy set of
routines. With PyQInt, the student can for example build their own Hartree-Fock
routine. Some common electronic structure routine, most notably the
Hartree-Fock algorithm, is also readily available.

> **Note**
> Although PyQInt connects to a C++ backend, it is certainly not optimized for
> speed and might be (too) slow for anything outside of the calculation of the
> electronic structure of simple molecules.

## Documentation

PyQInt comes with detailed documentation and examples, which can be found
at https://pyqint.imc-tue.nl.

## Features

The following molecular integrals are supported by PyQInt

- [x] Overlap integral
- [x] Kinetic integral
- [x] Dipole integral
- [x] Nuclear integral
- [x] Two-electron repulsion integral

as well as the following geometric derivatives

- [x] Overlap integral
- [x] Kinetic integral
- [x] Nuclear integral
- [x] Two-electron repulsion integral

PyQInt offers additional features such as
* Performing [restricted Hartree-Fock](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method)
  calculations using [DIIS](https://en.wikipedia.org/wiki/DIIS)
* Calculation of [Crystal Orbital Hamilton Population](http://www.cohp.de/) coefficients
* Construction of localized orbitals using the [Boys-Foster method](https://en.wikipedia.org/wiki/Localized_molecular_orbitals#Foster-Boys)
* Visualization of molecular orbitals

All routines are (automatically) tested and verified against several open-source
as well as commercial programs that use cartesian Gaussian orbitals. Nevertheless,
if you spot any mistake, please kindly open an [issue](https://github.com/ifilot/pyqint/issues)
in this Github repository.

In the image below, the (canonical) molecular orbitals as found using a restricted
Hartree-Fock calculation for the CO molecule are shown.

![Molecular orbitals of CO](img/co.jpg)
