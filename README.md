# PyQInt

[![status](https://jose.theoj.org/papers/2a73fa24200e8c1ec47fc6e37f818a54/status.svg)](https://jose.theoj.org/papers/2a73fa24200e8c1ec47fc6e37f818a54)
[![Conda pkg](https://github.com/ifilot/pyqint/actions/workflows/build_conda.yml/badge.svg)](https://github.com/ifilot/pyqint/actions/workflows/build_conda.yml)
[![PyPI pkg](https://github.com/ifilot/pyqint/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/ifilot/pyqint/actions/workflows/build_wheels.yml)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/version.svg)](https://anaconda.org/ifilot/pyqint)
[![PyPI](https://img.shields.io/pypi/v/pyqint?style=flat-square)](https://pypi.org/project/pyqint/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Purpose

PyQInt is a Python-based, teaching-oriented implementation of the Hartree-Fock
method, designed to make the inner workings of electronic structure theory
accessible and transparent. It provides a clear, readable interface to
fundamental components such as molecular integrals over Gaussian basis
functions, SCF procedures (with DIIS acceleration), orbital localization, and
geometry optimization.

What sets PyQInt apart is its educational design philosophy: all matrices,
intermediate results, and algorithmic steps are exposed—allowing students,
educators, and developers to inspect, understand, and experiment with every part
of the computation. Whether you are learning how Hartree-Fock works, developing
your own extensions, or teaching a course in computational chemistry, PyQInt
offers a hands-on, exploratory platform.

For students interested in the theoretical foundations and algorithmic
implementation of Hartree–Fock, we recommend the open-access textbook Elements
of [Electronic Structure Theory](https://ifilot.pages.tue.nl/elements-of-electronic-structure-theory/index.html).

> [!NOTE] 
> PyQInt connects to a C++ backend for core numerical routines, but it
> is not optimized for performance. It is best suited for learning, prototyping,
> and small molecule calculations - **not** production-scale quantum chemistry.

> [!TIP]  
> Interested in other **education** quantum chemical codes? Have a look at the
> packages below.
> * [PyDFT](https://github.com/ifilot/pydft) is a pure-Python density functional
>   theory code, built on top of PyQInt.
> * [HFCXX](https://github.com/ifilot/hfcxx) is a full C++ code for performing
>   Hartree-Fock calculations.
> * [DFTCXX](https://github.com/ifilot/dftcxx) is a full C++ code for performing
>   Density Functional Theory Calculations.

## Documentation

PyQInt comes with detailed documentation and examples, which can be found
at https://ifilot.github.io/pyqint/.

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
* Performing [restricted
  Hartree-Fock](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method)
  calculations using [DIIS](https://en.wikipedia.org/wiki/DIIS)
* Calculation of [Crystal Orbital Hamilton Population](http://www.cohp.de/)
  coefficients
* Construction of localized orbitals using the [Foster-Boys
  method](https://en.wikipedia.org/wiki/Localized_molecular_orbitals#Foster-Boys)
* Geometry optimization using Conjugate Gradient
* Visualization of molecular orbitals

All routines are (automatically) tested and verified against several open-source
as well as commercial programs that use cartesian Gaussian orbitals.
Nevertheless, if you spot any mistake, please kindly open an
[issue](https://github.com/ifilot/pyqint/issues) in this Github repository.

In the image below, the (canonical) molecular orbitals as found using a
restricted Hartree-Fock calculation for the CO molecule are shown.

![Molecular orbitals of CO](img/co.jpg)

## Community guidelines

* Contributions to PyQInt are always welcome and appreciated. Before doing so,
  please first read the [CONTRIBUTING](CONTRIBUTING.md) guide.
* For reporting issues or problems with the software, you are kindly invited to
  to open a [new issue with the bug
  label](https://github.com/ifilot/pyqint/issues/new?labels=bug).
* If you seek support in using PyQInt, please [open an issue with the
  question](https://github.com/ifilot/pyqint/issues/new?labels=question) label.
* If you wish to contact the developers, please send an e-mail to ivo@ivofilot.nl.

## License

Unless otherwise stated, all code in this repository is provided under the GNU
General Public License version 3.