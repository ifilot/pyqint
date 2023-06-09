# PyQInt

[![build](https://github.com/ifilot/pyqint/actions/workflows/build.yml/badge.svg)](https://github.com/ifilot/pyqint/actions/workflows/build.yml)
[![docs](https://github.com/ifilot/pyqint/actions/workflows/docs.yml/badge.svg)](https://github.com/ifilot/pyqint/actions/workflows/docs.yml)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/version.svg)](https://anaconda.org/ifilot/pyqint)
[![PyPI](https://img.shields.io/pypi/v/pyqint?style=flat-square)](https://pypi.org/project/pyqint/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

__Table of Contents__

* [Purpose](#purpose)
* [Documentation](#documentation)
* [Examples](#documentation)
    - [Contour plots of the molecular orbitals of CO](#contour-plots-of-the-molecular-orbitals-of-co)

## Purpose

PyQInt is a Python package for calculating one- and two-electron integrals as
encountered in electronic structure calculations. Since integral evaluation can
be quite computationally intensive, they are programmed in C++ and connected to
Python using Cython. PyQInt offers additional features such as performed restricted **Hartree-Fock
calculations**, calculation of **Crystal Orbital Hamilton Population** coefficients
and constructing localized orbitals using the **Boys-Foster** method.

## Documentation

PyQInt comes with detailed documentation and examples, which can be found
at https://pyqint.imc-tue.nl.

## Examples

### Contour plots of the molecular orbitals of CO

Below, the molecular orbitals for CO using a `sto-3g` basis set are shown.

![Molecular orbitals of CO](img/co.jpg)
