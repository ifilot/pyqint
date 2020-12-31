# PyQInt

[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/version.svg)](https://anaconda.org/ifilot/pyqint)
[![PyPI](https://img.shields.io/pypi/v/pyqint?style=flat-square)](https://pypi.org/project/pyqint/)

__Table of Contents__

* [Purpose](#purpose)
* [Installation](#installation)
    - [Anaconda](#anaconda)
    - [PyPi](#pypi)
* [Usage](#usage)
    - [Overlap integrals](#overlap-integrals)
    - [Kinetic integrals](#kinetic-integrals)
    - [Nuclear attraction integrals](#nuclear-attraction-integrals)
    - [Two-electron integrals](#two-electron-integrals)

## Purpose

PyQInt is a Python package for calculating one- and two-electron integrals as encountered in electronic structure calculations. Since integral evaluation can be quite computationally intensive, they are programmed in C++ and connected to Python using Cython.

## Installation

### Anaconda

[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/version.svg)](https://anaconda.org/ifilot/pyqint)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/platforms.svg)](https://anaconda.org/ifilot/pyqint)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/downloads.svg)](https://anaconda.org/ifilot/pyqint)
[![Anaconda-Server Badge](https://anaconda.org/ifilot/pyqint/badges/installer/conda.svg)](https://conda.anaconda.org/ifilot)

Open Anaconda prompt and type

```
conda install -c ifilot pyqint
```

### PyPi

[![PyPI](https://img.shields.io/pypi/v/pyqint?color=green&style=flat-square)](https://pypi.org/project/pyqint/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/pypi?style=flat-square)](https://pypi.org/project/pyqint/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyqint?style=flat-square)


Open a terminal and type

```
pip install pyqint
```

## Usage

### Overlap integrals
```python
import pyqint
import numpy as np
from copy import deepcopy

# construct integrator object
integrator = PyQInt()

# build cgf for hydrogen separated by 1.4 a.u.
cgf1 = cgf([0.0, 0.0, 0.0])

cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

# create a copy of the CGF
cgf2 = deepcopy(cgf1)
cgf2.p[2] = 1.4

# construct empty matrix
S = np.zeros((2,2))
S[0,0] = integrator.overlap(cgf1, cgf1)
S[0,1] = S[1,0] = integrator.overlap(cgf1, cgf2)
S[1,1] = integrator.overlap(cgf2, cgf2)

# output result
print(S)
```

### Kinetic integrals
```python
import pyqint
import numpy as np
from copy import deepcopy

# construct integrator object
integrator = PyQInt()

# build cgf for hydrogen separated by 1.4 a.u.
cgf1 = cgf([0.0, 0.0, 0.0])

cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

# create a copy of the CGF
cgf2 = deepcopy(cgf1)
cgf2.p[2] = 1.4

# construct empty matrix
T = np.zeros((2,2))
T[0,0] = integrator.kinetic(cgf1, cgf1)
T[0,1] = T[1,0] = integrator.kinetic(cgf1, cgf2)
T[1,1] = integrator.kinetic(cgf2, cgf2)

# output result
print(T)
```

### Nuclear attraction integrals
```python
import pyqint
import numpy as np
from copy import deepcopy

# construct integrator object
integrator = PyQInt()

# build cgf for hydrogen separated by 1.4 a.u.
cgf1 = cgf([0.0, 0.0, 0.0])

cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

# create a copy of the CGF
cgf2 = deepcopy(cgf1)
cgf2.p[2] = 1.4

# Build nuclear attraction integrals
V1 = np.zeros((2,2))
V1[0,0] = integrator.nuclear(cgf1, cgf1, cgf1.p, 1)
V1[0,1] = V1[1,0] = integrator.nuclear(cgf1, cgf2, cgf1.p, 1)
V1[1,1] = integrator.nuclear(cgf2, cgf2, cgf1.p, 1)

V2 = np.zeros((2,2))
V2[0,0] = integrator.nuclear(cgf1, cgf1, cgf2.p, 1)
V2[0,1] = V2[1,0] = integrator.nuclear(cgf1, cgf2, cgf2.p, 1)
V2[1,1] = integrator.nuclear(cgf2, cgf2, cgf2.p, 1)

# print result
print(V1,V2)
```

### Two-electron integrals

```python
import pyqint
import numpy as np
from copy import deepcopy

# construct integrator object
integrator = PyQInt()

# build cgf for hydrogen separated by 1.4 a.u.
cgf1 = cgf([0.0, 0.0, 0.0])

cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)

# create a copy of the CGF
cgf2 = deepcopy(cgf1)
cgf2.p[2] = 1.4

T1111 = integrator.repulsion(cgf1, cgf1, cgf1, cgf1)
T1122 = integrator.repulsion(cgf1, cgf1, cgf2, cgf2)
T1112 = integrator.repulsion(cgf1, cgf1, cgf1, cgf2)
T2121 = integrator.repulsion(cgf2, cgf1, cgf2, cgf1)
T1222 = integrator.repulsion(cgf1, cgf2, cgf2, cgf2)
T2211 = integrator.repulsion(cgf2, cgf2, cgf1, cgf1)

print(T1111)
print(T1122)
print(T1112)
print(T2121)
print(T1222)
print(T2211)
```
