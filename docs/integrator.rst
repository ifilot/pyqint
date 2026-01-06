.. index:: integrator

Integrator object
=================

The :code:`PyQInt` class is the low-level numerical integrator used internally
throughout the ``pyqint`` package. It provides efficient evaluation of
Gaussian-type orbital (GTO) and contracted Gaussian function (CGF) integrals
via a C++ backend with optional OpenMP parallelization.

Most users do **not** interact with :code:`PyQInt` directly. Instead, higher-
level routines construct and manage an integrator instance internally. Direct
use of :code:`PyQInt` is intended for advanced users who need fine-grained
control over integral evaluation or wish to experiment with custom workflows.

Overview
--------

The integrator implements core quantum-chemical building blocks, including:

* Overlap integrals
* Kinetic energy integrals
* Nuclear attraction integrals
* Two-electron (electron–repulsion) integrals
* First derivatives with respect to nuclear coordinates
* Grid-based evaluation of basis functions and molecular orbitals

The implementation is optimized for repeated evaluation and makes use of
internal caching where possible.

Creating an Integrator
----------------------

An integrator instance is created by specifying limits on angular momentum
and the maximum order of the Boys function:

.. code-block:: python

   from pyqint import PyQInt

   qint = PyQInt(lmax=4, nu_max=12)

The parameters control internal cache sizes and directly affect both memory
usage and performance. Reusing a single integrator instance is recommended
when evaluating many integrals.

Typical Usage
-------------

Most integrals are evaluated using contracted Gaussian functions (CGFs):

.. code-block:: python

   S = qint.overlap(cgf1, cgf2)
   T = qint.kinetic(cgf1, cgf2)
   V = qint.nuclear(cgf1, cgf2, rc, zc)

Two-electron integrals can be computed either individually or in batches:

.. code-block:: python

   eri = qint.repulsion(cgf1, cgf2, cgf3, cgf4)

For performance-critical workflows, the integrator provides routines that
build full matrices and tensors in a single call using OpenMP:

.. code-block:: python

   S, T, V, tei = qint.build_integrals_openmp(cgfs, nuclei)

Parallel Execution
------------------

If OpenMP support was available at build time, the integrator will
automatically parallelize suitable workloads. The number of threads is
controlled via standard OpenMP environment variables such as ``OMP_NUM_THREADS``.

The current thread count can be queried at runtime:

.. code-block:: python

   nthreads = qint.get_num_threads()

Build Information
-----------------

For diagnostic and reproducibility purposes, the integrator exposes
build-time information:

.. code-block:: python

   info = qint.get_compile_info()

The returned dictionary includes compiler details, OpenMP availability,
and compilation date and time.

Grid-Based Evaluation
---------------------

The integrator also supports evaluation of basis functions, molecular
orbitals, and their gradients on three-dimensional grids:

.. code-block:: python

   grid = qint.build_rectgrid3d(-5.0, 5.0, 50)
   values = qint.plot_wavefunction(grid, coeff, cgfs)

These routines are primarily intended for visualization and analysis.

Notes
-----

* The :code:`PyQInt` object is **not thread-safe** across Python threads.
* Creating multiple integrators increases memory usage due to internal caches.
* Most users should rely on higher-level interfaces provided by ``pyqint``.
