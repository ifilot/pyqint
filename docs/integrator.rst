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

The ``eri_engine`` keyword selects the electron-repulsion kernel for primitive
integrals, contracted integrals, batch builds, and their geometric derivatives.
The default is ``"hellsing"`` (cached Hellsing). Use ``"rys"`` to opt into Rys
quadrature. The reference variants are ``"hellsing_direct"``,
``"huzinaga"`` (block Boys functions), and ``"huzinaga_scalar"``.

Use ``eri_engine="auto"`` to opt into the ERI orchestrator. It records primitive
counts and verifies shared centers from the primitive coordinates once per
basis build. Batch ERI builds reuse prepared Gaussian products for eligible
Hellsing quartets. For verified Cartesian contractions, the orchestrator
skips odd-parity integrals when all four centers lie on the same Cartesian
reflection plane. It selects Rys for (ss|ss), for sufficiently large batches
of s/p shell components, and for quartets with all functions d-type or higher
along the same Cartesian axis. Verified shell components share roots,
weights, and axis recurrences.
Geometric derivatives and standalone ERIs use cached Hellsing in automatic
mode.

Call ``build_tei_openmp(cgfs, profile=True)`` to get the tensor and a dictionary
with times in seconds for pair setup, quartet planning, integral evaluation,
dense tensor assembly, and NumPy conversion, along with the route counts.

.. code-block:: python

   qint = PyQInt(eri_engine="rys")

Typical Usage
-------------

Most integrals are evaluated using contracted Gaussian functions (CGFs):

.. code-block:: python

   S = qint.overlap(cgf1, cgf2)
   T = qint.kinetic(cgf1, cgf2)
   V = qint.nuclear(cgf1, cgf2, rc, zc)

Use ``qint.build_dipole_openmp(cgfs, origin=(0.0, 0.0, 0.0))`` for the
three Cartesian position integral matrices in bohr. The returned array
has shape ``(3, N, N)``.

Two-electron integrals can be computed either individually or in batches:

.. code-block:: python

   eri = qint.repulsion(cgf1, cgf2, cgf3, cgf4)

For performance-critical workflows, the integrator provides routines that
build full matrices and tensors in a single call using OpenMP:

.. code-block:: python

   S, T, V, tei = qint.build_integrals_openmp(cgfs, nuclei)

Pass ``packed=True`` to ``build_integrals_openmp`` or ``build_tei_openmp``
to obtain only the unique ERIs. The result is a one-dimensional array in
``teindex`` order. For example, ``tei[qint.teindex(i, j, k, l)]`` returns
the integral ``(ij|kl)`` for any ordering of the four indices.

.. code-block:: python

   tei = qint.build_tei_openmp(cgfs, packed=True)
   J, K = qint.contract_packed_tei(tei, density)
   J, Ka, Kb = qint.contract_packed_tei(tei, density_alpha, density_beta)

With two spin densities, the Coulomb matrix uses their sum and each
exchange matrix uses its corresponding spin density.

Memory preflight
----------------

Bulk ERI and geometric derivative builds check their expected memory use
against currently available RAM. On Linux, the check also considers readable
cgroup limits. It accounts for a process address-space limit where available.
If the estimate exceeds the available memory, PyQInt raises ``MemoryError``
with the estimated requirement. Automatic ERI routing also checks its
workspace as the quartet plan grows.

The check uses a snapshot of available memory. When no memory figures are
available, the build proceeds using the system allocator. Selecting packed
storage reduces the ERI memory requirement for energy calculations.
Hartree–Fock jobs that request forces also check their derivative workspace
before SCF.

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

To evaluate *all* basis functions (or their gradients) on the same grid, for
example to build numerical integration grids for density functional theory,
use the batched routines. These evaluate all basis functions within a single
parallel loop over the grid points, which is considerably faster than calling
:code:`plot_basis_function` once per basis function:

.. code-block:: python

   amps  = qint.plot_basis_functions(grid, cgfs)           # (nbf, npts)
   grads = qint.plot_basis_function_gradients(grid, cgfs)  # (nbf, npts, 3)

Notes
-----

* The :code:`PyQInt` object is **not thread-safe** across Python threads.
* Creating multiple integrators increases memory usage due to internal caches.
* Most users should rely on higher-level interfaces provided by ``pyqint``.
