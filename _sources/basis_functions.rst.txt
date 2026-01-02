.. index:: basis-functions

Basis functions
###############

.. contents:: Table of Contents
    :depth: 3

Gaussian Type Orbitals (GTOs)
=============================

:program:`PyQInt` employs *Cartesian Gaussian Type Orbitals* (GTOs), defined as

.. math::

    \Phi(\alpha, l, m, n, \vec{R}) =
    N (x - X)^{l} (y - Y)^{m} (z - Z)^{n}
    \exp \left(- \alpha \lvert \vec{r} - \vec{R} \rvert^{2} \right)

Here, :math:`\alpha` denotes the orbital exponent,
:math:`\vec{R} = (X, Y, Z)` the position of the orbital center,
:math:`(l, m, n)` the orders of the Cartesian polynomial prefactor,
and :math:`N` a normalization constant chosen such that

.. math::

    \langle \Phi \mid \Phi \rangle = 1

GTO Construction
****************

Gaussian Type Orbitals form the fundamental building blocks of *Contracted
Gaussian Functions* (CGFs; see below). In typical workflows, users interact
with CGFs rather than individual GTOs. Nevertheless, GTOs can be constructed
explicitly when needed via::

    from pyqint import GTO

    coeff = 1.0      # linear expansion coefficient
    alpha = 0.5      # exponent
    l, m, n = 0, 0, 0
    position = (0.0, 0.0, 0.0)

    gto = GTO(coeff, position, alpha, l, m, n)

.. note::

    - The normalization constant :math:`N` is computed automatically by
      :program:`PyQInt` based on the values of :math:`\alpha` and
      :math:`(l, m, n)` and must not be supplied explicitly.
    - When working with individual GTOs (i.e. not as part of a CGF), the
      coefficient should be set to ``1.0``. This parameter represents the
      linear expansion coefficient used internally by CGFs.

Retrieving GTO Data Members
***************************

The :code:`GTO` class is intentionally designed as a *flat and accessible*
Python object. All defining parameters of the orbital are stored as public
data members and can be accessed directly.

The following attributes are available:

- ``c`` — linear expansion coefficient
- ``p`` — orbital center :math:`(X, Y, Z)`
- ``alpha`` — Gaussian exponent
- ``l``, ``m``, ``n`` — Cartesian polynomial orders

For example::

    coeff = gto.c
    position = gto.p
    alpha = gto.alpha
    l = gto.l
    m = gto.m
    n = gto.n

    print("Coefficient:", coeff)
    print("Position:", position)
    print("Exponent:", alpha)
    print("Orders:", (l, m, n))

These attributes fully define the mathematical form of the primitive Gaussian
orbital and may be inspected or reused when constructing higher-level objects,
such as CGFs.

Evaluating a GTO
****************

In addition to direct data access, a GTO provides methods for numerical
evaluation. For example, the orbital amplitude at a point
:math:`(x, y, z)` can be obtained using::

    value = gto.get_amp(x, y, z)

The normalization constant used internally may be retrieved via::

    norm = gto.get_norm()

This separation between *data members* and *evaluation routines* allows users
to both inspect the orbital parameters and efficiently compute values when
needed.

Contracted Gaussian Functions (CGF)
===================================

Several GTOs can be combined to produce a so-called Contracted Gaussian Functional which
is esentially a linear combination of GTOs as given by

.. math::

    \phi = \sum_{i} c_{i} \Phi_{i}(\alpha,l,m,n,\vec{R})

To build a CGF, we first have to produce the CGF object and then
add GTOs to it::

    from pyqint import PyQInt, CGF

    cgf = CGF([0.0, 0.0, 0.0])

    cgf.add_gto(0.154329, 3.425251, 0, 0, 0)
    cgf.add_gto(0.535328, 0.623914, 0, 0, 0)
    cgf.add_gto(0.444635, 0.168855, 0, 0, 0)

.. note::
    The first argument of the :code:`add_gto` function is the linear expansion coefficient
    :math:`c_{i}` and the second argument is :math:`\alpha`.

Position Handling and Advanced Usage
************************************

A :code:`CGF` object carries an explicit spatial position :math:`\vec{R}`,
which is specified when the object is constructed. This position is used as
the *default center* for all Gaussian Type Orbitals (GTOs) added to the CGF via
the standard :code:`add_gto` method.

When using

.. code-block:: python

    cgf.add_gto(c, alpha, l, m, n)

the position of the CGF itself is implicitly copied to the newly created GTO.
This is the typical and intended usage pattern for atom-centered basis sets,
where all primitive Gaussians contributing to a contracted function are located
on the same atomic center. Under this assumption, the CGF represents a single
atomic basis function composed of multiple primitives sharing a common center.

Non-Atom-Centered GTOs and Mixed-Center CGFs
********************************************

In more advanced scenarios, it may be desirable to construct contracted
functions composed of GTOs located at *different spatial positions*. This can
arise, for example, when building symmetry-adapted basis functions or other
non-standard orbital constructions.

For such cases, the :code:`CGF` class provides the method
:code:`add_gto_with_position`, which allows explicit specification of the GTO
center:

.. code-block:: python

    cgf.add_gto_with_position(c, p, alpha, l, m, n)

where :code:`p` is a three-component array specifying the Cartesian coordinates
of the GTO center.

When this method is used, the position stored in the CGF object itself is
*ignored* for that particular GTO, and the explicitly provided position is
stored with the primitive instead.

Integral Evaluation and Position Resolution
*******************************************

It is important to emphasize that **all integral evaluations are performed using
the positions of the individual GTOs**, not the position stored in the CGF
object. The CGF position serves only as a convenient default when adding
atom-centered primitives via :code:`add_gto`. Consequently, CGFs containing GTOs
at mixed centers are handled correctly by the integral engine, provided that
:code:`add_gto_with_position` is used where appropriate.

Users are therefore encouraged to:

- Use :code:`add_gto` for standard atom-centered basis functions
- Use :code:`add_gto_with_position` when constructing non-standard or
  symmetry-adapted contracted functions
- Rely on GTO positions - **not CGF positions** — when reasoning about integral
  evaluation
