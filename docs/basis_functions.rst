.. index:: basis-functions

Basis functions
###############

.. contents:: Table of Contents
    :depth: 3

Gaussian type orbitals (GTO)
============================

:program:`PyQInt` uses cartesian Gaussian type orbitals as given by

.. math::

    \Phi(\alpha,l,m,n,\vec{R}) = N (x - X)^{l} (y - Y)^{m} (z - Z)^{n} \exp \left(- \alpha |\vec{r} - \vec{R}|^{2} \right)

wherein :math:`\alpha` is the exponent, :math:`\vec{R} = \left(X,Y,Z\right)` the
position of the orbital, :math:`(l,m,n)` the orders of the pre-exponential
polynomial, and :math:`N` a normalization constant such that

.. math::

    \left< \Phi | \Phi \right> = 1

GTOs are a fundamental building block of CGF (see below) and typically a user
would not directly work with them (a notable exception is provided below).
Nevertheless, GTO objects can be constructed as follows::

    from pyqint import PyQInt, CGF, GTO

    coeff = 1.0    # coefficients only have meaning for GTOs within a CGF
    alpha = 0.5
    l,m,n = 0,0,0
    p = (0,0,0)
    G = GTO(coeff, p, alpha, l, m, n)

.. note::
    * The normalization constant is automatically calculated by `PyQInt` based
      on the value of :math:`\alpha` and :math:`(l,m,n)` and does not have
      to be supplied by the user.
    * If you work with individual GTOs, the first parameter to construct the GTO
      should have a value of 1.0. This first parameter corresponds to the linear
      expansion coefficient used in the formulation of Contracted Gaussian Functions
      (see below).

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
