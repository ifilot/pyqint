Integral evaluation acceleration routines
=========================================

.. contents:: Table of Contents
    :depth: 3

Two-electron integrals are evaluated in :program:`PyQInt` using a
Hellsing-type expansion combined with a tabulated evaluation of the
Boys function. To reduce the computational cost associated with repeated
integral evaluations, both the Hellsing kernels and the Boys function
values are cached.

This section describes how these caches are constructed, how they are
used during integral evaluation, and how their behavior can be controlled
by the user.

.. note::

    This section is intended for advanced users familiar with molecular
    integral evaluation in quantum chemistry. For most users, the default
    settings are sufficient and require no modification.

Overview
--------

The :code:`PyQInt` integrator constructs internal lookup tables for

* one-dimensional Hellsing expansion kernels
* Boys function values up to a maximum order.

These tables are created during construction of the integrator object
and are reused across all subsequent integral evaluations.

Hellsing Kernel Cache
---------------------

The Hellsing expansion expresses Cartesian electron-electron repulsion integrals
as finite sums over one-dimensional kernels. 

.. math::

    \begin{aligned}
        K^{(1\mathrm{D})}_{l_1 l_2 l_3 l_4}
        &=
        \sum_{\substack{
        i_1,i_2,i_3,i_4 \\
        o_1,o_2,o_3,o_4 \\
        r_1,r_2,u
        }}
        (-1)^{l_1+l_2+o_2+r_1+o_3+r_2+u}
        \, l_1! \, l_2! \, l_3! \, l_4!
        \nonumber\\[0.6em]
        &\quad\times
        \frac{(o_1+o_2)!}
            {4^{\,i_1+i_2+r_1}
            \, i_1! \, i_2!
            \, o_1! \, o_2!
            \, r_1!
            \, (l_1-2i_1-o_1)!
            \, (l_2-2i_2-o_2)!
            \, (o_1+o_2-2r_1)!}
        \nonumber\\[0.8em]
        &\quad\times
        \frac{(o_3+o_4)!}
            {4^{\,i_3+i_4+r_2}
            \, i_3! \, i_4!
            \, o_3! \, o_4!
            \, r_2!
            \, (l_3-2i_3-o_3)!
            \, (l_4-2i_4-o_4)!
            \, (o_3+o_4-2r_2)!}
        \nonumber\\[0.8em]
        &\quad\times
        \frac{\mu!}
            {4^{\,u}
            \, u!
            \, (\mu-2u)!} \, .
    \end{aligned}

with

.. math::

    \begin{aligned}
        \mu &= l_1 + l_2 + l_3 + l_4
            - 2(i_1+i_2+i_3+i_4)
            - (o_1+o_2+o_3+o_4),
        \\
        0 &\le i_k \le \left\lfloor \frac{l_k}{2} \right\rfloor,
        \qquad
        0 \le o_k \le l_k - 2 i_k,
        \\
        0 &\le r_1 \le \left\lfloor \frac{o_1+o_2}{2} \right\rfloor,
        \qquad
        0 \le r_2 \le \left\lfloor \frac{o_3+o_4}{2} \right\rfloor,
        \\
        0 &\le u \le \left\lfloor \frac{\mu}{2} \right\rfloor.
    \end{aligned}

For a given angular momentum quartet :math:`(l_1, l_2, l_3, l_4)`, the
corresponding kernel depends only on these angular momenta and not on the
molecular geometry. For this reason, all required Hellsing kernels are
precomputed and stored in a lookup table indexed by the angular momenta. The
size of this table is controlled by the parameter :math:`l_{\max}`, which
defines the maximum angular momentum supported by the cache.

Boys Function Lookup Table and Recurrence
-----------------------------------------

The Boys function

.. math::

    F_\nu(T) = \int_0^1 t^{2\nu} e^{-T t^2} \, \mathrm{d}t

appears ubiquitously in the evaluation of electron–electron repulsion
integrals. Direct numerical evaluation of this integral is expensive
when performed repeatedly for many values of the order :math:`\nu`
and argument :math:`T`. To reduce this cost, :program:`PyQInt` employs
a hybrid strategy combining a precomputed lookup table with stable
recurrence relations and exact fallback evaluation.

Lookup Table Construction
~~~~~~~~~~~~~~~~~~~~~~~~~

Upon construction of the integrator, the Boys function is tabulated on a
uniform grid :math:`T_i = i \, \Delta T` with :math:`\Delta T = 0.05`,
covering :math:`0 \leq T < T_{\mathrm{table}} = 40`. For every grid point,
the orders :math:`\nu = 0, \dots, \nu_{\max} + 7` are stored.

The table entries are computed from the series of positive terms

.. math::

    F_\nu(T)
    =
    e^{-T}
    \sum_{k=0}^{\infty}
    \frac{(2T)^k}{(2\nu+1)(2\nu+3)\cdots(2\nu+2k+1)},

which converges for all :math:`T` without cancellation, for the highest
order only; the lower orders follow from the downward recurrence given below,
which is numerically stable.

Taylor Expansion
~~~~~~~~~~~~~~~~

Because the derivative of the Boys function is again a Boys function,

.. math::

    \frac{\mathrm{d}F_\nu(T)}{\mathrm{d}T} = -F_{\nu+1}(T),

the value at an arbitrary :math:`T` follows from a Taylor expansion around
the nearest grid point :math:`T_i`:

.. math::

    F_\nu(T)
    =
    \sum_{k=0}^{7}
    \frac{F_{\nu+k}(T_i) \, (T_i - T)^k}{k!}.

As :math:`|T - T_i| \leq \Delta T / 2`, the truncation error of these eight
terms is below :math:`10^{-16}` relative to :math:`F_\nu(T)`, such that the
Boys function is obtained to machine precision at the cost of a single table
lookup and a short polynomial evaluation.

Recursive Evaluation of Boys Function Blocks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In practical integral evaluation, entire blocks of Boys function values

.. math::

    \{ F_0(T), F_1(T), \dots, F_{\nu_{\max}}(T) \}

are required. To compute these efficiently, :program:`PyQInt` evaluates only
the highest order and obtains the others by the downward recurrence

.. math::

    F_{\nu-1}(T)
    =
    \frac{2T F_\nu(T) + e^{-T}}{2\nu - 1},

which is stable for all :math:`T`, as every term is positive.

For large arguments, :math:`T \geq T_{\mathrm{table}}` and
:math:`T > 2\nu_{\max}`, the value of :math:`F_0(T)` is computed from the
error function,

.. math::

    F_0(T) = \frac{1}{2} \sqrt{\frac{\pi}{T}} \, \mathrm{erf}\left(\sqrt{T}\right),

and the higher orders are obtained via upward recurrence:

.. math::

    F_{\nu+1}(T)
    =
    \frac{(2\nu+1) F_\nu(T) - e^{-T}}{2T}.

This recurrence subtracts two numbers of similar size and is therefore only
used where :math:`T` is large compared to :math:`\nu`.

Exact and Fallback Evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the requested order exceeds the tabulated maximum, or if :math:`T`
lies outside the tabulated range, the Boys function is evaluated without the
table: by the series of positive terms for :math:`T < 2\nu + 40`, and by the
error function and upward recurrence otherwise.

This layered strategy yields the Boys function to near machine precision
across the full range of parameters encountered during integral evaluation.

Integrator Construction and Defaults
------------------------------------

Upon construction of a :code:`PyQInt` integrator object, default cache
sizes are selected:

* :math:`l_{\max} = 4`
* :math:`\nu_{\max} = 12`

These values are sufficient for basis sets composed of first- and second-row
elements. If integrals are requested that exceed the current cache limits, the
integrator will automatically expand the Hellsing kernel cache or fall back to
direct Boys function evaluation.

Advanced users may override the default cache sizes by explicitly specifying
:math:`l_{\max}` and :math:`\nu_{\max}` when constructing the integrator.
Increasing these parameters may improve performance for systems with
high-angular-momentum basis functions at the cost of increased memory usage.

.. code:: python

    from pyqint import PyQInt

    # construct integrator class
    integrator = PyQInt(lmax=6, nu_max=12)

    # use it
    # ...

All wrapper classes, e.g. :code:`HF`, :code:`GeometryOptimization`, etc. will
automatically determine optimal caching sizes and adjust the Kernel cache if
needed.