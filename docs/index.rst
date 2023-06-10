PyQInt: a Python package for Gaussian integrals
===============================================

.. image:: https://img.shields.io/github/v/tag/ifilot/pyqint?label=version
   :alt: GitHub tag (latest SemVer)
.. image:: https://github.com/ifilot/pyqint/actions/workflows/build.yml/badge.svg
   :target: https://github.com/ifilot/pyqint/actions/workflows/build.yml
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

:program:`PyQInt` is a Python package for calculating one- and two-electron
integrals as encountered in electronic structure calculations. Since
integral evaluation can be quite computationally intensive, the evaluation
is programmed in C++ and connected to Python using Cython.

:program:`PyQInt` mainly serves as an educational package to teach students
how to perform (simple) electronic structure calculations wherein the most
difficult task, i.e. the integral evaluation, is already encapsulated in
a handy set of routines. With :program:`PyQInt`, the student can for example
build their own Hartree-Fock routine. Some common electronic structure routine,
most notably the Hartree-Fock algorithm, is also readily available.

:program:`PyQInt` offers supporting scripts for facile visualization of result
such as producing contour plots for the molecular orbitals. Below, an example
is shown for the molecular orbitals of the CO molecule.

.. figure:: _static/img/co.jpg

    Canonical molecular orbitals of CO visualized using contour plots.

:program:`PyQInt` not only supports calculation of the canonical molecular orbitals
via the (restricted) Hartee-Fock procedure, but can also be used to construct
the localized molecular orbitals which is relevant for showing the similarity
between modern electronic structure methods and classical Lewis theory. In the
image below, one can observe the **canonical** molecular orbitals for the CO molecule
as well as the **localized** molecular orbitals.

.. figure:: _static/img/co_isosurfaces.jpg

    Canonical molecular orbitals of CO visualized using isosurfaces with an
    isovalue of +/-0.1.

.. figure:: _static/img/co_isosurfaces_foster_boys.jpg

    Localized molecular orbitals of CO visualized using isosurfaces with an
    isovalue of +/-0.1. Note that the localization procedure has only been
    applied to the occupied molecular orbitals. Observe that the localized
    orbitals contain a triple-degenerate state corresponding to the triple
    bond and two lone pairs for C and O.

:program:`PyQInt` has been developed at the Eindhoven University of Technology,
Netherlands. :program:`PyQInt` and its development are hosted on `github
<https://github.com/ifilot/pyqint>`_.  Bugs and feature
requests are ideally submitted via the `github issue tracker
<https://github.com/ifilot/pyqint/issues>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   user_interface
   community_guidelines

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
