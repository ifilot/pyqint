PyQInt: a Python package for Gaussian integrals
===============================================

.. image:: https://img.shields.io/github/v/tag/ifilot/pyqint?label=version
   :alt: GitHub tag (latest SemVer)
.. image:: https://github.com/ifilot/pyqint/actions/workflows/build.yml/badge.svg
   :target: https://github.com/ifilot/pyqint/actions/workflows/build.yml
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

:program:`PyQInt` is a Python package for calculating one- and two-electron integrals as
encountered in electronic structure calculations. Since integral evaluation can
be quite computationally intensive, they are programmed in C++ and connected to
Python using Cython.

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
