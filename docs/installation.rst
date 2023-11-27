.. _installation:
.. index:: Installation

Installation
============

.. tip::
    For Windows users with relatively little experience with Python, we warmly
    recommend to use the `Anaconda distribution <https://www.anaconda.com/products/distribution>`_.
    Anaconda is an all-in-one package containing the Python compiler,
    an integrated desktop environment (Spyder) and plenty of useful Python
    packages such as numpy and matplotlib.

:program:`PyQInt` is distributed via both Anaconda package as well as PyPI. For
Windows, it is recommended to install :program:`PyQInt` via Anaconda, while
for Linux, we recommend to use PyPI.

Windows / Anaconda
------------------

To install :program:`PyQInt` under Windows, open an Anaconda Prompt window
and run::

    conda install -c ifilot pyqint

.. note::
    Sometimes Anaconda is unable to resolve the package dependencies. This can
    be caused by a broken environment. An easy solution is to create a new
    environment. See the "Troubleshooting" section at the end of this page
    for more information.

Linux / PyPI
------------

To install :program:`PyQInt` systemwide, run::

    sudo pip install pyqint

or to install :program:`PyQInt` only for the current user, run::

    pip install pyqint

Troubleshooting
---------------

The Anaconda packaging system can sometimes be quite finicky and sometimes
packages conflict with each other. A way to work around this issue is to create
a separate environment and only use that environment for the electronic
resources associated with this project.

To create the new environment (called eoesc-env), run::

    conda create -n eoesc-env python=3.11

Next, open the environment with::

    conda activate eoesc-env

and install the required packages::

    conda install -c ifilot pyqint pytessel

Finally, you can install the IDE Spyder using::

    conda install spyder matplotlib scipy pandas openpyxl
