.. _installation:
.. index:: Installation

Installation
============

.. tip::
    For Windows users with relatively little experience with Python, we warmly
    recommend the `Anaconda distribution <https://www.anaconda.com/products/distribution>`_.
    Anaconda provides Python, an integrated desktop environment (Spyder), and
    widely used packages such as NumPy and Matplotlib.

:program:`PyQInt` is available via both Conda and PyPI. We recommend **using an
isolated environment** rather than installing into the global Python or Conda
``base`` environment. This avoids package conflicts and makes your setup easier
to reproduce.

Why environments?
-----------------
An environment keeps the packages for a project separate from the rest of your
system (and from other projects). This prevents accidental upgrades or conflicts
that can break existing setups. It also makes it easy to remove the environment
later without touching your system-wide Python. For these reasons, installing
into Conda’s ``base`` environment is discouraged.

Windows / Conda (recommended)
-----------------------------

Create and activate a fresh Conda environment, then install :program:`PyQInt`::

    conda create -n eoesc-env -c ifilot pyqint
    conda activate eoesc-env

.. note::
    If you would like to use :program:`PyQInt`’s **isosurface functionality**,
    you will also need the optional package ``pytessel``. You can install it
    together with :program:`PyQInt`::

        conda create -n eoesc-env -c ifilot pyqint pytessel
        conda activate eoesc-env

    or add it later with::

        conda install -c ifilot pytessel

Optionally install the Spyder IDE and common scientific packages::

    conda install spyder matplotlib scipy pandas openpyxl

Linux / macOS / PyPI (recommended via virtual environment)
----------------------------------------------------------

Use Python’s built-in virtual environments to isolate your installation::

    python3 -m venv venv
    source venv/bin/activate
    pip install --upgrade pip
    pip install pyqint

Optional extra for **isosurface functionality**::

    pip install pytessel

To leave the environment, run::

    deactivate

Alternative (single-user install without a virtual environment)
---------------------------------------------------------------

If you prefer not to create a virtual environment, you can install for the
current user only (no system-wide changes)::

    pip install --user pyqint

.. warning::
    We **do not recommend** using ``sudo pip install ...``. Installing packages
    with administrative privileges can overwrite or conflict with system Python
    components and may break tools your operating system relies on.

Troubleshooting
---------------

- If Conda reports conflicts during installation, it usually means the active
  environment has incompatible packages. Creating a fresh environment (as shown
  above) is the most reliable fix.

- If ``pip`` cannot find the package or fails due to permissions, ensure your
  virtual environment is activated, or use ``--user`` for a local user install.

- If Spyder does not see the environment, launch Spyder **from within** the
  activated environment (``conda activate eoesc-env`` then ``spyder``), or
  configure the interpreter path in Spyder’s preferences.