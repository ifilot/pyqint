.. _installation:
.. index:: Installation

Installation
============

:program:`PyQInt` is available via PyPI.

Why environments?
-----------------
An environment keeps the packages for a project separate from the rest of your
system (and from other projects). This prevents accidental upgrades or conflicts
that can break existing setups. It also makes it easy to remove the environment
later without touching your system-wide Python. For these reasons, installing
into Conda's ``base`` environment is discouraged.

PyPI (recommended via virtual environment)
----------------------------------------------------------

Use Python's built-in virtual environments to isolate your installation::

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