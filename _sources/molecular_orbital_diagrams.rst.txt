.. index:: molecular-orbital-diagrams

Molecular Orbital Diagrams
==========================

While :program:`PyQInt` provides access to molecular orbital energies and
coefficients, the construction of molecular orbital (MO) diagrams is handled
externally using the :program:`PyMoDia` package. PyMoDia is a dedicated tool for
building clear and customizable MO diagrams from quantum-chemical data. The
PyMoDia package is documented separately at: https://ifilot.github.io/pymodia

In this section, we demonstrate how results obtained from :program:`PyQInt` can be
exported and visualized using PyMoDia.

.. warning::

    PyMoDia currently only support **restricted** Hartree-Fock calculations.

Overview
--------

PyMoDia represents molecular orbital diagrams as energy-level schemes,
optionally including fragment orbitals, interaction arrows, and symmetry-adapted
labels. To simplify interoperability, PyMoDia provides helper routines that can
automatically extract the relevant information from :program:`PyQInt`
Hartree-Fock results.

:program:`PyMoDia` can be easily installed via

.. code:: bash

    pip install pymodia

Canonical Molecular Orbital Diagram for CO
-------------------------------------------

The following example constructs a canonical molecular orbital diagram for
carbon monoxide (CO) using Hartree-Fock/STO-3G orbitals computed with
:program:`PyQInt`.

.. code-block:: python

    import os
    from pymodia import MoDia, MoDiaData, autobuild_from_pyqint, MoDiaSettings
    from pyqint import MoleculeBuilder, HF, FosterBoys

    # Perform PyQInt calculations for CO
    mol = MoleculeBuilder.from_name('co')
    res = HF(mol, 'sto3g').rhf()

    # Adjust PyMoDia visualization settings
    settings = MoDiaSettings()
    settings.orbc_color = '#555555'
    settings.arrow_color = '#CC0000'

    # Automatically construct molecule and fragment definitions
    mol, f1, f2 = autobuild_from_pyqint(res, name='co')

    # Retrieve molecular orbital energies
    moe = res['orbe']

    # Slightly shift one orbital level to improve visual separation
    moe[6] += 0.1

    # Build PyMoDia data object
    data = MoDiaData(mol, f1, f2)
    data.set_moe(moe)

    # Construct and export the MO diagram
    diagram = MoDia(
        data,
        draw_level_labels=True,
        level_labels_style='mo_ao',
        mo_labels=['1σ', '2σ', '3σ', '4σ', '1π', '1π', '5σ', '2π', '2π', '6σ'],
        settings=settings
    )

    diagram.export_svg(
        os.path.join(os.path.dirname(__file__), "mo_co_canonical.svg")
    )

.. figure:: _static/img/mo_co_canonical.svg

    Molecular orbital diagram for the canonical molecular orbitals of CO.

.. note::

   Atomic orbital (AO) energy levels shown in PyMoDia diagrams are obtained from
   an internal element-specific dictionary and are intended for qualitative
   visualization purposes only. They do not necessarily represent the true
   atomic orbital energies for a given calculation and may be adjusted manually
   (see the PyMoDia manual). Quantitative AO energies can, in principle, be
   obtained from unrestricted Hartree–Fock calculations on the isolated atoms,
   which are not yet supported in :program:`PyQInt`.

.. tip::

    For complex molecules or closely spaced orbital manifolds, it is sometimes
    necessary to apply small *visual offsets* to selected molecular orbital
    energy levels in order to avoid overlaps in the diagram. In the example
    above, the energy of the :math:`5\sigma` orbital is shifted slightly upward
    to separate it from the degenerate :math:`2\pi` levels.

    It is important to emphasize that:

    - These adjustments affect **only the vertical placement** of the orbital
      levels in the diagram
    - The **numerical energy values displayed** in the figure correspond to the
      *original, unmodified* Hartree-Fock orbital energies
    - The underlying quantum-chemical results remain unchanged

    Such adjustments are therefore purely cosmetic and serve only to improve the
    clarity and readability of the resulting MO diagram.

Localized Orbital Diagrams
--------------------------

:program:`PyMoDia` can equally well be used to visualize localized molecular
orbitals, such as those obtained via Foster-Boys localization. Since orbital
localization corresponds to a unitary transformation of the occupied subspace,
the localized orbitals can be passed directly to :program:`PyMoDia` using the
same workflow.

Localized MO diagrams often provide a more chemically intuitive picture, making
bonding patterns and fragment interactions easier to interpret, while retaining
a direct connection to the underlying :program:`PyQInt` calculation.

.. code-block:: python

    import os
    from pymodia import MoDia, MoDiaData, autobuild_from_pyqint, MoDiaSettings
    from pyqint import MoleculeBuilder, HF, FosterBoys

    # Perform PyQInt calculations for CO and its localization
    mol = MoleculeBuilder.from_name('co')
    res = HF(mol, 'sto3g').rhf()

    # adjust settings
    settings = MoDiaSettings()
    settings.orbc_color = '#555555'
    settings.arrow_color = '#CC0000'

    # making diagram for localized orbitals
    resfb = FosterBoys(res).run()
    resfb['nuclei'] = res['nuclei'] # no longer required from PyQInt >= 1.2.0
    mol, f1, f2 = autobuild_from_pyqint(resfb, name='co')

    # we make here a small adjustment to the height of the third orbital to avoid
    # overlap with the triple degenerate state of the localized MOs of CO
    moe = resfb['orbe']
    moe[2] -= 0.5
    moe[6] += 0.1

    # build data object
    data = MoDiaData(mol, f1, f2)
    data.set_moe(moe)

    # construct diagram
    labels = [''] * len(resfb['orbe'])
    diagram = MoDia(data, draw_level_labels=True, level_labels_style='mo_ao',
                    mo_labels=labels,
                    settings=settings)
    diagram.export_svg(os.path.join(os.path.dirname(__file__), "mo_co_localized.svg"))

.. figure:: _static/img/mo_co_localized.svg

    Molecular orbital diagram for the localized molecular orbitals of CO. Note
    that we observe a threefold degenerate state correspond to the C-O triple
    bond.