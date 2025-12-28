.. index:: population-analysis

Population Analysis
===================

.. contents:: Table of Contents
    :depth: 3

.. warning::

    Population Analysis methods currently only support **restricted**
    Hartree-Fock calculations.

Single Atom Population Analysis
-------------------------------

Single-atom population analysis aims to assign an effective electronic
population (and corresponding atomic charge) to each atom in a molecule by
partitioning the total electron density among atom-centered basis functions.
While such atomic charges are not quantum-mechanical observables, they are
widely used as qualitative indicators of charge transfer, polarity, and
chemical bonding.

In :program:`PyQInt`, two closely related single-atom population analysis
techniques are implemented:

* **Mulliken Population Analysis**
* **Löwdin Population Analysis**

Both methods operate on the one-particle density matrix expressed in a
localized atomic orbital basis.

Mulliken Population Analysis
****************************

In Mulliken population analysis, the electronic population associated with
atom :math:`A` is obtained by summing the diagonal elements of the
*overlap-weighted density matrix* corresponding to basis functions centered on
that atom.

The overlap-weighted density matrix is defined as

.. math::

    \mathbf{P}^{(\mathrm{M})} = \mathbf{P}\mathbf{S}

where :math:`\mathbf{P}` is the density matrix and :math:`\mathbf{S}` is the
overlap matrix.

The Mulliken population of atom :math:`A` is then given by

.. math::

    N_A^{(\mathrm{M})} =
    \sum_{i \in A} P^{(\mathrm{M})}_{ii}

and the corresponding Mulliken atomic charge is

.. math::

    q_A^{(\mathrm{M})} = Z_A - N_A^{(\mathrm{M})}

where :math:`Z_A` is the nuclear charge of atom :math:`A`.

Mulliken population analysis is straightforward to compute but is known to be
highly sensitive to the choice of basis set, particularly in the presence of
diffuse or highly overlapping basis functions.

Löwdin Population Analysis
**************************

Löwdin population analysis mitigates some of the basis-set dependence of
Mulliken analysis by first transforming the density matrix into an
orthonormalized basis using Löwdin symmetric orthogonalization.

The Löwdin orthogonalization matrix is defined as

.. math::

    \mathbf{X} = \mathbf{S}^{-1/2}

where :math:`\mathbf{S}` is the overlap matrix. The density matrix in the
orthonormalized basis is then given by

.. math::

    \mathbf{P}' = \mathbf{X}\mathbf{P}\mathbf{X}

The Löwdin population of atom :math:`A` is obtained as

.. math::

    N_A^{(\mathrm{L})} =
    \sum_{i \in A} P'_{ii}

with the corresponding Löwdin atomic charge

.. math::

    q_A^{(\mathrm{L})} = Z_A - N_A^{(\mathrm{L})}

Because the Löwdin transformation produces an orthonormal basis, Löwdin
populations are generally more stable with respect to basis-set variations and
are often regarded as more physically meaningful than Mulliken populations.

Example: Atomic Charge Analysis
*******************************

The following example demonstrates Mulliken and Löwdin atomic charge analysis
for methane (CH\ :sub:`4`) using a minimal STO-3G basis set.

.. code-block:: python

    from pyqint import MoleculeBuilder, HF, PopulationAnalysis

    mol = MoleculeBuilder.from_name('CH4')
    res = HF(mol, 'sto3g').rhf()
    pa = PopulationAnalysis(res)

    charges_mulliken = [pa.mulliken(n) for n in range(len(mol))]
    charges_lowdin   = [pa.lowdin(n)   for n in range(len(mol))]

    print('Symbol  Lowdin        Mulliken')
    for a, cl, cm in zip(mol, charges_lowdin, charges_mulliken):
        print('%2s  %12.6f  %12.6f' % (a[0], cl, cm))

Example output::

    Symbol  Lowdin        Mulliken
     C     -0.134084     -0.248559
     H      0.033521      0.062140
     H      0.033521      0.062140
     H      0.033521      0.062140
     H      0.033521      0.062140

In both population analysis schemes, carbon carries a net negative charge,
while the four hydrogen atoms carry equal positive charges, reflecting the
tetrahedral symmetry of the molecule. The Löwdin charges are smaller in
magnitude than the Mulliken charges, consistent with the reduced basis-set
sensitivity of Löwdin population analysis.

For neutral molecules, the sum of all atomic charges obtained from either scheme
is exactly zero. 

Molecular Orbital Hamilton and Overlap Population Analysis
----------------------------------------------------------

.. note::

    In earlier versions of :program:`PyQInt`, this functionality was
    erroneously referred to as *Crystal Orbital Hamilton Population* (COHP).
    While the underlying idea is closely related, COHP is formally defined
    for *crystal orbitals* (Bloch functions) in periodic systems.

Background of MOHP, MOOP, and MOBI
**********************************

Within the scope of chemical bonding analysis, molecular orbitals can be
classified as bonding, antibonding, or non-bonding with respect to any given
pair of atoms. When working with localized basis functions, this classification
can be made explicit by projecting molecular orbitals onto atomic subspaces.

Two closely related population analysis techniques are implemented in *pyqint*:

* **MOHP** - Molecular Orbital *Hamilton* Population
* **MOOP** - Molecular Orbital *Overlap* Population
* **MOBI** - Molecular Orbital *Bond* Index

All three analyses quantify the contribution of a given molecular orbital to the
interaction between two atoms.

In a localized basis representation, the **MOHP coefficient** of molecular
orbital :math:`k` with respect to atoms :math:`A` and :math:`B` is defined as

.. math::

    \mathrm{MOHP}_k =
    2 \sum_{i \in A} \sum_{j \in B}
    C_{ik} \, H_{ij} \, C_{jk}

where:

- :math:`C_{ik}` and :math:`C_{jk}` are elements of the molecular orbital
  coefficient matrix :math:`\mathbf{C}`
- :math:`H_{ij}` is an element of the Fock (Hamiltonian) matrix
- the factor of 2 accounts for spin degeneracy in restricted Hartree-Fock theory

Analogously, the **MOOP coefficient** is defined as

.. math::

    \mathrm{MOOP}_k =
    2 \sum_{i \in A} \sum_{j \in B}
    C_{ik} \, S_{ij} \, C_{jk}

where :math:`S_{ij}` is an element of the overlap matrix :math:`\mathbf{S}`.

.. note::

    Both MOHP and MOOP can be evaluated for virtual (unoccupied) orbitals.
    However, the interpretation of such values should be made with caution,
    as virtual orbitals do not correspond to occupied electronic states.

Lastly, the **MOBI coefficient** is defined as

.. math::

    \mathrm{MOBI}_k =
    2 \sum_{i \in A} \sum_{j \in B}
    C_{ik} \, P_{ij} \, C_{jk}

where :math:`P_{ij}` is an element of the density matrix :math:`\mathbf{P}`.

.. note::

    MOHP, MOOP, and MOBI can be evaluated for virtual (unoccupied) orbitals.
    However, the interpretation of such values should be made with caution,
    as virtual orbitals do not correspond to occupied electronic states.

Procedure of MOHP, MOOP, and MOBI
*********************************

MOHP, MOOP, and MOBI calculations are performed using the `PopulationAnalysis` class, 
which takes the output of a Hartree-Fock calculation as input.

The example below demonstrates MOHP, MOOP, and MOBI analysis for the CO molecule.

.. code-block:: python

    from pyqint import MoleculeBuilder, HF, PopulationAnalysis

    mol = MoleculeBuilder.from_name('CO')
    res = HF(mol, 'sto3g').rhf()
    mopa = PopulationAnalysis(res)

    mohp = mopa.mohp(0, 1)
    moop = mopa.moop(0, 1)
    mobi = mopa.mobi(0, 1)

    print('MOHP, MOOP, and MOBI values of canonical Hartree-Fock orbitals')
    for i, (e, h, o, b) in enumerate(zip(res['orbe'], mohp, moop, mobi)):
        print('%3i %12.4f %12.4f %12.4f %12.4f' % (i+1, e, h, o, b))

Example output::

  MOHP, MOOP, and MOBI values of canonical Hartree-Fock orbitals
    1     -20.3914       0.0319      -0.0017      +0.0005
    2     -11.0902       0.0094      -0.0009      +0.0003
    3      -1.4047      -0.4347       0.2937      -0.0219
    4      -0.6899       0.1977      -0.0663      +0.3657
    5      -0.5094      -0.2736       0.1562      +0.5121
    6      -0.5094      -0.2736       0.1562      +0.5121
    7      -0.4409       0.0813      -0.0375      +0.2263
    8       0.2865       0.4489      -0.2562      -0.8401
    9       0.2865       0.4489      -0.2562      -0.8401
   10       0.9253       5.1204      -2.6744      -1.4625

MOHP, MOOP, and MOBI analysis of Foster-Boys localized orbitals
***************************************************************

It is often insightful to perform population analysis on *localized*
molecular orbitals. Since Foster-Boys localization corresponds to a unitary
transformation of the occupied orbital subspace, the output of a
Foster-Boys localization can be passed directly to the `PopulationAnalysis` class.

The example below compares MOHP values obtained from canonical and
Foster-Boys localized orbitals.

.. code-block:: python

    from pyqint import MoleculeBuilder, HF, FosterBoys, PopulationAnalysis
    import numpy as np

    mol = MoleculeBuilder.from_name('CO')
    res = HF(mol, 'sto3g').rhf()
    mopa = PopulationAnalysis(res)
    mohp = mopa.mohp(0, 1)
    moop = mopa.moop(0, 1)
    mobi = mopa.mobi(0, 1)

    res_fb = FosterBoys(res).run()
    mopa_fb = PopulationAnalysis(res_fb)
    mohp_fb = mopa_fb.mohp(0, 1)
    moop_fb = mopa_fb.moop(0, 1)
    mobi_fb = mopa_fb.mobi(0, 1)

    print('Sum of MOHP (canonical orbitals):',
          np.sum(mohp[:7]))
    print('Sum of MOHP (Foster-Boys orbitals):',
          np.sum(mohp_fb[:7]))

    print('Sum of MOOP (canonical orbitals):',
          np.sum(moop[:7]))
    print('Sum of MOOP (Foster-Boys orbitals):',
          np.sum(moop_fb[:7]))

    print('Sum of MOBI (canonical orbitals):',
          np.sum(mobi[:7]))
    print('Sum of MOBI (Foster-Boys orbitals):',
          np.sum(mobi_fb[:7]))

Example output::

    Sum of MOHP (canonical orbitals): -0.6617486641766972
    Sum of MOHP (Foster-Boys orbitals): -0.6617486641766972
    Sum of MOOP (canonical orbitals): 0.49960515685531653
    Sum of MOOP (Foster-Boys orbitals): 0.49960515685531626
    Sum of MOBI (canonical orbitals): 1.5951431033478687
    Sum of MOBI (Foster-Boys orbitals): 1.5951431033478685

The results demonstrate that while individual MOHP, MOOP, and MOBI contributions may
differ significantly between canonical and localized orbitals, the **sum over
the occupied subspace is invariant under unitary transformations**. This
invariance reflects the fact that the total bonding character, electron density,
and total energy of the system are preserved under orbital localization
procedures.