"""Atomic density guesses for molecular SCF calculations."""

import warnings
from functools import lru_cache
from math import prod, sqrt

import numpy as np

from ..basis.cgf import CGF
from ..basis.spherical_harmonics import spherical_harmonics
from ..pyqint_core import PyQInt


# Canonical orthogonalization drops overlap eigenvalues below this value.
LINDEP_THRESHOLD = 1e-7

# Angular momentum of each subshell in Aufbau order (1s, 2s, 2p, ...).
_AUFBAU_ORDER = (0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1)
# Neutral transition atoms with an s-to-d promotion in the built-in bases.
_SD_PROMOTIONS = {24: 1, 29: 1, 41: 1, 42: 1, 44: 1, 45: 1, 46: 2, 47: 1}


def sad_density(cgfs, nuclei):
    """Return the spin-summed superposition of neutral atomic densities.

    ``cgfs`` and ``nuclei`` use the same format as
    :meth:`PyQInt.build_integrals_openmp`. Each atomic density is placed in
    the block for that atom's basis functions. The density can have
    fractional occupations, and its trace with the overlap matrix is based
    on the neutral atoms. The molecular SCF sets the actual occupations.

    Return ``None`` if a basis function cannot be assigned to one nucleus,
    as with multicenter symmetry-adapted functions or off-nucleus functions.
    """
    centers = {tuple(position): i for i, (position, _) in enumerate(nuclei)}
    blocks = [[] for _ in nuclei]
    for i, cgf in enumerate(cgfs):
        if not cgf.gtos:
            return None
        center = tuple(cgf.gtos[0].p)
        if center not in centers or any(tuple(gto.p) != center for gto in cgf.gtos):
            return None
        blocks[centers[center]].append(i)

    density = np.zeros((len(cgfs), len(cgfs)))
    for (_, charge), indices in zip(nuclei, blocks):
        if indices:
            # Primitive parameters identify the basis independently of geometry.
            basis = tuple(tuple((g.c, g.alpha, g.l, g.m, g.n) for g in cgfs[i].gtos)
                          for i in indices)
            density[np.ix_(indices, indices)] = _atomic_density(charge, basis)
    return density


def _shell_electrons(charge):
    """Count the neutral atom's Aufbau electrons of each angular momentum."""
    remaining = charge
    counts = [0, 0, 0, 0]
    for angular in _AUFBAU_ORDER:
        filled = min(remaining, 2 * (2 * angular + 1))
        counts[angular] += filled
        remaining -= filled
    promoted = _SD_PROMOTIONS.get(charge, 0)
    counts[0] -= promoted
    counts[2] += promoted
    return counts


def _cartesian_norm(powers):
    """Angular part of the normalization for a Cartesian Gaussian."""
    return 1 / sqrt(prod(prod(range(1, 2 * power, 2)) for power in powers))


def _spherical_polynomial(angular, magnetic, radial_power):
    """Express r**(2*radial_power) times a solid harmonic in Cartesian form."""
    coefficients = {tuple(powers): value * _cartesian_norm(powers)
                    for value, powers in spherical_harmonics[angular][magnetic]}
    for _ in range(radial_power):
        expanded = {}
        for powers, value in coefficients.items():
            for axis in range(3):
                term = powers[:axis] + (powers[axis] + 2,) + powers[axis + 1:]
                expanded[term] = expanded.get(term, 0.0) + value
        coefficients = expanded
    return {powers: value / _cartesian_norm(powers) for powers, value in coefficients.items()}


def _spherical_basis(basis, overlap, shell_electrons):
    """Orthogonalize the radial functions of each angular momentum channel."""
    shells = {}
    for index, primitives in enumerate(basis):
        powers = primitives[0][2:]
        if any(primitive[2:] != powers for primitive in primitives):
            return None
        key = (sum(powers), tuple(primitive[:2] for primitive in primitives))
        components = shells.setdefault(key, {})
        if powers in components:
            return None
        components[powers] = index

    channels = {}
    for (degree, _), components in shells.items():
        if degree not in spherical_harmonics or len(components) != (degree + 1) * (degree + 2) // 2:
            return None
        # A Cartesian shell also supplies lower angular momenta, e.g. the
        # sixth d function has s symmetry and a radial factor r**2.
        for angular in range(degree, -1, -2):
            radial = np.zeros((len(basis), 2 * angular + 1))
            for i, magnetic in enumerate(range(-angular, angular + 1)):
                polynomial = _spherical_polynomial(angular, magnetic, (degree - angular) // 2)
                for powers, weight in polynomial.items():
                    radial[components[powers], i] = weight
            channels.setdefault(angular, []).append(radial)

    transforms = []
    blocks = []
    for angular, functions in channels.items():
        degeneracy = 2 * angular + 1
        radial = np.stack(functions, axis=1)
        metric = np.einsum('umi,uv,vni->mn', radial, overlap, radial, optimize=True) / degeneracy
        values, vectors = np.linalg.eigh(metric)
        keep = values > LINDEP_THRESHOLD
        orthogonalizer = vectors[:, keep] / np.sqrt(values[keep])
        transform = np.einsum('upi,pq->uqi', radial, orthogonalizer).reshape(len(basis), -1)
        electrons = shell_electrons[angular] if angular < len(shell_electrons) else 0
        if electrons > 2 * transform.shape[1]:
            return None
        transforms.append(transform)
        blocks.append((degeneracy, transform.shape[1], electrons))
    if sum(block[2] for block in blocks) < sum(shell_electrons):
        return None
    return np.column_stack(transforms), blocks


def _spherical_density(fock, blocks):
    density = np.zeros_like(fock)
    first = 0
    for degeneracy, size, electrons in blocks:
        if electrons:
            block = np.s_[first:first + size, first:first + size]
            radial_size = size // degeneracy
            angular_fock = fock[block].reshape(radial_size, degeneracy, radial_size, degeneracy)
            # Project onto a central field so roundoff cannot split a shell.
            radial_fock = np.einsum('piqi->pq', angular_fock) / degeneracy
            radial_density = _aufbau_density(radial_fock, electrons / degeneracy)
            density[block] = np.kron(radial_density, np.eye(degeneracy))
        first += size
    return density


def _aufbau_density(fock, electrons):
    """Fill energy-ordered orbitals, sharing electrons among degenerate ones."""
    energies, orbitals = np.linalg.eigh(fock)
    occupations = np.zeros(len(energies))
    starts = np.r_[0, np.flatnonzero(np.diff(energies) > 1e-7) + 1, len(energies)]
    for first, last in zip(starts[:-1], starts[1:]):
        width = last - first
        filled = min(electrons, 2 * width)
        occupations[first:last] = filled / width
        electrons -= filled
    return (orbitals * occupations) @ orbitals.T


@lru_cache(maxsize=128)
def _atomic_density(charge, basis):
    """Solve an isolated atom by restricted SCF with fractional occupations."""
    origin = np.zeros(3)
    cgfs = []
    for primitives in basis:
        cgf = CGF(origin)
        for primitive in primitives:
            cgf.add_gto(*primitive)
        cgfs.append(cgf)

    integrator = PyQInt(eri_engine="auto")
    overlap, kinetic, attraction, eri = integrator.build_integrals_openmp(
        cgfs, [(origin, charge)], packed=True)
    spherical = _spherical_basis(basis, overlap, _shell_electrons(charge))
    if spherical is not None:
        transform, blocks = spherical
    else:
        # Custom atom-centered bases can lack members of a Cartesian shell.
        eigenvalues, eigenvectors = np.linalg.eigh(overlap)
        keep = eigenvalues > LINDEP_THRESHOLD
        if not np.any(keep):
            return np.zeros_like(overlap)
        transform = eigenvectors[:, keep] / np.sqrt(eigenvalues[keep])
        blocks = None

    def density_from_fock(fock):
        fock = transform.T @ fock @ transform
        return _aufbau_density(fock, charge) if blocks is None else _spherical_density(fock, blocks)

    hcore = kinetic + attraction
    density = density_from_fock(hcore)

    for _ in range(128):
        ao_density = transform @ density @ transform.T
        coulomb, exchange = integrator.contract_packed_tei(eri, ao_density)
        updated = density_from_fock(hcore + coulomb - 0.5 * exchange)
        if np.linalg.norm(updated - density) < 1e-8:
            return transform @ updated @ transform.T
        density = 0.5 * (density + updated)

    warnings.warn(f"SAD atomic SCF for Z={charge} did not converge; using its last density.",
                  RuntimeWarning, stacklevel=3)
    return transform @ density @ transform.T
