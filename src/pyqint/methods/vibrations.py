"""Harmonic analysis of a Cartesian energy Hessian."""

import numpy as np
from scipy import constants

from ..basis import CGF


WAVENUMBER_FACTOR = np.sqrt(
    constants.value("Hartree energy") / (constants.atomic_mass * constants.value("Bohr radius") ** 2)
) / (2 * np.pi * constants.c * 100)


def displace_basis(cgfs, coordinates, displacement):
    """Copy custom CGFs with each primitive following its nucleus."""
    def move(position, *, required):
        position = np.asarray(position, dtype=float)
        # Use the same squared-distance criterion as the nuclear derivatives.
        matches = np.flatnonzero(np.sum((coordinates - position) ** 2, axis=1) < 1e-4)
        if matches.size == 1:
            return position + displacement[matches[0]]
        if required or matches.size:
            raise ValueError("Vibrational analysis requires each primitive to follow a unique nucleus.")
        return position.copy()

    basis = []
    for cgf in cgfs:
        displaced = CGF(move(cgf.p, required=False))
        for gto in cgf.gtos:
            displaced.add_gto_with_position(gto.c, move(gto.p, required=True), gto.alpha, gto.l, gto.m, gto.n)
        basis.append(displaced)
    return basis


def harmonic_analysis(hessian, coordinates, masses):
    """Project rigid motion and return signed wavenumbers and Cartesian modes."""
    centered = coordinates - np.average(coordinates, axis=0, weights=masses)
    sqrt_mass = np.sqrt(masses)
    inverse_mass = np.repeat(1 / sqrt_mass, 3)
    weighted = inverse_mass[:, None] * hessian * inverse_mass[None, :]

    # The span of the infinitesimal rigid motions has rank 3, 5, or 6.
    axes = np.eye(3)
    translations = np.broadcast_to(axes, (len(masses), 3, 3))
    rotations = np.cross(axes[:, None, :], centered[None, :, :]).transpose(1, 2, 0)
    radius = np.max(np.linalg.norm(centered, axis=1))
    if radius:
        rotations /= radius
    rigid = np.concatenate((translations, rotations), axis=2) * sqrt_mass[:, None, None]
    vectors, singular_values, _ = np.linalg.svd(rigid.reshape(-1, 6), full_matrices=True)
    rank = np.count_nonzero(singular_values > singular_values[0] * 1e-10)
    vibrations = vectors[:, rank:]

    eigenvalues, modes = np.linalg.eigh(vibrations.T @ weighted @ vibrations)
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * WAVENUMBER_FACTOR
    cartesian = inverse_mass[:, None] * (vibrations @ modes)
    cartesian /= np.linalg.norm(cartesian, axis=0)
    return frequencies, cartesian.T.reshape(len(frequencies), len(masses), 3)
