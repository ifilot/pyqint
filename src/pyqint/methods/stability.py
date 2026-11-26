"""Internal stability of real restricted Hartree–Fock orbitals."""

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import LinearOperator, eigsh


def lowest_rhf_mode(orbitals, fock, nocc, coulomb_exchange):
    """Return the lowest eigenpair of A + B in occupied–virtual coordinates."""
    occupied = orbitals[:, :nocc]
    virtual = orbitals[:, nocc:]
    shape = (virtual.shape[1], nocc)
    size = shape[0] * shape[1]
    if size == 0:
        return None, None

    f_oo = occupied.T @ fock @ occupied
    f_vv = virtual.T @ fock @ virtual

    def hessian_action(vector):
        rotation = vector.reshape(shape)
        # A rotation changes the density by 2(Cv κ Coᵀ + Co κᵀ Cvᵀ).
        # The resulting Fock response gives the two-electron part of A + B.
        transition = virtual @ rotation @ occupied.T
        coulomb, exchange = coulomb_exchange(2 * (transition + transition.T))
        response = virtual.T @ (coulomb - 0.5 * exchange) @ occupied
        return (f_vv @ rotation - rotation @ f_oo + response).ravel()

    if size <= 64:
        hessian = np.column_stack([hessian_action(unit) for unit in np.eye(size)])
        values, vectors = np.linalg.eigh((hessian + hessian.T) * 0.5)
        return float(values[0]), vectors[:, 0].reshape(shape)

    gaps = (np.diag(f_vv)[:, None] - np.diag(f_oo)).ravel()
    value, vector = _davidson_lowest(hessian_action, gaps)
    return value, vector.reshape(shape)


def _davidson_lowest(action, gaps):
    """Find the lowest mode using orbital-energy gaps as a preconditioner."""
    # Sample all symmetry blocks, with more weight on low-energy rotations.
    start = np.random.default_rng(0).standard_normal(gaps.size)
    start /= gaps - gaps.min() + 0.2
    trial = start / np.linalg.norm(start)
    basis = np.empty((gaps.size, 0))
    images = np.empty_like(basis)

    for _ in range(100):
        basis = np.column_stack((basis, trial))
        images = np.column_stack((images, action(trial)))
        projected = basis.T @ images
        values, vectors = np.linalg.eigh((projected + projected.T) * 0.5)
        value = values[0]
        vector = basis @ vectors[:, 0]
        residual = images @ vectors[:, 0] - value * vector
        if np.linalg.norm(residual) < 1e-7:
            return float(value), vector

        if basis.shape[1] == 24:
            basis = basis @ vectors[:, :3]
            images = images @ vectors[:, :3]

        trial = residual / np.maximum(np.abs(gaps - value), 0.1)
        for _ in range(2):
            trial -= basis @ (basis.T @ trial)
        norm = np.linalg.norm(trial)
        if norm < 1e-12:
            break
        trial /= norm

    # Difficult spectra can defeat the gap preconditioner. Retry with a
    # matrix-free symmetric Lanczos solver before reporting a failure.
    hessian = LinearOperator((gaps.size, gaps.size), matvec=action, dtype=np.float64)
    values, vectors = eigsh(hessian, k=1, which="SA", v0=start, tol=1e-7)
    return float(values[0]), vectors[:, 0]


def rotate_rhf_occupied(orbitals, rotation):
    """Apply an orthogonal occupied–virtual rotation and return occupied MOs."""
    nvirt, nocc = rotation.shape
    generator = np.zeros((nocc + nvirt, nocc + nvirt))
    generator[nocc:, :nocc] = rotation
    generator[:nocc, nocc:] = -rotation.T
    return orbitals @ expm(generator)[:, :nocc]
