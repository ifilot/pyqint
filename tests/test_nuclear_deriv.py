import unittest
from copy import deepcopy
import numpy as np
from pyqint import MoleculeBuilder


class TestNuclearDeriv(unittest.TestCase):
    """
    Unit tests for analytic nuclear-nuclear repulsion gradients.

    The analytic gradients are compared against finite-difference
    reference values obtained from the nuclear repulsion energy E_nn.
    """

    def test_nuclear_derivatives_h2(self):
        """
        Test analytic nuclear-nuclear repulsion gradients for H2.

        This is a minimal sanity check: two nuclei with a single
        internuclear distance. The gradient must be equal and opposite
        on the two atoms.
        """
        mol = MoleculeBuilder.from_name('h2')
        _, nuclei = mol.build_basis('sto3g')

        # Finite-difference reference gradient
        forces_fd = calculate_deriv_nuclear_repulsion_finite_difference(mol)

        # Analytic gradient
        forces = np.zeros(forces_fd.shape)
        for i in range(len(mol.get_atoms())):      # loop over nuclei
            for j in range(3):                     # loop over x, y, z
                forces[i, j] = deriv_nuclear_repulsion(nuclei, i, j)

        np.testing.assert_allclose(forces, forces_fd, rtol=1e-5, atol=1e-8)

    def test_nuclear_derivatives_h2o(self):
        """
        Test analytic nuclear-nuclear repulsion gradients for H2O
        using the STO-3G basis set.

        The analytic gradients are compared against a 5-point
        finite-difference reference. Agreement within 1e-5 relative
        tolerance confirms correctness of the implementation.
        """
        mol = MoleculeBuilder.from_name('h2o')
        _, nuclei = mol.build_basis('sto3g')

        # Finite-difference reference gradient
        forces_fd = calculate_deriv_nuclear_repulsion_finite_difference(mol)

        # Analytic gradient
        forces = np.zeros(forces_fd.shape)
        for i in range(len(mol.get_atoms())):      # loop over nuclei
            for j in range(3):                     # loop over x, y, z
                forces[i, j] = deriv_nuclear_repulsion(nuclei, i, j)

        np.testing.assert_allclose(forces, forces_fd, rtol=1e-5, atol=1e-8)


def energy_nuclear_repulsion(nuclei):
    """
    Compute the nuclear-nuclear repulsion energy:

        E_nn = sum_{i<j} Z_i Z_j / |R_i - R_j|

    Parameters
    ----------
    nuclei : list of (position, charge)
        Nuclear positions and charges.

    Returns
    -------
    float
        Nuclear repulsion energy.
    """
    energy = 0.0
    for i in range(len(nuclei)):
        ri = np.array(nuclei[i][0], dtype=float)
        Zi = float(nuclei[i][1])
        for j in range(i + 1, len(nuclei)):
            rj = np.array(nuclei[j][0], dtype=float)
            Zj = float(nuclei[j][1])
            Rij = np.linalg.norm(ri - rj)
            energy += Zi * Zj / Rij
    return energy


def deriv_nuclear_repulsion(nuclei, nucid, coord):
    """
    Analytic derivative of the nuclear-nuclear repulsion energy E_nn
    with respect to a single Cartesian coordinate of a nucleus.

    This returns the **energy gradient** component:

        dE_nn / dR_{nucid, coord}

    Note
    ----
    This is the gradient, not the force.
    Forces are given by: F = -∇E.

    Parameters
    ----------
    nuclei : list of (position, charge)
        Nuclear positions and charges.
    nucid : int
        Index of the nucleus with respect to which the derivative is taken.
    coord : int
        Cartesian direction (0 = x, 1 = y, 2 = z).

    Returns
    -------
    float
        Gradient component dE_nn / dR.
    """
    pc = np.array(nuclei[nucid][0], dtype=float)
    Zc = float(nuclei[nucid][1])

    dE = 0.0
    for i in range(len(nuclei)):
        if i == nucid:
            continue
        pi = np.array(nuclei[i][0], dtype=float)
        Zi = float(nuclei[i][1])

        r = pi - pc
        rnorm = np.linalg.norm(r)
        dE += Zc * Zi * (pc[coord] - pi[coord]) / (rnorm ** 3)

    # Return gradient (not force)
    return -dE


def calculate_deriv_nuclear_repulsion_finite_difference(mol):
    """
    Finite-difference reference for the nuclear-nuclear repulsion gradient.

    A 5-point central stencil is used to approximate the derivative
    of E_nn with respect to each nuclear coordinate.

    Parameters
    ----------
    mol : Molecule
        Molecule object.

    Returns
    -------
    ndarray, shape (N_atoms, 3)
        Finite-difference nuclear gradients.
    """
    forces = np.zeros((len(mol.get_atoms()), 3), dtype=float)

    h = 1e-2  # finite-difference step size
    for i in range(len(mol.get_atoms())):
        for coord in range(3):

            mol_m2 = deepcopy(mol)
            mol_m2.get_atoms()[i][1][coord] -= 2 * h
            _, nuclei_m2 = mol_m2.build_basis('sto3g')
            e_m2 = energy_nuclear_repulsion(nuclei_m2)

            mol_m1 = deepcopy(mol)
            mol_m1.get_atoms()[i][1][coord] -= h
            _, nuclei_m1 = mol_m1.build_basis('sto3g')
            e_m1 = energy_nuclear_repulsion(nuclei_m1)

            mol_p1 = deepcopy(mol)
            mol_p1.get_atoms()[i][1][coord] += h
            _, nuclei_p1 = mol_p1.build_basis('sto3g')
            e_p1 = energy_nuclear_repulsion(nuclei_p1)

            mol_p2 = deepcopy(mol)
            mol_p2.get_atoms()[i][1][coord] += 2 * h
            _, nuclei_p2 = mol_p2.build_basis('sto3g')
            e_p2 = energy_nuclear_repulsion(nuclei_p2)

            # 5-point central difference formula
            forces[i, coord] = (
                -e_p2 + 8.0 * e_p1 - 8.0 * e_m1 + e_m2
            ) / (12.0 * h)

    return forces


if __name__ == "__main__":
    unittest.main()