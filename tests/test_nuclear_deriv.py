import unittest
from pyqint import PyQInt, Molecule
from copy import deepcopy
import numpy as np
import os

class TestNuclearDeriv(unittest.TestCase):

    def testNuclearRepulsionDerivatives(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)

        cgfs, nuclei = mol.build_basis('sto3g')

        forces_fd = calculate_deriv_nuclear_repulsion_finite_difference(mol)

        forces = np.zeros(forces_fd.shape)
        for i in range(0, len(mol.get_atoms())): # loop over nuclei
            for j in range(0, 3): # loop over directions
                forces[i,j] = deriv_nuclear_repulsion(nuclei, i, j)

        np.testing.assert_almost_equal(forces, forces_fd, 5)

    def test_derivatives_h2o_fulltest(self):
        """
        Test Derivatives of water
        """
        # build integrator object
        integrator = PyQInt()

        # build hydrogen molecule
        mol = Molecule('H2O')
        mol.add_atom('O', 0.00000, -0.07579, 0.00000, unit='angstrom')
        mol.add_atom('H', 0.86681, 0.60144, 0.00000, unit='angstrom')
        mol.add_atom('H',  -0.86681, 0.60144, 0.00000, unit='angstrom')
        cgfs, nuclei = mol.build_basis('sto3g')

        # get position and charge
        probe_atom = 0
        O = nuclei[probe_atom][0]
        Ochg = nuclei[probe_atom][1]

        # load results from file
        fname = os.path.join(os.path.dirname(__file__), 'results', 'nuclear_deriv_h2o.txt')
        vals = np.loadtxt(fname).reshape((len(cgfs), len(cgfs), 3, 3))
        for i in range(0, len(cgfs)): # loop over cgfs
            for j in range(0, len(cgfs)): # loop over cgfs
                for k in range(0,3):  # loop over nuclei
                    for l in range(0,3):  # loop over directions
                        force = integrator.nuclear_deriv(cgfs[i], cgfs[j], O, Ochg, nuclei[k][0], l)
                        np.testing.assert_almost_equal(force, vals[i,j,k,l], 5)

def energy_nuclear_attraction(nuclei):
    energy = 0.0
    for i in range(0, len(nuclei)):
        for j in range(i+1, len(nuclei)):
            r = np.linalg.norm(np.array(nuclei[i][0]) - np.array(nuclei[j][0]))
            energy += nuclei[i][1] * nuclei[j][1] / r
    return energy

def deriv_nuclear_repulsion(nuclei, nucid, coord):
    Vnn = 0.0
    pc = nuclei[nucid][0]
    for i in range(0, len(nuclei)):
        if nucid != i:
            pi = nuclei[i][0]
            Vnn += nuclei[nucid][1] * nuclei[i][1] * (pi[coord] - pc[coord]) / np.linalg.norm(pi - pc)**3

    return Vnn

def calculate_deriv_nuclear_repulsion_finite_difference(mol):
    forces = np.zeros((len(mol.get_atoms()), 3))

    h = 1e-2

    for i in range(len(mol.get_atoms())):      # loop over nuclei
        for j in range(3):                     # loop over x, y, z

            # x - 2h
            mol_m2 = deepcopy(mol)
            mol_m2.get_atoms()[i][1][j] -= 2 * h
            _, nuclei_m2 = mol_m2.build_basis('sto3g')
            e_m2 = energy_nuclear_attraction(nuclei_m2)

            # x - h
            mol_m1 = deepcopy(mol)
            mol_m1.get_atoms()[i][1][j] -= h
            _, nuclei_m1 = mol_m1.build_basis('sto3g')
            e_m1 = energy_nuclear_attraction(nuclei_m1)

            # x + h
            mol_p1 = deepcopy(mol)
            mol_p1.get_atoms()[i][1][j] += h
            _, nuclei_p1 = mol_p1.build_basis('sto3g')
            e_p1 = energy_nuclear_attraction(nuclei_p1)

            # x + 2h
            mol_p2 = deepcopy(mol)
            mol_p2.get_atoms()[i][1][j] += 2 * h
            _, nuclei_p2 = mol_p2.build_basis('sto3g')
            e_p2 = energy_nuclear_attraction(nuclei_p2)

            # 5-point stencil
            forces[i, j] = (
                -e_p2
                + 8.0 * e_p1
                - 8.0 * e_m1
                + e_m2
            ) / (12.0 * h)

    return forces

if __name__ == '__main__':
    unittest.main()
