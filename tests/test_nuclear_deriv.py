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

        np.testing.assert_almost_equal(forces, forces_fd, 4)

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
                        np.testing.assert_almost_equal(force, vals[i,j,k,l], 4)

def calculate_force_finite_difference(cgf_id1, cgf_id2, nuc_id, coord, probe_atom):
    # build integrator object
    integrator = PyQInt()

    # distance
    diff = 0.00001

    vals = np.zeros(2)
    for i,v in enumerate([-1,1]):
        # build hydrogen molecule
        mol = Molecule('H2O')
        mol.add_atom('O', 0.00000, -0.07579, 0.00000, unit='angstrom')
        mol.add_atom('H', 0.86681, 0.60144, 0.00000, unit='angstrom')
        mol.add_atom('H',  -0.86681, 0.60144, 0.00000, unit='angstrom')

        # adjust molecule
        mol.get_atoms()[nuc_id][1][coord] += v * diff / 2

        # build basis
        cgfs, nuclei = mol.build_basis('sto3g')

        # calculate values
        O = nuclei[probe_atom][0]
        Ochg = nuclei[probe_atom][1]
        vals[i] = integrator.nuclear(cgfs[cgf_id1], cgfs[cgf_id2], O, Ochg)

    return (vals[1] - vals[0]) / diff

def energy_nuclear_repulsion(nuclei):
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
    forces = np.zeros((3,3))

    sz = 0.0001

    for i in range(0, len(mol.get_atoms())): # loop over nuclei
        for j in range(0, 3): # loop over directions
            mol1 = deepcopy(mol)
            mol1.get_atoms()[i][1][j] -= sz / 2
            mol2 = deepcopy(mol)
            mol2.get_atoms()[i][1][j] += sz / 2

            cgfs, nuclei1= mol1.build_basis('sto3g')
            cgfs, nuclei2 = mol2.build_basis('sto3g')

            energy1 = energy_nuclear_repulsion(nuclei1)
            energy2 = energy_nuclear_repulsion(nuclei2)

            forces[i,j] = (energy2 - energy1) / sz

    return forces

if __name__ == '__main__':
    unittest.main()
