import unittest
from pyqint import PyQInt, cgf, gto, Molecule, HF
from copy import deepcopy
import numpy as np
import multiprocessing
import os
from nose.tools import nottest
from scipy.stats import linregress

class TestHFDeriv(unittest.TestCase):

    def testHartreeFockForcesCore(self):
        """
        Test Hartree-Fock calculation on water using STO-3G basis set
        """
        mol = Molecule()
        mol.add_atom('O', 0.0, 0.0, 0.0)
        mol.add_atom('H', 0.7570, 0.5860, 0.0)
        mol.add_atom('H', -0.7570, 0.5860, 0.0)
        basis = 'sto3g'

        # calculate forces using analytical derivatives
        solver = HF()
        res = solver.rhf(mol, basis, calc_forces=False)

        # get density matrix
        C = res['orbc']
        P = res['density']
        e = res['orbe']

        # calculate forces via routines
        forces = np.zeros((len(mol.atoms), 3))
        for i in range(0, len(mol.atoms)): # loop over nuclei
            for j in range(0, 3): # loop over directions
                forces[i,j] = rhf_force_nuc_dir(mol, basis, C, P, e, i, j)

        # calculate forces using finite difference
        forces_fd = calculate_forces_finite_difference(mol)

        np.testing.assert_almost_equal(forces, forces_fd)

    # @nottest
    # def testHartreeFockForces(self):
    #     """
    #     Test Hartree-Fock calculation on water using STO-3G basis set
    #     """
    #     mol = Molecule()
    #     mol.add_atom('O', 0.0, 0.0, 0.0)
    #     mol.add_atom('H', 0.7570, 0.5860, 0.0)
    #     mol.add_atom('H', -0.7570, 0.5860, 0.0)

    #     # calculate forces using analytical derivatives
    #     solver = HF()
    #     res = solver.rhf(mol, 'sto3g', calc_forces=True)

    #     # calculate forces using finite difference
    #     forces = calculate_forces_finite_difference(mol)

    #     np.testing.assert_almost_equal(res['forces'], forces)

def perform_hf(mol):
    sol = HF().rhf(mol, 'sto3g')
    return sol

def calculate_forces_finite_difference(mol):
    forces = np.zeros((3,3))

    d = 0.1
    npoints = 5

    for i in range(0, len(mol.atoms)): # loop over nuclei
        for j in range(0, 3): # loop over directions

            # build five point stencil
            energies = np.zeros(npoints)
            x = np.linspace(-d, d, 5)
            for k,z in enumerate(x):
                molnew = deepcopy(mol)
                molnew.atoms[i][1][j] += z
                energies[k] = perform_hf(molnew)['energy']
                
            res = linregress(x, energies)

            forces[i,j] = res[0]
            print(i,j,"R2",res[2]**2)

    return forces

def calculate_forces_core_finite_difference(mol):
    forces = np.zeros((3,3))

    d = 0.05
    npoints = 5

    print('Assess fit quality')
    for i in range(0, len(mol.atoms)): # loop over nuclei
        for j in range(0, 3): # loop over directions
            
            # build five point stencil
            energies = np.zeros(npoints)
            x = np.linspace(-d, d, 5)
            for k,z in enumerate(x):
                molnew = deepcopy(mol)
                molnew.atoms[i][1][j] += z
                energies[k] = perform_hf(molnew)['ecore']
                
            res = linregress(x, energies)

            forces[i,j] = res[0]
            print(i,j,"R2",res[2]**2)

    return forces

def rhf_force_nuc_dir(mol, basis, C, P, e, nucleus, direction):
    # build cgfs, nuclei and calculate nr of electrons
    cgfs, nuclei = mol.build_basis(basis)
    nelec = int(np.sum([at[1] for at in nuclei]))

    # build integrator object
    integrator = PyQInt()

    # build overlap and kinetic derivatives
    S = np.zeros((len(cgfs), len(cgfs)))
    for i in range(0, len(cgfs)):
        for j in range(0, len(cgfs)):
            S[i,j] = integrator.overlap_deriv(cgfs[i], cgfs[j], nuclei[nucleus][0], direction)

    # build kinetic derivatives
    T = np.zeros((len(cgfs), len(cgfs)))
    for i in range(0, len(cgfs)):
        for j in range(0, len(cgfs)):
            T[i,j] = integrator.kinetic_deriv(cgfs[i], cgfs[j], nuclei[nucleus][0], direction)

    # build nuclear derivatives
    V = np.zeros((len(cgfs), len(cgfs)))
    for i in range(0, len(cgfs)):
        for j in range(0, len(cgfs)):
            for k in range(0, len(nuclei)):
                V[i,j] += integrator.nuclear_deriv(cgfs[i], cgfs[j], nuclei[k][0], nuclei[k][1], nuclei[nucleus][0], direction)

    # build Q matrix
    Q = np.zeros(S.shape)
    for i in range(S.shape[0]):
        for j in range(S.shape[0]):
            for k in range(0,int(nelec/2)):
                Q[j,i] += 2.0 * e[k] * C[i,k] * C[j,k]

    # build two-electron derivatives
    N = len(cgfs)
    teint_calc = np.multiply(np.ones(integrator.teindex(N,N,N,N)), -1.0)
    teint = np.zeros(integrator.teindex(N,N,N,N))
    for i, cgf1 in enumerate(cgfs):
        for j, cgf2 in enumerate(cgfs):
            ij = i*(i+1)/2 + j
            for k, cgf3 in enumerate(cgfs):
                for l, cgf4 in enumerate(cgfs):
                    kl = k * (k+1)/2 + l
                    if ij <= kl:
                        idx = integrator.teindex(i,j,k,l)
                        if teint_calc[idx] < 0:
                            teint_calc[idx] = 1
                            teint[idx] = integrator.repulsion_deriv(cgfs[i], cgfs[j], cgfs[k], cgfs[l], nuclei[nucleus][0], direction)

    # build H-core derivatives
    Hcore = T + V

    # calculate electronic derivate
    deriv = 0.0
    for i in range(0, len(cgfs)):
        for j in range(0, len(cgfs)):
            deriv += P[j,i] * Hcore[i,j]
            for k in range(0, len(cgfs)):
                for l in range(0, len(cgfs)):
                    idx_rep = integrator.teindex(i,j,k,l)
                    idx_exc = integrator.teindex(i,l,k,j)
                    deriv += 0.5 * P[j,i] * P[k,l] * (teint[idx_rep] - 0.5 * teint[idx_exc])
            deriv -= Q[j,i] * S[i,j]

    # calculate nuclear derivative
    Vnn = 0.0
    pc = nuclei[nucleus][0]
    for i in range(0, len(nuclei)):
        if nucleus != i:
            pi = nuclei[i][0]
            Vnn += nuclei[nucleus][1] * nuclei[i][1] * (pi[direction] - pc[direction]) / np.linalg.norm(pi - pc)**3

    return deriv + Vnn

if __name__ == '__main__':
    unittest.main()
