import unittest
from pyqint import PyQInt, cgf, gto, Molecule, HF
from copy import deepcopy
import numpy as np
import multiprocessing
import os
from nose.tools import nottest
from scipy.stats import linregress

def main():
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
    integrator = PyQInt()
    
    P = res['density']
    
    S = res['overlap']
    T = res['kinetic']
    V = res['nuclear']
    teint = res['teint']
    C = res['orbc']
    orbe = res['orbe']
    
    # calculate G
    G = np.zeros(S.shape)
    for i in range(S.shape[0]):
        for j in range(S.shape[0]):
            for k in range(S.shape[0]):
                for l in range(S.shape[0]):
                    idx_rep = integrator.teindex(i,j,l,k)
                    idx_exc = integrator.teindex(i,k,l,j)
                    G[i,j] += P[k,l] * (teint[idx_rep] - 0.5 * teint[idx_exc])

    # build Fock matrix
    F = T + V + G
    
    print(np.einsum('ij,ji->i', T+V, C))
    
if __name__ == '__main__':
    main()