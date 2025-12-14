# -*- coding: utf-8 -*-

from pyqint import MoleculeBuilder, HF, PyQInt
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def main():
    mol = MoleculeBuilder().from_name('ethylene')
    cgfs, nuclei = mol.build_basis('sto3g')
    from pyqint import CGF, GTO
    res = HF().rhf(mol, 'sto3g')
    
    # build transformation matrix that casts the original basis onto a
    # symmetry adapted basis; no normalization is applied
    B = np.zeros((14,14))
    for i in range(0,5):
        B[i*2,i] = 1
        B[i*2,i+7] = 1
        B[i*2+1,i] = 1
        B[i*2+1,i+7] = -1
    
    B[10,5] = 1
    B[10,6] = 1
    B[10,-2] = 1
    B[10,-1] = 1
    
    B[11,5] = 1
    B[11,6] = -1
    B[11,-2] = 1
    B[11,-1] = -1
    
    B[12,5] = 1
    B[12,6] = 1
    B[12,-2] = -1
    B[12,-1] = -1
    
    B[13,5] = 1
    B[13,6] = -1
    B[13,-2] = -1
    B[13,-1] = 1

    # build yet another transformation basis that re-orders the functions such
    # that all irreps belonging to the same symmetry group are consecutive; this
    # yields block-diagonal matrices
    order = []
    overlap = B @ res['overlap'] @ B.T
    for i in range(len(overlap)):
        for j in range(len(overlap)):
            if abs(overlap[i,j]) > 0.1:
                if j not in order:
                    order.append(j)
    n = len(order)
    P = np.zeros((n, n), dtype=int)
    P[np.arange(n), order] = 1
    
    symfuncs = {
        'A$_{g}$'  : 4,
        'B$_{3u}$' : 4,
        'B$_{2g}$' : 1,
        'B$_{1u}$' : 1,
        'B$_{1g}$' : 2,
        'B$_{2u}$' : 2
    }
    symlabels= []
    for k,v in symfuncs.items():
        for i in range(v):
            symlabels.append(k + '(%i)' % (i+1))
    overlap = P @ overlap @ P.T

    basislabels= []
    counts = {}
    for n in nuclei:
        if n[1] not in counts.keys():
            counts[n[1]] = 1
        else:
            counts[n[1]] += 1
        
        if n[1] == 6:
            for s in ['1s', '2s', '2p_{x}', '2p_{y}', '2p_{z}']:
                basislabels.append('C$_{%s}^{(%i)}$' % (s, counts[n[1]]))
        else:
            basislabels.append('H$_{1s}^{(%i)}$' % (counts[n[1]]))

    overlap = P @ overlap @ P.T

    fig, ax = plt.subplots(1,1, dpi=144, figsize=(7,7))
    plot_matrix(ax, P @ B, basislabels, symlabels)
    plt.show()

    # construct new basis
    B = P @ B
    cgfs_symad = [CGF() for i in range(len(B))]
    for i in range(len(B)): # loop over new basis functions
        for j in range(len(cgfs)): # loop over old basis functions
            if abs(B[i,j]) > 0.01:  # verify non-negligble contribution
                for g in cgfs[j].gtos:
                    cgfs_symad[i].gtos.append(GTO(g.c*B[i,j], g.p, g.alpha, g.l, g.m, g.n))

    # re-perform Hartree-Fock calculation using the symmetry adapted basis
    res = HF().rhf(mol, cgfs_symad, verbose=True)
    fig, ax = plt.subplots(1,1, dpi=144, figsize=(7,7))
    plot_matrix(ax, res['overlap'], symlabels, symlabels)
    plt.show()

    fig, ax = plt.subplots(1,1, dpi=144, figsize=(7,7))
    plot_matrix(ax, res['fock'], symlabels, symlabels)
    plt.show()

def plot_matrix(ax, mat, xlabels, ylabels, xlabelrot = 0):
    """
    Produce plot of matrix
    """
    mv = np.max(np.abs(mat))
    ax.imshow(mat, vmin=-mv, vmax=mv, cmap='PiYG')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(i, j, '%.2f' % mat[j,i], ha='center', va='center',
                    fontsize=7, color='black' if abs(mat[j,i]) < (0.7 * mv) else 'white')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.hlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    ax.vlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    
    # add basis functions as axes labels
    ax.set_xticks(np.arange(0, mat.shape[0]))
    ax.set_xticklabels(xlabels, rotation=xlabelrot)
    ax.set_yticks(np.arange(0, mat.shape[0]))
    ax.set_yticklabels(ylabels, rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=7)
    
if __name__ == '__main__':
    main()