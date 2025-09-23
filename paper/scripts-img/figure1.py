# -*- coding: utf-8 -*-

from pyqint import PyQInt, MoleculeBuilder, HF
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

#
# Figure 2 of the JOSE paper: PyQInt: A Teaching-Oriented Hartree–Fock 
# Implementation in Python
#

def main():
    # calculate sto3g coefficients for CO
    cgfs, coeff, orbe = calculate_co()
    fig, ax = plt.subplots(1,1,dpi=300, figsize=(5,5))
    plot_matrix(ax, coeff, ['O-1s', 'O-2s', 'O-2p$_{x}$', 'O-2p$_{y}$', 'O-2p$_{z}$',
                            'C-1s', 'C-2s', 'C-2p$_{x}$', 'C-2p$_{y}$', 'C-2p$_{z}$'])
    ax.set_xticks(np.arange(0,len(orbe)), 
                  ["$\psi_{%i}$ (%6.4f Ht)" % (i+1,e) for i,e in enumerate(orbe)], 
                  rotation=-90)
    ax.set_xlabel('Molecular orbital')
    ax.set_ylabel('Atomic orbital')
    ax.set_title('Coefficient matrix $\mathbf{C}$')
    plt.tight_layout()
    plt.savefig('../img/figure1-co-coefficient-matrix.pdf')

def calculate_co():
    mol = MoleculeBuilder().from_name('CO')
    result = HF().rhf(mol, 'sto3g')
    return result['cgfs'], result['orbc'], result['orbe']

def plot_matrix(ax, mat, labels, title = None, xlabelrot = 0, **kwargs):
    """
    Produce plot of matrix
    """
    ax.imshow(mat, vmin=-1, vmax=1, cmap='PiYG')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(i, j, '%.2f' % mat[j,i], ha='center', va='center',
                    fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.hlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    ax.vlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    
    # add basis functions as axes labels
    ax.set_xticks(np.arange(0, mat.shape[0]))
    ax.set_xticklabels(labels, rotation=xlabelrot)
    ax.set_yticks(np.arange(0, mat.shape[0]))
    ax.set_yticklabels(labels, rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=7)
    
    # add title if supplied
    if title:
        if 'titlefontsize' in kwargs:
            titlefontsize = kwargs['titlefontsize']
            ax.set_title(title, fontsize=titlefontsize)
        else:
            ax.set_title(title)

if __name__ == '__main__':
    main()