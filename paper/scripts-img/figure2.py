# -*- coding: utf-8 -*-

from pyqint import PyQInt, MoleculeBuilder, HF
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#
# Figure 2 of the JOSE paper: PyQInt: A Teaching-Oriented Hartree–Fock 
# Implementation in Python
#

def main():
    # calculate sto3g coefficients for CO
    cgfs, coeff, orbe = calculate_co()

    # visualize orbitals
    fig, ax = plt.subplots(2,5, figsize=(12, 5), dpi=300)
    sz = 3
    for i in range(0,2):
        for j in range(0,5):
            dens = plot_wavefunction(cgfs, coeff[:,i*5+j], sz=sz)
            limit = max(abs(np.min(dens)), abs(np.max(dens)) )
            im = ax[i,j].contourf(dens, origin='lower',
              extent=[-sz, sz, -sz, sz], cmap='PiYG', vmin=-limit, vmax=limit,
              levels=15)
            im = ax[i,j].contour(dens, origin='lower', colors='black',
              extent=[-sz, sz, -sz, sz], vmin=-limit, vmax=limit,
              levels=15)
            ax[i,j].set_xlabel('x [Bohr]')
            ax[i,j].set_ylabel('z [Bohr]')
            ax[i,j].set_aspect('equal', adjustable='box')
            ax[i,j].set_xticks(np.linspace(-3,3, 7))
            ax[i,j].set_yticks(np.linspace(-3,3, 7))
            ax[i,j].grid(linestyle='--', alpha=0.5)
            ax[i,j].set_title('$\psi_{%i}$ ($\epsilon_{%i}$ = %6.4f Ht)' %
                              (i*5+j+1, i*5+j+1, orbe[i*5+j]))
    plt.tight_layout()
    plt.savefig('../img/figure2-orbitals-co-contour.pdf')

def calculate_co():
    mol = MoleculeBuilder().from_name('CO')

    result = HF().rhf(mol, 'sto3g')

    return result['cgfs'], result['orbc'], result['orbe']

def plot_wavefunction(cgfs, coeff, sz=3.5):
    # build integrator
    integrator = PyQInt()

    # build grid
    x = np.linspace(-sz, sz, 150)
    z = np.linspace(-sz, sz, 150)
    xx, zz = np.meshgrid(x,z)
    yy = np.zeros(len(x) * len(z))
    grid = np.vstack([xx.flatten(), yy, zz.flatten()]).reshape(3,-1).T
    res = integrator.plot_wavefunction(grid, coeff, cgfs).reshape((len(z), len(x)))

    return res

if __name__ == '__main__':
    main()