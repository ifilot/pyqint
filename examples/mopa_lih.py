from pyqint import Molecule, HF
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

#
# Perform a (simplified) geometry optimization on the CO molecule and
# calculate the MOHP coefficients for the molecule.
#

def main():
    # find optimum for CO molecule
    res = scipy.optimize.minimize(optimize_lih, 1.14, tol=1e-4)
    print('Optimal distance found at d = %f' % res.x[0])
    
    # calculate sto-3g coefficients for h2o
    result = calculate_lih(res.x)
    
    energies = result['orbe']
    coeff = result['orbc']
    H = result['fock']
    
    mohp = np.zeros(len(coeff))
    for k in range(0, len(coeff)): # loop over molecular orbitals
        for i in range(0,len(coeff)-1): # loop over orbitals on Li
            for j in range(len(coeff)-1,len(coeff)): # loop over orbitals on H
                mohp[k] += 2.0 * H[i,j] * coeff[i,k] * coeff[j,k]
                
    for i,(chi,e) in enumerate(zip(mohp, energies)):
        print('%02i %+8.4f %+8.4f' % (i+1, e, chi))

    # labels = [
    #     'Li,1s',
    #     'Li,2s',
    #     'Li,2p$_{x}$',
    #     'Li,2p$_{y}$',
    #     'Li,2p$_{z}$',
    #     'H,1s'
    # ]

    # # fig, ax = plt.subplots(1,2, dpi=300)
    # # plot_matrix(ax[0], H, xlabels=labels, ylabels=labels,
    # #             title='Hamiltonian matrix')
    # # plot_matrix(ax[1], coeff, ylabels=labels, title='Coefficient matrix',
    # #             orbe=energies)

def optimize_lih(d):
    """
    Optimization function for scipy.optimize.minimize
    """
    mol = Molecule()
    mol.add_atom('Li', 0.0, 0.0, -d[0]/2, unit='angstrom')
    mol.add_atom('H', 0.0, 0.0,  d[0]/2, unit='angstrom')
    
    result = HF().rhf(mol, 'sto3g')
    
    return result['energy']

def calculate_lih(d):
    """
    Full function for evaluation
    """
    mol = Molecule()
    mol.add_atom('Li', 0.0, 0.0, -d[0]/2, unit='angstrom')
    mol.add_atom('H', 0.0, 0.0,  d[0]/2, unit='angstrom')
    
    result = HF().rhf(mol, 'sto3g')
    
    return result

def plot_matrix(ax, mat, xlabels=None, ylabels=None, title = None, orbe=None):
    """
    Produce plot of matrix
    """
    ax.imshow(mat, vmin=-1, vmax=1, cmap='PiYG')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(i, j, '%.3f' % mat[j,i], ha='center', va='center',
                    fontsize=7, zorder=5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.hlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    ax.vlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    
    # add xlabels if available
    if xlabels:
        ax.set_xticks(np.arange(0, mat.shape[0]))
        ax.set_xticklabels(ylabels)
        ax.tick_params(axis='both', which='major', labelsize=7)
    
    # add basis functions as axes labels
    if ylabels:
        ax.set_yticks(np.arange(0, mat.shape[0]))
        ax.set_yticklabels(ylabels)
        ax.tick_params(axis='both', which='major', labelsize=7)
    
    if orbe is not None:
        orbe = ['(%i) / %.4f' % ((i+1),e) for i,e in enumerate(orbe)]
        ax.set_xticks(np.arange(0, mat.shape[0]))
        ax.set_xticklabels(orbe, rotation=90)
    
    # add title if supplied
    if title:
        ax.set_title(title)

if __name__ == '__main__':
    main()
