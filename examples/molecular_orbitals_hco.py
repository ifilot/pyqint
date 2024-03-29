from pyqint import Molecule, HF, PyQInt, GeometryOptimization
import numpy as np
import matplotlib.pyplot as plt

#
# Plot the isosurfaces for the CO molecule
#

def main():
    # calculate sto-3g coefficients for co
    result = optimize_hco()
    energies = result['orbe']
    coeff = result['orbc']
        
    fig, ax = plt.subplots(2, 5, dpi=144, figsize=(12,5))           
    for i,(chi,e) in enumerate(zip(coeff.transpose(), energies)):
        
        if i >= 10:
            break
        
        res, x, y = build_contourplot(result['cgfs'], chi, sz=3, plane='xz')
        vmax = np.max(np.abs(res))
        vmin = -vmax
        ax[i//5,i%5].contourf(x, y, res, levels=15, cmap='PiYG',
                              vmin=vmin, vmax=vmax)
        ax[i//5,i%5].contour(x, y, res, levels=15, colors='black',
                             vmin=vmin, vmax=vmax)
        ax[i//5,i%5].set_xlabel('x [Bohr]')
        ax[i//5,i%5].set_ylabel('z [Bohr]')
        ax[i//5,i%5].set_title('Energy = %.4f Ht' % e)
        ax[i//5,i%5].grid(linestyle='--', color='black', alpha=0.5)
        
    plt.tight_layout()
    plt.savefig('hco.jpg')

def optimize_hco():
    """
    Optimization function for scipy.optimize.minimize
    """
    mol = Molecule()
    mol.add_atom('O', -0.12770367,  -0.00000000,   1.23419284)
    mol.add_atom('C', -0.50625665,   0.00000000,  -1.09149431)
    mol.add_atom('H',  1.57882331,  -0.00000000,  -2.06681794)
    mol.set_charge(-1)
    
    res = GeometryOptimization().run(mol, 'p321')
    
    return res['data']

def build_contourplot(cgfs, coeff, sz=2, npts=50, plane='xy'):
    integrator = PyQInt()
    
    # build grid
    x = np.linspace(-sz, sz, 50)
    y = np.linspace(-sz, sz, 50)
    xx, yy = np.meshgrid(x,y)
    zz = np.zeros(len(x) * len(y))
    
    if plane == 'xy':
        points = [xx.flatten(), yy.flatten(), zz]
    elif plane == 'xz':
        points = [xx.flatten(), zz, yy.flatten()]
    elif plane == 'yz':
        points = [zz, xx.flatten(), yy.flatten()]
    else:
        raise Exception('Unknown plane: %s' % plane)
    
    grid = np.vstack(points).reshape(3,-1).T
    res = integrator.plot_wavefunction(grid, coeff, cgfs).reshape((len(y), len(x)))
    
    return res, x, y

if __name__ == '__main__':
    main()