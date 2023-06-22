from pyqint import Molecule, HF
import numpy as np
import scipy.optimize

#
# Perform a (simplified) geometry optimization on the CO molecule and
# calculate the COHP coefficients for the molecule.
#

def main():
    # find optimum for CO molecule
    res = scipy.optimize.minimize(optimize_lih, [1.14], tol=1e-4)
    print('Optimal distance found at d = %f' % res.x)
    
    # calculate sto-3g coefficients for h2o
    result = calculate_lih(res.x)
    
    energies = result['orbe']
    coeff = result['orbc']
    H = result['fock']
    
    cohp = np.zeros(len(coeff))
    for k in range(0, len(coeff)):
        for i in range(0,len(coeff)//2): # loop over orbitals on Li
            for j in range(len(coeff)//2,len(coeff)): # loop over orbitals on H
                cohp[k] += 2.0 * H[i,j] * coeff[i,k] * coeff[j,k]
                
    for i,(chi,e) in enumerate(zip(cohp, energies)):
        print('%02i %+8.4f %+8.4f' % (i+1, e, chi))

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

if __name__ == '__main__':
    main()
