# -*- coding: utf-8 -*-

from .pyqint import PyQInt
import numpy as np
import random
import scipy.optimize

class FosterBoys:
    """
    Routine for constructing localized orbitals using Foster-Boys procedure
    """
    def __init__(self, res):
        # copy objects from Hartree-Fock result dictionary
        self.orbc_canonical = res['orbc']
        self.orbe_canonical = res['orbe']
        self.nelec = res['nelec']
        self.H = res['fock']
        self.cgfs = res['cgfs']
        self.maxiter = 1000
        self.occ = [1 if i < self.nelec//2 else 0 for i in range(0, len(self.cgfs))]
        
        # construct dipole tensor
        self.dipol = self.__build_dipole_tensor(self.cgfs)
    
    def run(self):
        """
        Perform the Foster-Boys optimization routine
        """
        # start with a random unitary transformation
        C = self.__construct_random_orthogonal_transformation(self.orbc_canonical, 
                                                              self.nelec)

        old_r2 = 0.0 # start assuming orbital centroids lie at origin
        diff = 1
        nriter = 0
        while diff > 1e-7 and nriter < self.maxiter:
            C, r2 = self.__perform_foster_boys_mixing_scipy(C)
            diff = np.abs(r2 - old_r2)
            old_r2 = r2
            nriter += 1
        
        if nriter == self.maxiter:
            raise Exception('Foster-Boys procedure did not converge.')

        orbe, orbc = self.__calculate_molecular_orbital_energies(C)

        result = {
            'orbe': orbe,
            'orbc': orbc,
            'nriter': nriter,
            'r2start': self.__calculate_r2(self.orbc_canonical),
            'r2final': self.__calculate_r2(orbc)
        }

        return result
    
    def __perform_foster_boys_mixing_scipy(self, C):
        """
        Successively perform n(n-1)/2 optimization among the occupied
        molecular orbitals using the scipy minimize function
        """
        r2start = self.__calculate_r2(C)
        r2final = r2start
        
        N = self.nelec // 2
        for i in range(N):
            for j in range(i+1, N):
                
                res = scipy.optimize.minimize(self.__evaluate_orbital_centroids, 
                                              np.pi/4, 
                                              args=(C, i, j))
                
                alpha = res.x[0]
                Cnew = C.copy()
                Cnew[:,i] = np.cos(alpha) * C[:,i] + np.sin(alpha) * C[:,j]
                Cnew[:,j] = -np.sin(alpha) * C[:,i] + np.cos(alpha) * C[:,j]
                
                r2 = self.__calculate_r2(Cnew)
                if(r2 > r2start):
                    C = Cnew.copy()
                    r2final = r2start
        
        return C, r2final
    
    def __evaluate_orbital_centroids(self, alpha, C, i, j):
        """
        Evaluate the result of r2 after a 2x2 rotation
        between solutions i and j
        """
        
        Cnew = C.copy()
        Cnew[:,i] = np.cos(alpha) * C[:,i] + np.sin(alpha) * C[:,j]
        Cnew[:,j] = -np.sin(alpha) * C[:,i] + np.cos(alpha) * C[:,j]
        r2 = self.__calculate_r2(Cnew)
        
        return -r2
    
    def __calculate_r2(self, C):
        """
        Calculate r2 from coefficient matrix
        """
        # use efficient Einstein summation
        # i -> molecular orbitals
        # j -> basis functions 1
        # k -> basis functions 2
        # l -> cartesian directions dipole moment
        dipolest = np.einsum('ji,ki,jkl->il', C, C, self.dipol)
        r2 = np.einsum('ij,i->', dipolest**2, self.occ)
        
        return r2
    
    def __construct_random_orthogonal_transformation(self, C, nelec, nops=100):
        """
        Construct unitary transformation matrix via a series of two-dimensional
        unitary transformations among the occupied molecular orbitals
        """
        for i in range(0,nops):
            Cnew = C.copy()
            # randomly pick two numbers among the occupied orbitals
            n = random.sample(range(nelec//2),2)
            
            gamma = random.random() * 2.0 * np.pi
            
            for j in range(0,len(C)):
                Cnew[j,n[0]] = np.cos(gamma) * C[j,n[0]] + np.sin(gamma) * C[j,n[1]]
                Cnew[j,n[1]] = -np.sin(gamma) * C[j,n[0]] + np.cos(gamma) * C[j,n[1]]

            C = Cnew.copy()

        return C
    
    def __build_dipole_tensor(self, cgfs):
        """
        Build a dipole tensor
        """
        # create new cgfs
        N = len(cgfs)
        mat = np.zeros((N,N,3))
        integrator = PyQInt()
        for i,cgf1 in enumerate(cgfs):
            for j,cgf2 in enumerate(cgfs):
                for k in range(0,3):
                    mat[i,j,k] = integrator.dipole(cgf1, cgf2, k, 0.0)
        
        return mat

    def __calculate_molecular_orbital_energies(self, C):
        """
        Calculate the one-electron MO energies from the Hamiltonian matrix
        and the coefficient matrix, both in their original basis

        Return *ordered* list of eigenvalue and -vector pairs
        """
        orbe = np.zeros(len(C))
        for i in range(len(C)):
            orbe[i] = C[:,i].dot(self.H.dot(C[:,i]))

        # produce list of indices for eigenvalues in ascending order
        oidx = np.argsort(orbe)

        return orbe[oidx], C[:,oidx]
