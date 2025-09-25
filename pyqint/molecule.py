# -*- coding: utf-8 -*-

import json
import os
import numpy as np
from .cgf import CGF
from .element import Element

class Molecule:
    """
    Molecule class
    """
    def __init__(self, _name='unknown'):
        self.__atoms = []
        self.__charges = []
        self.name = _name
        self.__charge = 0
        self.__nelec = None

    def __str__(self):
        res = "Molecule: %s\n" % self.name
        for atom in self.__atoms:
            res += " %s (%f,%f,%f)\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2])

        return res
    
    def __len__(self):
        return len(self.__atoms)

    def get_nelec(self):
        """
        Get the number of electrons
        """
        if self.__nelec == None:
            raise Exception('You need to use build_basis() or get_nuclei() before using this function.')
        return self.__nelec - self.__charge

    def get_atoms(self):
        return self.__atoms
    
    def get_charge(self):
        return self.__charge

    def set_charge(self, charge):
        """
        Set the charge of the molecule
        """
        self.__charge = charge

    def add_atom(self, atom, x, y, z, unit='bohr'):
        """
        Add an atom to the molecule
        """
        ang2bohr = 1.8897259886

        x = float(x)
        y = float(y)
        z = float(z)

        if unit == "bohr":
            self.__atoms.append([atom, np.array([x, y, z])])
        elif unit == "angstrom":
            self.__atoms.append([atom, np.array([x*ang2bohr, y*ang2bohr, z*ang2bohr])])
        else:
            raise RuntimeError("Invalid unit encountered: %s. Accepted units are 'bohr' and 'angstrom'." % unit)

        self.__charges.append(0)

    def build_basis(self, name):
        """
        Build a basis set from a label

        Returns list of CGFs and nuclei
        """
        basis_filename = os.path.join(os.path.dirname(__file__), 'basissets', '%s.json' % name)
        f = open(basis_filename, 'r')
        basis = json.load(f)
        f.close()
        
        self.__cgfs = []

        for aidx, atom in enumerate(self.__atoms):
            cgfs_template = basis[atom[0]]

            # store information about the nuclei
            self.__charges[aidx] = cgfs_template['atomic_number']

            for cgf_t in cgfs_template['cgfs']:
                # s-orbitals
                if cgf_t['type'] == 'S':
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 0, 0)

                # p-orbitals
                if cgf_t['type'] == 'P':
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 1, 0, 0)
                    
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 1, 0)
                    
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 0, 1)

                # d-orbitals
                if cgf_t['type'] == 'D':
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 2, 0, 0)
                    
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 2, 0)
                    
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 0, 2)
                    
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 1, 1, 0)
                    
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 1, 0, 1)
                    
                    self.__cgfs.append(CGF(atom[1]))
                    for gto in cgf_t['gtos']:
                        self.__cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 1, 1)

        # build nuclei objects
        self.get_nuclei()
        self.__nelec = np.sum(self.__charges)

        return self.__cgfs, self.__nuclei
    
    def get_nuclei(self):
        """
        Get the nuclei as a packed array
        """
        el = Element()

        # reset list
        self.__nuclei = []

        for aidx, atom in enumerate(self.__atoms):

            # store information about the nuclei
            self.__charges[aidx] = getattr(el, atom[0])
            self.__nuclei.append([atom[1], self.__charges[aidx]])

        # populate number of electrons
        self.__nelec = np.sum(self.__charges)

        return self.__nuclei