# -*- coding: utf-8 -*-

import json
import os
from .pyqint import cgf

class Molecule:
    """
    Molecule class
    """
    def __init__(self):
        self.atoms = []
        self.charges = []

    def add_atom(self, atom, x, y, z):
        self.atoms.append([atom, x, y, z])
        self.charges.append(0)

    def build_basis(self, name):
        self.cgfs = []
        self.nuclei = []

        basis_filename = os.path.join(os.path.dirname(__file__), 'basis', '%s.json' % name)
        f = open(basis_filename, 'r')
        basis = json.load(f)
        f.close()

        for aidx, atom in enumerate(self.atoms):
            cgfs_template = basis[atom[0]]

            # store information about the nuclei
            self.charges[aidx] = cgfs_template['atomic_number']
            self.nuclei.append([[atom[1], atom[2], atom[3]], cgfs_template['atomic_number']])

            for cgf_t in cgfs_template['cgfs']:
                # s-orbitals
                if cgf_t['type'] == 'S':
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 0, 0)

                # p-orbitals
                if cgf_t['type'] == 'P':
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 1, 0, 0)
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 1, 0)
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 0, 1)

                # d-orbitals
                if cgf_t['type'] == 'D':
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 2, 0, 0)
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 2, 0)
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 0, 2)
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 1, 1, 0)
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 1, 0, 1)
                    self.cgfs.append(cgf([atom[1], atom[2], atom[3]]))
                    for gto in cgf_t['gtos']:
                        self.cgfs[-1].add_gto(gto['coeff'], gto['alpha'], 0, 1, 1)

        return self.cgfs, self.nuclei
