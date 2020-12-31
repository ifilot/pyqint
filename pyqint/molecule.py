# -*- coding: utf-8 -*-

class Molecule:
    """
    Molecule class
    """
    def __init__(self):
        self.atoms = []

    def add_atom(self, atom, x, y, z):
        self.atoms.append([atom, x, y, z])
