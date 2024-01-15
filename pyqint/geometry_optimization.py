# -*- coding: utf-8 -*-

from . import HF, Molecule
import numpy as np
import scipy.optimize

class GeometryOptimization:
    """
    Class to perform geometry optimizaton using the Conjugate Gradient algorithm
    as implemented in the scipy library
    """
    def __init__(self, verbose=False):
        self.cinit = None
        self.P = None
        self.orbe = None
        self.verbose = verbose
        self.iter = 0
        self.energies_history = []
        self.forces_history = []
        self.coordinates_history = []
        self.coord = None

    def run(self, mol, basis, gtol=1e-5):
        """
        Perform geometry optimization
        """
        x0 = self.__unpack(mol) # unpack coordinates from molecule class
        self.iter = 0

        # reset arrays
        self.energies_history = []
        self.forces_history = []
        self.coordinates_history = []

        if self.verbose:
            self.__print_break(newline=False)
            print("START GEOMETRY OPTIMIZATION")
            print("USING CONJUGATE GRADIENT PROCEDURE")
            self.__print_break(newline=True)

        res_opt = scipy.optimize.minimize(self.energy, x0, args=(mol, basis), method='CG',
                                          jac=self.jacobian, options={'gtol':gtol})

        res = {
            'opt': res_opt,
            'energies': self.energies_history,
            'forces': self.forces_history,
            'coordinates': self.coordinates_history,
            'data': self.last_energy_run
        }

        return res

    def energy(self, x, mol, basis):
        self.iter += 1

        mol = self.__pack(mol, x) # pack positions into new molecule class

        if self.verbose:
            self.__print_break(newline=False)
            print('  START GEOMETRY OPTIMIZATION STEP %03i' % self.iter)
            self.__print_break()
            self.__print_positions(mol)
            print() # newline

        res = HF().rhf(mol, basis, orbc_init=self.cinit, calc_forces=True)

        if self.verbose:
            self.__print_energies(res)
            print() # print newline

        # cache matrices
        self.cinit = res['orbc']
        self.P = res['density']
        self.orbe = res['orbe']
        self.forces = res['forces']
        self.coord = x
        self.last_energy_run = res

        # append history
        self.forces_history.append(self.forces)
        self.coordinates_history.append(self.coord.reshape(len(mol.get_atoms()),3))
        self.energies_history.append(res['energies'][-1])

        if self.verbose:
            self.__print_forces(mol, self.forces)
            print() # print newline
            self.__print_break(newline=False)
            print('  END GEOMETRY OPTIMIZATION STEP %03i' % self.iter)
            self.__print_break(newline=True)

        return res['energies'][-1]

    def jacobian(self, x, mol, basis):
        """
        Calculate the forces. Note that the forces are typically already
        calculated in the energy step so here only a simple check is done
        to see if the coordinates match. If so, the forces are simply returned.
        """
        # check if forces are already calculated
        if self.coord is not None and np.max(x - self.coord) < 1e-5:
            return self.forces.flatten()
        else:
            # no forces have yet been calculated, indicative that no energy run has
            # yet been done
            res = HF().rhf(mol, basis, orbc_init=self.cinit, calc_forces=True)
            self.cinit = res['orbc']
            self.P = res['density']
            self.orbe = res['orbe']
            self.forces = res['forces']
            self.coord = x
            return res['forces'].flatten()

    def __unpack(self, mol):
        """
        Unpack coordinates from molecule class
        """
        coords = []
        for atom in mol.get_atoms():
            coords.append(atom[1])

        return np.array(coords).flatten()

    def __pack(self, mol, coords):
        """
        Pack coordinates into new molecule class
        """

        newmol = Molecule(mol.name)
        newmol.set_charge(mol.get_charge())
        coords = coords.reshape((len(mol.get_atoms()), 3))
        for i,atom in enumerate(mol.get_atoms()):
            newmol.add_atom(mol.get_atoms()[i][0], coords[i][0], coords[i][1], coords[i][2])

        return newmol

    def __print_positions(self, mol):
        """
        Print atomic positions in nice formatting
        """
        print('-------------')
        print('  POSITIONS  ')
        print('-------------')
        for atom in mol.get_atoms():
            print('  %2s %12.8f %12.8f %12.8f' % (atom[0],
                                                  atom[1][0],
                                                  atom[1][1],
                                                  atom[1][2]))

    def __print_forces(self, mol, forces):
        """
        Print forces using nice formatting
        """
        print('----------')
        print('  FORCES  ')
        print('----------')

        for atom,force in zip(mol.get_atoms(), forces):
            print('  %2s %12.4e %12.4e %12.4e' % (atom[0],
                                                  force[0],
                                                  force[1],
                                                  force[2]))

    def __print_energies(self, res):
        """
        Print the energy terms in each iteration
        """
        print('------------')
        print('  ENERGIES  ')
        print('------------')

        print('  Kinetic:                     %12.8f' % res['ekin'])
        print('  Nuclear:                     %12.8f' % res['enuc'])
        print('  Electron-electron repulsion: %12.8f' % res['erep'])
        print('  Exchange:                    %12.8f' % res['ex'])
        print('  Nuclear repulsion:           %12.8f' % res['enucrep'])
        print('  TOTAL:                       %12.8f' % res['energies'][-1])

    def __print_break(self, newline=True):
        print('=' * 80)
        if newline:
            print()
