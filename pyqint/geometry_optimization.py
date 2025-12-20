# -*- coding: utf-8 -*-

from . import HF, Molecule
import numpy as np
import scipy.optimize
import time

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
        self.mol = mol
        self.iter = 0

        # reset arrays
        self.energies_history = []
        self.forces_history = []
        self.coordinates_history = []

        if self.verbose:
            self.__print_break(ch='#', newline=False)
            print(" START GEOMETRY OPTIMIZATION: USING CONJUGATE GRADIENT PROCEDURE")
            self.__print_break(ch='#', newline=True)

        res_opt = scipy.optimize.minimize(self.energy, x0, args=(mol, basis), method='CG',
                                          jac=self.jacobian, options={'gtol':gtol})

        res = {
            'opt': res_opt,
            'energies': self.energies_history,
            'forces': self.forces_history,
            'coordinates': self.coordinates_history,
            'data': self.last_energy_run,
            'mol': self.mol
        }

        return res

    def energy(self, x, mol, basis):
        st = time.perf_counter()
        self.iter += 1
        mol = self.__pack(mol, x) # pack positions into new molecule class

        if self.verbose:
            self.__print_break(ch='=', newline=False)
            print('  START GEOMETRY OPTIMIZATION STEP %03i' % self.iter)
            self.__print_break(ch='=')

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

        # keep track of time
        elapsed = (time.perf_counter() - st)

        if self.verbose:
            self.__print_atoms(mol, self.forces)
            print()
            print('Elapsed time: %.4f seconds' % elapsed)
            print()
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

    def __print_atoms(self, mol, forces):
        """
        Print atomic positions in nice formatting
        """
        self.__print_break(ch='-', n=80, newline=False)
        print('    POSITIONS AND FORCES')
        self.__print_break(ch='-', n=80, newline=False)
        for atom,force in zip(mol.get_atoms(), forces):
            print('  %2s | %10.6f %10.6f %10.6f | %+10.4e %+10.4e %+10.4e' % 
                  (atom[0], atom[1][0], atom[1][1], atom[1][2],
                   force[0], force[1], force[2]))

    def __print_energies(self, res):
        """
        Print the energy terms in each iteration
        """
        self.__print_break(ch='-', n=80, newline=False)
        print('    ENERGIES')
        self.__print_break(ch='-', n=80, newline=False)

        print('  Kinetic:                     %12.8f' % res['ekin'])
        print('  Nuclear:                     %12.8f' % res['enuc'])
        print('  Electron-electron repulsion: %12.8f' % res['erep'])
        print('  Exchange:                    %12.8f' % res['ex'])
        print('  Nuclear repulsion:           %12.8f' % res['enucrep'])
        print('  TOTAL:                       %12.8f' % res['energies'][-1])

    def __print_break(self, ch='=', n = 80, newline=True):
        print(ch * n)
        if newline:
            print()

    def write_multiframe_xyz(self, filename, comment_fmt=None):
        """
        Write a multi-frame XYZ file from a list of coordinate arrays.

        Parameters
        ----------
        filename : str
            Output XYZ file name.
        comment_fmt : callable, optional
            Function that takes (frame_index, coords) and returns a comment string
            for the XYZ frame header.
            Example:
                lambda i, c: f"step={i}"
        """

        BOHR_TO_ANGSTROM = 0.52917721092

        atoms = self.mol.get_atoms()
        symbols = [atom[0] for atom in atoms]
        natoms = len(symbols)
        coords_list = self.coordinates_history

        # Basic validation
        for i, coords in enumerate(coords_list):
            if coords.shape != (natoms, 3):
                raise ValueError(
                    f"Frame {i} has shape {coords.shape}, expected ({natoms}, 3)"
                )

        with open(filename, "w") as f:
            for iframe, coords in enumerate(coords_list):
                coords = coords * BOHR_TO_ANGSTROM
                f.write(f"{natoms}\n")

                if comment_fmt is None:
                    f.write(f"frame={iframe}\n")
                else:
                    f.write(comment_fmt(iframe, coords) + "\n")

                # Atomic coordinates
                for sym, (x, y, z) in zip(symbols, coords):
                    f.write(f"{sym:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")
