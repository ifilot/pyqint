# -*- coding: utf-8 -*-

"""
Geometry optimization using Hartree–Fock energies and forces.

This module defines a stateful geometry optimizer that operates on
a single Molecule and basis set, using SciPy's conjugate-gradient
minimizer.
"""

from __future__ import annotations

import time
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import numpy.typing as npt
import scipy.optimize

from . import HF, Molecule


BOHR_TO_ANGSTROM: float = 0.52917721092

Vec = npt.NDArray[np.float64]
Coords = npt.NDArray[np.float64]


class GeometryOptimization:
    """
    Perform geometry optimization for a single molecule.

    The optimizer owns:
      - the molecule being optimized
      - the basis set
      - cached SCF quantities used for acceleration
      - histories of energies, forces, and coordinates

    A GeometryOptimization instance is intended to be used for
    *one molecule only*.
    """

    def __init__(
        self,
        mol: Molecule,
        basis: str,
        *,
        verbose: bool = False,
    ) -> None:
        self.mol: Molecule = mol
        self.basis: str = basis
        self.verbose: bool = verbose

        # SCF cache (used to accelerate convergence)
        self.cinit: Optional[Vec] = None
        self.P: Optional[Vec] = None
        self.orbe: Optional[Vec] = None

        # Optimization bookkeeping
        self.iter: int = 0
        self.coord: Optional[Vec] = None
        self.forces: Optional[Coords] = None
        self.last_energy_run: Optional[Dict[str, Any]] = None

        # History tracking
        self.energies_history: List[float] = []
        self.forces_history: List[Coords] = []
        self.coordinates_history: List[Coords] = []

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, gtol: float = 1e-5) -> Dict[str, Any]:
        """
        Run the geometry optimization.

        Parameters
        ----------
        gtol : float
            Gradient convergence threshold passed to SciPy.

        Returns
        -------
        dict
            Optimization results and history.
        """
        x0 = self._unpack_coords(self.mol)
        self.iter = 0

        # Reset histories
        self.energies_history.clear()
        self.forces_history.clear()
        self.coordinates_history.clear()

        if self.verbose:
            self._print_break("#", newline=False)
            print(" START GEOMETRY OPTIMIZATION (CONJUGATE GRADIENT)")
            self._print_break("#")

        res_opt = scipy.optimize.minimize(
            self._energy,
            x0,
            method="CG",
            jac=self._jacobian,
            options={"gtol": gtol},
        )

        return {
            "opt": res_opt,
            "energies": self.energies_history,
            "forces": self.forces_history,
            "coordinates": self.coordinates_history,
            "data": self.last_energy_run,
            "mol": self.mol,
        }

    # ------------------------------------------------------------------
    # SciPy callbacks
    # ------------------------------------------------------------------

    def _energy(self, x: Vec) -> float:
        """
        Compute total electronic energy for the given coordinates.

        This function is called repeatedly by SciPy.
        """
        start = time.perf_counter()
        self.iter += 1

        # Build updated molecule geometry
        mol = self._pack_coords(self.mol, x)

        if self.verbose:
            self._print_break("=", newline=False)
            print(f"  GEOMETRY OPTIMIZATION STEP {self.iter:03d}")
            self._print_break("=")

        res = HF(mol,self.basis,).rhf(
            orbc_init=self.cinit,
            calc_forces=True,
        )

        # Cache SCF results
        self.cinit = res["orbc"]
        self.P = res["density"]
        self.orbe = res["orbe"]
        self.forces = res["forces"]
        self.coord = x.copy()
        self.last_energy_run = res

        # Record history
        self.energies_history.append(res["energies"][-1])
        self.forces_history.append(self.forces)
        self.coordinates_history.append(
            x.reshape(len(mol.get_atoms()), 3)
        )

        if self.verbose:
            self._print_energies(res)
            self._print_atoms(mol, self.forces)
            elapsed = time.perf_counter() - start
            print(f"\nElapsed time: {elapsed:.4f} s\n")
            self._print_break()

        return res["energies"][-1]

    def _jacobian(self, x: Vec) -> Vec:
        """
        Return forces (negative gradient) for the optimizer.

        If forces were already computed at these coordinates,
        they are reused to avoid an extra SCF call.
        """
        if self.coord is not None and np.max(np.abs(x - self.coord)) < 1e-5:
            assert self.forces is not None
            return self.forces.flatten()

        # Fallback: compute energy + forces explicitly
        mol = self._pack_coords(self.mol, x)
        res = HF(mol, self.basis,).rhf(
            orbc_init=self.cinit,
            calc_forces=True,
        )

        self.cinit = res["orbc"]
        self.P = res["density"]
        self.orbe = res["orbe"]
        self.forces = res["forces"]
        self.coord = x.copy()

        return self.forces.flatten()

    # ------------------------------------------------------------------
    # Coordinate handling
    # ------------------------------------------------------------------

    def _unpack_coords(self, mol: Molecule) -> Vec:
        """
        Extract Cartesian coordinates from a Molecule into a flat array.
        """
        coords = [pos for _, pos in mol.get_atoms()]
        return np.asarray(coords, dtype=np.float64).flatten()

    def _pack_coords(self, mol: Molecule, coords: Vec) -> Molecule:
        """
        Create a new Molecule with updated Cartesian coordinates.
        """
        newmol = Molecule(mol.get_name())
        newmol.set_charge(mol.get_charge())

        coords = coords.reshape(len(mol.get_atoms()), 3)
        for (symbol, _), (x, y, z) in zip(mol.get_atoms(), coords):
            newmol.add_atom(symbol, x, y, z)

        return newmol

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def write_multiframe_xyz(
        self,
        filename: str,
        comment_fmt: Optional[Callable[[int, Coords], str]] = None,
    ) -> None:
        """
        Write the optimization trajectory to a multi-frame XYZ file.
        """
        atoms = self.mol.get_atoms()
        symbols = [sym for sym, _ in atoms]
        natoms = len(symbols)

        with open(filename, "w") as f:
            for iframe, coords in enumerate(self.coordinates_history):
                coords_ang = coords * BOHR_TO_ANGSTROM
                f.write(f"{natoms}\n")
                comment = (
                    comment_fmt(iframe, coords_ang)
                    if comment_fmt
                    else f"frame={iframe}"
                )
                f.write(comment + "\n")

                for sym, (x, y, z) in zip(symbols, coords_ang):
                    f.write(f"{sym:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")

    # ------------------------------------------------------------------
    # Pretty-printing helpers
    # ------------------------------------------------------------------

    def _print_atoms(self, mol: Molecule, forces: Coords) -> None:
        self._print_break("-", n=80, newline=False)
        print("    POSITIONS AND FORCES")
        self._print_break("-", n=80)

        for (sym, pos), force in zip(mol.get_atoms(), forces):
            pos_ang = pos * BOHR_TO_ANGSTROM
            force_ang = force * BOHR_TO_ANGSTROM
            print(
                f"  {sym:2s} | "
                f"{pos_ang[0]:10.6f} {pos_ang[1]:10.6f} {pos_ang[2]:10.6f} | "
                f"{force_ang[0]:+10.4e} {force_ang[1]:+10.4e} {force_ang[2]:+10.4e}"
            )

    def _print_energies(self, res: Dict[str, Any]) -> None:
        self._print_break("-", n=80, newline=False)
        print("    ENERGIES")
        self._print_break("-", n=80)
        print(f"  Kinetic:                     {res['ekin']:12.8f}")
        print(f"  Nuclear:                     {res['enuc']:12.8f}")
        print(f"  Electron-electron repulsion: {res['erep']:12.8f}")
        print(f"  Exchange:                    {res['ex']:12.8f}")
        print(f"  Nuclear repulsion:           {res['enucrep']:12.8f}")
        print(f"  TOTAL:                       {res['energies'][-1]:12.8f}")
        print() # add newline

    def _print_break(self, ch: str = "=", n: int = 80, newline: bool = True) -> None:
        print(ch * n)
        if newline:
            print()