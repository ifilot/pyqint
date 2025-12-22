# -*- coding: utf-8 -*-

"""
Molecule representation and basis construction.

This module defines the Molecule class, which stores atomic geometry, handles
unit conversion, constructs Gaussian basis functions from JSON basis-set
definitions, and exposes nuclear and electronic properties required for
quantum-chemical calculations.
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple, Optional

import numpy as np
import numpy.typing as npt

from .cgf import CGF
from .element import Element


# --- Type aliases ------------------------------------------------------------

Vec3 = npt.NDArray[np.float64]
AtomEntry = Tuple[str, Vec3]        # (element symbol, position in bohr)
NucleusEntry = Tuple[Vec3, int]     # (position, nuclear charge)


class Molecule:
    """
    Container class describing a molecule.

    A Molecule consists of:
      - a list of atoms with positions (in bohr)
      - nuclear charges
      - an optional total molecular charge
      - a set of contracted Gaussian functions (CGFs) once a basis is built

    The class itself is intentionally lightweight; most heavy computation
    happens in CGF/GTO objects and downstream integrals.
    """

    def __init__(self, name: str = "unknown") -> None:
        # Atomic symbols and positions
        self.__atoms: List[AtomEntry] = []

        # Nuclear charges (one per atom, populated later)
        self.__charges: List[int] = []

        # Human-readable molecule name
        self.__name: str = name

        # Net molecular charge (used when computing electron count)
        self.__charge: int = 0

        # Total number of electrons (computed after basis/nuclei setup)
        self.__nelec: Optional[int] = None

        # Populated lazily
        self.__cgfs: List[CGF]
        self.__nuclei: List[NucleusEntry]

    def __str__(self) -> str:
        """
        Human-readable string representation.
        """
        res = f"Molecule: {self.__name}\n"
        for symbol, pos in self.__atoms:
            res += f" {symbol} ({pos[0]:.6f},{pos[1]:.6f},{pos[2]:.6f})\n"
        return res

    def __len__(self) -> int:
        """
        Number of atoms in the molecule.
        """
        return len(self.__atoms)

    # --- Basic properties ----------------------------------------------------

    def get_nelec(self) -> int:
        """
        Return the number of electrons in the molecule.

        This value is only defined after nuclear charges have been
        determined (via build_basis() or get_nuclei()).
        """
        if self.__nelec is None:
            raise RuntimeError(
                "Number of electrons not initialized. "
                "Call build_basis() or get_nuclei() first."
            )
        return self.__nelec - self.__charge

    def get_atoms(self) -> List[AtomEntry]:
        """
        Return the list of atoms as (symbol, position) tuples.
        """
        return self.__atoms

    def get_charge(self) -> int:
        """
        Return the net molecular charge.
        """
        return self.__charge
    
    def get_name(self) -> str:
        """
        Return the name of the molecule
        """
        return self.__name

    def set_charge(self, charge: int) -> None:
        """
        Set the net molecular charge.

        Positive values correspond to cations, negative values to anions.
        """
        self.__charge = charge

    # --- Geometry handling ---------------------------------------------------

    def add_atom(
        self,
        atom: str,
        x: float,
        y: float,
        z: float,
        unit: str = "bohr",
    ) -> None:
        """
        Add an atom to the molecule.

        Parameters
        ----------
        atom : str
            Atomic symbol (e.g. "H", "C", "O").
        x, y, z : float
            Cartesian coordinates.
        unit : {"bohr", "angstrom"}
            Unit of the supplied coordinates.

        Notes
        -----
        Internally, all coordinates are stored in bohr.
        """
        ang2bohr = 1.8897259886

        x, y, z = float(x), float(y), float(z)

        if unit == "bohr":
            pos = np.array([x, y, z], dtype=np.float64)
        elif unit == "angstrom":
            pos = np.array(
                [x * ang2bohr, y * ang2bohr, z * ang2bohr],
                dtype=np.float64,
            )
        else:
            raise RuntimeError(
                f"Invalid unit '{unit}'. "
                "Accepted units are 'bohr' and 'angstrom'."
            )

        # Nuclear charge is populated later
        self.__atoms.append((atom, pos))
        self.__charges.append(0)

    # --- Basis construction --------------------------------------------------

    def build_basis(self, name: str) -> Tuple[List[CGF], List[NucleusEntry]]:
        """
        Build a Gaussian basis set from a JSON basis definition.

        Parameters
        ----------
        name : str
            Basis-set label (corresponding to a JSON file in basissets/).

        Returns
        -------
        cgfs : list[CGF]
            Contracted Gaussian functions for the molecule.
        nuclei : list[(Vec3, int)]
            Nuclear positions and charges.
        """
        basis_filename = os.path.join(
            os.path.dirname(__file__),
            "basissets",
            f"{name}.json",
        )

        with open(basis_filename, "r") as f:
            basis = json.load(f)

        self.__cgfs = []

        # Loop over atoms and attach basis functions
        for aidx, (symbol, position) in enumerate(self.__atoms):
            atom_basis = basis[symbol]

            # Store nuclear charge
            self.__charges[aidx] = atom_basis["atomic_number"]

            for cgf_def in atom_basis["cgfs"]:

                # s orbitals
                if cgf_def["type"] == "S":
                    cgf = CGF(position)
                    for gto in cgf_def["gtos"]:
                        cgf.add_gto(gto["coeff"], gto["alpha"], 0, 0, 0)
                    self.__cgfs.append(cgf)

                # p orbitals (x, y, z)
                elif cgf_def["type"] == "P":
                    for l, m, n in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
                        cgf = CGF(position)
                        for gto in cgf_def["gtos"]:
                            cgf.add_gto(gto["coeff"], gto["alpha"], l, m, n)
                        self.__cgfs.append(cgf)

                # d orbitals (Cartesian representation)
                elif cgf_def["type"] == "D":
                    for l, m, n in [
                        (2, 0, 0),
                        (0, 2, 0),
                        (0, 0, 2),
                        (1, 1, 0),
                        (1, 0, 1),
                        (0, 1, 1),
                    ]:
                        cgf = CGF(position)
                        for gto in cgf_def["gtos"]:
                            cgf.add_gto(gto["coeff"], gto["alpha"], l, m, n)
                        self.__cgfs.append(cgf)

        # Finalize nuclear information and electron count
        self.get_nuclei()
        self.__nelec = int(np.sum(self.__charges))

        return self.__cgfs, self.__nuclei

    # --- Nuclear information -------------------------------------------------

    def get_nuclei(self) -> List[NucleusEntry]:
        """
        Return nuclear positions and charges.

        This also initializes the total electron count, assuming
        neutral atoms prior to applying molecular charge.
        """
        el = Element()
        self.__nuclei = []

        for aidx, (symbol, position) in enumerate(self.__atoms):
            self.__charges[aidx] = int(getattr(el, symbol))
            self.__nuclei.append((position, self.__charges[aidx]))

        self.__nelec = int(np.sum(self.__charges))
        return self.__nuclei