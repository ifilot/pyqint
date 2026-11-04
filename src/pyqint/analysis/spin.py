# -*- coding: utf-8 -*-

"""
Helpers for handling restricted (RHF) and unrestricted (UHF) Hartree-Fock
result dictionaries in a uniform way.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

SPINS = ("alpha", "beta")


def is_unrestricted(res: Dict[str, Any]) -> bool:
    """
    Return whether a result dictionary originates from a UHF calculation.
    """
    return "orbc_alpha" in res


def get_spin_channel(res: Dict[str, Any], spin: Optional[str] = None) -> Dict[str, Any]:
    """
    Extract the orbital data of a single spin channel.

    Parameters
    ----------
    res
        Result dictionary of an RHF or UHF calculation.
    spin
        For UHF results, either ``'alpha'`` or ``'beta'``. Must be ``None``
        for RHF results.

    Returns
    -------
    dict
        Dictionary with keys

        - ``'orbc'``    : orbital coefficients
        - ``'orbe'``    : orbital energies
        - ``'fock'``    : Fock matrix
        - ``'density'`` : density matrix of this channel (total density
          for RHF)
        - ``'nocc'``    : number of occupied orbitals
        - ``'occ_factor'`` : electrons per occupied orbital (2 for RHF,
          1 for UHF)
    """
    if not is_unrestricted(res):
        if spin is not None:
            raise ValueError(
                "The 'spin' argument is only valid for unrestricted "
                "Hartree-Fock results."
            )
        # use .get() for quantities that are absent in user-constructed
        # result dictionaries (e.g. for plotting purposes only)
        return {
            "orbc": res["orbc"],
            "orbe": res.get("orbe"),
            "fock": res.get("fock"),
            "density": res.get("density"),
            "nocc": res["nelec"] // 2 if "nelec" in res else None,
            "occ_factor": 2.0,
        }

    if spin not in SPINS:
        raise ValueError(
            "For unrestricted Hartree-Fock results, 'spin' must be "
            "either 'alpha' or 'beta'."
        )

    return {
        "orbc": res["orbc_" + spin],
        "orbe": res["orbe_" + spin],
        "fock": res["fock_" + spin],
        "density": res["density_" + spin],
        "nocc": res["nalpha"] if spin == "alpha" else res["nbeta"],
        "occ_factor": 1.0,
    }
