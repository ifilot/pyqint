"""Construct symmetry-adapted linear combinations of Gaussian functions."""

from dataclasses import dataclass
from math import fsum, isfinite
from typing import Iterable, Optional

import numpy as np
import numpy.typing as npt

from ..pyqint_core import PyQInt
from .cgf import CGF


@dataclass
class SALCBasis:
    """Normalized CGFs, their effective AO transformation, and optional labels."""

    cgfs: list[CGF]
    transform: npt.NDArray[np.float64]
    labels: Optional[tuple[str, ...]]


def build_salcs(
    aos: Iterable[CGF],
    coefficients: npt.ArrayLike,
    *,
    labels: Optional[Iterable[str]] = None,
) -> SALCBasis:
    """Construct normalized CGFs from rows of real AO coefficients.

    Columns correspond to the normalized CGFs supplied in ``aos``. Each
    row defines one SALC. Optional labels follow the row order. Use the
    returned ``transform`` to transform AO matrices, for example
    ``transform @ overlap_ao @ transform.T``.
    """
    aos = tuple(aos)
    try:
        weights = np.asarray(coefficients)
    except (TypeError, ValueError) as exc:
        raise ValueError("SALC coefficients must form a rectangular matrix.") from exc
    if not aos or weights.ndim != 2 or weights.shape[0] == 0 or weights.shape[1] != len(aos):
        raise ValueError(f"SALC coefficients must have one row per SALC and {len(aos)} columns for the input CGFs.")
    if np.iscomplexobj(weights) or (weights.dtype == object and any(np.iscomplexobj(v) for v in weights.flat)):
        raise ValueError("SALC coefficients must be real and finite.")
    try:
        weights = np.asarray(weights, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("SALC coefficients must be real and finite.") from exc
    if not np.all(np.isfinite(weights)):
        raise ValueError("SALC coefficients must be real and finite.")

    if labels is not None:
        if isinstance(labels, str):
            raise ValueError("Supply one string label per SALC.")
        labels = tuple(labels)
        if len(labels) != len(weights) or any(not isinstance(label, str) for label in labels):
            raise ValueError("Supply one string label per SALC.")

    norms = [PyQInt._contraction_norm(ao) for ao in aos]
    for i, norm in enumerate(norms):
        if not isfinite(norm) or norm <= 0:
            raise ValueError(f"Input CGF {i} has zero or invalid norm.")

    cgfs = []
    transform = np.empty_like(weights)
    roundoff = 8 * np.finfo(np.float64).eps
    for i, weights_row in enumerate(weights):
        row_name = f"row {i}" if labels is None else f"row {i} ({labels[i]!r})"
        scale = np.max(np.abs(weights_row))
        if scale == 0:
            raise ValueError(f"SALC {row_name} has all-zero weights.")

        # Rescaling the whole row leaves the normalized SALC unchanged.
        row = weights_row / scale
        terms = {}
        for j in np.flatnonzero(row):
            weight = row[j] * norms[j]
            for gto in aos[j].gtos:
                primitive = (tuple(gto.p), gto.alpha, gto.l, gto.m, gto.n)
                terms.setdefault(primitive, []).append(weight * gto.c)

        # Merge identical primitives so exact cancellations are resolved
        # before computing the norm. Near cancellation in every primitive
        # indicates a row below floating-point resolution.
        primitive_weights = {primitive: fsum(values) for primitive, values in terms.items()}
        if not any(abs(value) > roundoff * fsum(map(abs, terms[primitive]))
                   for primitive, value in primitive_weights.items()):
            raise ValueError(f"SALC {row_name} cancels to zero.")

        nonzero = [(primitive, value) for primitive, value in primitive_weights.items() if value != 0]
        salc = CGF(list(nonzero[0][0][0]))
        for (position, alpha, l, m, n), value in nonzero:
            salc.add_gto_with_position(value, list(position), alpha, l, m, n)

        norm = salc.cgf._contraction_norm()
        if not isfinite(norm) or norm <= 0:
            raise ValueError(f"SALC {row_name} has zero or invalid norm.")
        cgfs.append(salc)
        transform[i] = norm * row

    return SALCBasis(cgfs, transform, labels)
