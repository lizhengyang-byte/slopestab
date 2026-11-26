"""Ordinary method of slices (Fellenius / Swedish circle method).

The oldest and simplest limit-equilibrium method: every slice is treated as an
independent block, interslice forces are ignored altogether.  It satisfies
moment equilibrium only and systematically under-estimates the factor of
safety (typically 5-15% low, more for high pore pressures), which makes it
useful as a conservative screening bound rather than a design value.
"""

from __future__ import annotations

import numpy as np

from ..geometry.slices import SliceSet
from .base import LEMResult, check_factor_of_safety, equilibrium_residuals


def fellenius(slices: SliceSet) -> LEMResult:
    """Analyse a prepared slice set with the ordinary method of slices."""
    kh = np.abs(slices.kh_force)
    n_eff = slices.w_eff * slices.cos_a - kh * slices.sin_a - slices.u * slices.l
    resist = slices.c * slices.l + n_eff * slices.tan_phi
    drive = slices.w_eff * slices.sin_a + kh * slices.cos_a

    total_drive = float(np.sum(drive))
    if total_drive <= 0.0:
        raise ArithmeticError(
            "fellenius: net driving force is not positive - the surface is not "
            "a kinematically admissible failure mechanism"
        )

    fs = check_factor_of_safety(float(np.sum(resist)) / total_drive, "fellenius")
    shear = resist / fs
    residuals = equilibrium_residuals(slices, n_eff, shear)

    return LEMResult(
        method="fellenius",
        fs=fs,
        slices=slices,
        normal_forces=n_eff,
        shear_forces=shear,
        residuals=residuals,
        converged=True,
        iterations=0,
        info={"note": "interslice forces neglected; conservative bound"},
    )
