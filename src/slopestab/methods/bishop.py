"""Bishop's simplified method (1955).

Moment equilibrium about the slip-circle centre, with interslice forces assumed
horizontal.  The mobilised normal force is obtained from the classic
``m_alpha`` factor, which makes the equation implicit in ``F``; it is solved by
fixed-point iteration with under-relaxation, the way it is done in practice.
"""

from __future__ import annotations

import numpy as np

from ..exceptions import ConvergenceError
from ..geometry.slices import SliceSet
from .base import LEMResult, equilibrium_residuals
from .fellenius import fellenius

M_ALPHA_FLOOR = 0.02


def bishop(slices: SliceSet, max_iter: int = 200, tol: float = 1e-4) -> LEMResult:
    """Analyse a prepared slice set with Bishop's simplified method."""
    kh = np.abs(slices.kh_force)
    n_base = slices.w_eff - slices.u * slices.b - kh * slices.sin_a
    drive = slices.w_eff * slices.sin_a + kh * slices.cos_a
    total_drive = float(np.sum(drive))
    if total_drive <= 0.0:
        raise ArithmeticError(
            "bishop: net driving force is not positive - the surface is not "
            "a kinematically admissible failure mechanism"
        )

    try:
        fs = fellenius(slices).fs
    except ArithmeticError:
        fs = 1.3
    if not np.isfinite(fs) or fs <= 0.0:
        fs = 1.3

    clamped = 0
    converged = False
    it = 0
    for step in range(1, max_iter + 1):
        it = step
        m_alpha = slices.cos_a + slices.sin_a * slices.tan_phi / fs
        bad = m_alpha < M_ALPHA_FLOOR
        clamped += int(np.sum(bad))
        m_alpha = np.where(bad, M_ALPHA_FLOOR, m_alpha)

        numerator = (slices.c * slices.b + n_base * slices.tan_phi) / m_alpha
        fs_new = float(np.sum(numerator)) / total_drive

        if abs(fs_new - fs) < tol:
            fs = fs_new
            converged = True
            break
        fs = 0.5 * (fs + fs_new)
        if fs <= 0.0 or not np.isfinite(fs):
            raise ConvergenceError(
                "bishop",
                "iteration reached a non-physical F: pore pressures probably exceed "
                "the soil weight on some slices - check that the water table lies "
                "below the ground surface",
                fs,
            )

    if not converged:
        raise ConvergenceError(
            "bishop", f"no convergence within {max_iter} iterations (tol={tol:g})", fs
        )

    m_alpha = np.maximum(slices.cos_a + slices.sin_a * slices.tan_phi / fs, M_ALPHA_FLOOR)
    n_eff = (n_base - slices.c * slices.l / fs * slices.sin_a) / m_alpha
    shear = (slices.c * slices.l + n_eff * slices.tan_phi) / fs
    residuals = equilibrium_residuals(slices, n_eff, shear)

    info = {"m_alpha_floor_hits": clamped} if clamped else {}
    if clamped:
        info["note"] = "m_alpha was floored on steep slices; check the slip surface"
    return LEMResult(
        method="bishop",
        fs=fs,
        slices=slices,
        normal_forces=n_eff,
        shear_forces=shear,
        residuals=residuals,
        converged=converged,
        iterations=it,
        info=info,
    )
