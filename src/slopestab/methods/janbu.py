"""Janbu's simplified method (1954, corrected 1973).

Horizontal force equilibrium with no interslice shear.  The correction factor
``f0`` accounts for the fact that the simplified assumption under-estimates the
factor of safety; it depends on the slip-surface shape and on the strength
mobilised, and in practice is read from Janbu's chart.

This implementation keeps ``f0`` explicit.  The default is ``f0 = 1.0`` (the
uncorrected method), and a chart-read value can be supplied by the user.  The
uncorrected result is typically 5-20% low, so the package recommends
Morgenstern-Price or Spencer for final design values.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import numpy as np

from ..exceptions import ConvergenceError
from ..geometry.slices import SliceSet
from .base import LEMResult, equilibrium_residuals


def janbu(
    slices: SliceSet,
    f0: float = 1.0,
    max_iter: int = 200,
    tol: float = 1e-4,
) -> LEMResult:
    """Analyse a prepared slice set with Janbu's simplified method.

    Parameters
    ----------
    slices:
        Slice set built by :func:`slopestab.geometry.slices.build_slices`.
    f0:
        Janbu correction factor read from the 1973 chart.  ``1.0`` reproduces
        the uncorrected simplified method.
    """
    if f0 <= 0.0:
        raise ValueError("f0 must be positive")

    kh = np.abs(slices.kh_force)
    n_base = slices.w_eff - slices.u * slices.b - kh * slices.sin_a
    tan_a = slices.sin_a / slices.cos_a
    denom = float(np.sum(slices.w_eff * tan_a + kh))
    if denom <= 0.0:
        raise ArithmeticError(
            "janbu: net driving force is not positive - the surface is not "
            "a kinematically admissible failure mechanism"
        )

    fs = 1.3
    converged = False
    it = 0
    for step in range(1, max_iter + 1):
        it = step
        n_alpha = slices.cos_a**2 * (1.0 + tan_a * slices.tan_phi / fs)
        n_alpha = np.maximum(n_alpha, 1e-6)
        numerator = (slices.c * slices.b + n_base * slices.tan_phi) / n_alpha
        fs_new = f0 * float(np.sum(numerator)) / denom
        if abs(fs_new - fs) < tol:
            fs = fs_new
            converged = True
            break
        fs = 0.5 * (fs + fs_new)
        if fs <= 0.0 or not np.isfinite(fs):
            raise ConvergenceError("janbu", "iteration diverged to a non-physical F", fs)

    if not converged:
        raise ConvergenceError(
            "janbu", f"no convergence within {max_iter} iterations (tol={tol:g})", fs
        )

    n_alpha = np.maximum(slices.cos_a**2 * (1.0 + tan_a * slices.tan_phi / fs), 1e-6)
    n_eff = n_base / n_alpha
    term = (slices.c * slices.b + n_base * slices.tan_phi) / n_alpha
    shear = term / fs
    residuals = equilibrium_residuals(slices, n_eff, shear)

    info: dict[str, object] = {"f0": f0}
    if abs(f0 - 1.0) < 1e-12:
        info["note"] = (
            "uncorrected Janbu simplified: typically 5-20% below Spencer / "
            "Morgenstern-Price; supply f0 for the corrected value"
        )
    return LEMResult(
        method="janbu",
        fs=fs,
        slices=slices,
        normal_forces=n_base,
        shear_forces=shear,
        residuals=residuals,
        converged=converged,
        iterations=it,
        info=info,
    )
