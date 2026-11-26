"""Shared equations for the rigorous methods (Spencer and Morgenstern-Price).

Both methods assume an interslice force resultant whose inclination is given by
``X = lambda * f(x) * E``.  With that assumption each slice can be put in force
equilibrium independently, which yields an explicit expression for the
effective normal force.  The two methods then differ only in how ``lambda`` is
determined:

* Spencer: ``f(x) = 1`` and ``(F, lambda)`` are solved simultaneously from
  horizontal force and moment equilibrium;
* Morgenstern-Price: ``f(x)`` is a prescribed shape and ``lambda`` is found so
  that the factor of safety from force equilibrium equals the one from moment
  equilibrium.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from ..geometry.slices import SliceSet

FShape = Callable[[np.ndarray], np.ndarray]


def normal_forces(slices: SliceSet, fs: float, lam: float, f: np.ndarray) -> np.ndarray:
    """Effective normal forces for a given ``F`` and interslice ratio ``lambda``.

    Derived from per-slice vertical equilibrium with ``dX = lambda*f*dE``.
    """
    sd = slices.slide_dir
    tp = slices.tan_phi
    kh = np.abs(slices.kh_force)
    cohesion = slices.c * slices.l / fs

    # Pore water acts along the inward normal: it reduces the effective vertical
    # load by u*b and adds u*l*sin_a to the tangential balance.
    u_l = slices.u * slices.l
    num = (
        slices.w_eff
        - slices.u * slices.b
        + lam * f * sd * kh
        + lam * f * sd * u_l * slices.sin_a
        + cohesion * (-slices.sin_a - lam * f * sd * slices.cos_a)
    )
    den = (
        slices.cos_a
        - lam * f * sd * slices.sin_a
        + (tp / fs) * (slices.sin_a + lam * f * sd * slices.cos_a)
    )
    if np.any(np.abs(den) < 1e-9):
        raise ArithmeticError(
            "the slice force equations are singular (denominators near zero); "
            "check the slip surface geometry"
        )
    return num / den


def shear_forces(slices: SliceSet, fs: float, n_force: np.ndarray) -> np.ndarray:
    """Mobilised shear forces ``S = (c'l + N tan phi') / F``."""
    return (slices.c * slices.l + n_force * slices.tan_phi) / fs


def force_imbalance(slices: SliceSet, fs: float, lam: float, f: np.ndarray) -> float:
    """Net horizontal force imbalance (kN/m); zero when force equilibrium holds."""
    n = normal_forces(slices, fs, lam, f)
    s = shear_forces(slices, fs, n)
    n_total = n + slices.u * slices.l
    de = -np.abs(slices.kh_force) - n_total * slices.sin_a + s * slices.cos_a
    return float(np.sum(de))


def moment_imbalance(slices: SliceSet, fs: float, lam: float, f: np.ndarray) -> float:
    """Net moment about the reference point (kN*m/m); zero at moment equilibrium."""
    n = normal_forces(slices, fs, lam, f)
    s = shear_forces(slices, fs, n)
    n_total = n + slices.u * slices.l
    base_dy = n_total * slices.ny - s * slices.uy
    base_dx = n_total * slices.nx - s * slices.ux
    m_base = np.sum(
        (slices.x_mid - slices.x_ref) * base_dy - (slices.y_base - slices.y_ref) * base_dx
    )
    m_load = np.sum(
        (slices.x_mid - slices.x_ref) * (-slices.w_eff)
        - (slices.y_c - slices.y_ref) * slices.kh_force
    )
    return float(m_base + m_load)


def constant_f(n: int) -> np.ndarray:
    return np.ones(n, dtype=float)


def f_shape(kind: str, slices: SliceSet) -> np.ndarray:
    """Build the interslice force function ``f(x)`` from a named shape.

    Supported shapes: ``constant``, ``half_sine``, ``sine``, ``trapezoid``.
    """
    x = np.asarray(slices.x_mid, dtype=float)
    x0, x1 = float(slices.x_left[0]), float(slices.x_right[-1])
    t = (x - x0) / max(x1 - x0, 1e-12)
    kind = kind.lower().replace("-", "_").replace(" ", "_")
    if kind in ("constant", "const", "spencer"):
        return np.ones_like(t)
    if kind in ("half_sine", "halfsine", "sine_half"):
        return np.sin(np.pi * t)
    if kind in ("sine", "full_sine"):
        return np.sin(2.0 * np.pi * t) * -1.0
    if kind in ("trapezoid", "trapezoidal"):
        taper = np.minimum(1.0, np.minimum(4.0 * t, 4.0 * (1.0 - t)))
        return np.clip(taper, 0.05, 1.0)
    raise ValueError(
        f"unknown interslice force function {kind!r}; "
        "use 'constant', 'half_sine', 'sine', 'trapezoid' or pass an array"
    )


def bisect(
    fn: Callable[[float], float], lo: float, hi: float, tol: float = 1e-6, max_iter: int = 120
) -> float:
    """Robust bisection on a sign-changing bracket (no SciPy dependency)."""
    f_lo, f_hi = fn(lo), fn(hi)
    if not (np.isfinite(f_lo) and np.isfinite(f_hi)):
        raise ArithmeticError("bisection: non-finite bracket values")
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    if f_lo * f_hi > 0.0:
        raise ArithmeticError("bisection: bracket does not contain a root")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = fn(mid)
        if not np.isfinite(f_mid):
            raise ArithmeticError("bisection: non-finite value inside the bracket")
        if abs(f_mid) <= tol or 0.5 * (hi - lo) <= tol:
            return mid
        if f_lo * f_mid <= 0.0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)
