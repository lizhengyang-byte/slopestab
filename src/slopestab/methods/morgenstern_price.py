"""Morgenstern-Price method (1965).

The general limit-equilibrium formulation: interslice forces follow
``X = lambda * f(x) * E`` for a prescribed shape ``f(x)``.  For each ``lambda``
the factor of safety is obtained from moment equilibrium, and ``lambda`` is
adjusted until the factor of safety from force equilibrium coincides with it.
This is the method most commercial software defaults to, and it reproduces
Spencer's result when ``f(x)`` is constant.

A particular ``f(x)`` shape may admit no ``lambda`` at all.  When that happens
the solver reports the closest value together with the force/moment gap rather
than pretending convergence.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..exceptions import ConvergenceError, GeometryError
from ..geometry.slices import SliceSet
from .base import LEMResult, check_factor_of_safety, equilibrium_residuals
from .gle import f_shape, force_imbalance, moment_imbalance, normal_forces, shear_forces

_FS_BRACKETS = ((0.2, 10.0), (0.05, 40.0), (0.01, 120.0))


def _bisect(fn, lo: float, hi: float, tol: float = 1e-6, max_iter: int = 100) -> float:
    """Bisection on a sign-changing bracket."""
    f_lo, f_hi = fn(lo), fn(hi)
    if not (np.isfinite(f_lo) and np.isfinite(f_hi)):
        raise ArithmeticError("non-finite bracket")
    if f_lo * f_hi > 0.0:
        raise ArithmeticError("bracket does not contain a root")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = fn(mid)
        if not np.isfinite(f_mid):
            raise ArithmeticError("non-finite value inside the bracket")
        if abs(f_mid) <= tol or 0.5 * (hi - lo) <= tol:
            return mid
        if f_lo * f_mid <= 0.0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


def _fs_from(fn, slices: SliceSet, f: np.ndarray, lam: float) -> float:
    """Factor of safety that zeroes ``fn`` for a fixed ``lambda``."""
    for lo, hi in _FS_BRACKETS:
        try:
            return _bisect(lambda fs: fn(slices, fs, lam, f), lo, hi, tol=1e-5)
        except ArithmeticError:
            continue
    raise ArithmeticError(f"no sign change in F for lambda={lam:.4f}")


def _false_position(delta, a: float, fa: float, b: float, fb: float, tol: float) -> float:
    """Secant/bisection hybrid that tolerates undefined trial values."""
    for _ in range(60):
        denom = fb - fa
        x = (a * fb - b * fa) / denom if abs(denom) > 1e-12 else 0.5 * (a + b)
        if not a < x < b:
            x = 0.5 * (a + b)
        fx = delta(x)
        if np.isnan(fx):
            x = 0.5 * (a + b) if abs(fa) < abs(fb) else 0.5 * (x + b)
            fx = delta(x)
            if np.isnan(fx):
                raise ConvergenceError(
                    "morgenstern_price", "lambda iteration left the solvable domain"
                )
        if abs(fx) <= tol or abs(b - a) <= tol:
            return x
        if fa * fx <= 0.0:
            b, fb = x, fx
        else:
            a, fa = x, fx
    return 0.5 * (a + b)


def morgenstern_price(
    slices: SliceSet,
    f: str | Sequence[float] | None = "half_sine",
    lam_range: tuple[float, float] = (-0.2, 1.2),
    n_scan: int = 31,
    tol: float = 1e-5,
) -> LEMResult:
    """Analyse a prepared slice set with the Morgenstern-Price method.

    Parameters
    ----------
    f:
        Interslice force function: ``"constant"``, ``"half_sine"``, ``"sine"``,
        ``"trapezoid"``, or an explicit array of length ``n_slices``.
    lam_range, n_scan:
        Bracket and resolution of the ``lambda`` scan.
    """
    if slices.surface_kind != "circle":
        raise GeometryError(
            "morgenstern_price: the moment equation is only defined for circular "
            "surfaces in v0.1; use bishop/janbu for a prescribed non-circular surface"
        )

    n = len(slices)
    if isinstance(f, str):
        f_values = f_shape(f, slices)
        f_name = f
    elif f is None:
        f_values = f_shape("half_sine", slices)
        f_name = "half_sine"
    else:
        f_values = np.asarray(f, dtype=float)
        if f_values.shape != (n,):
            raise ValueError(f"f must have length {n} (one value per slice)")
        f_name = "user"

    def delta(lam: float) -> float:
        try:
            fs_m = _fs_from(moment_imbalance, slices, f_values, lam)
            fs_f = _fs_from(force_imbalance, slices, f_values, lam)
        except ArithmeticError:
            return float("nan")
        return fs_m - fs_f

    lams = np.linspace(lam_range[0], lam_range[1], n_scan)
    finite = [(float(lam), float(delta(float(lam)))) for lam in lams]
    finite = [(lam, v) for lam, v in finite if np.isfinite(v)]
    if len(finite) < 2:
        raise ConvergenceError("morgenstern_price", "could not evaluate enough trial lambda values")

    pair = None
    for i in range(len(finite) - 1):
        if finite[i][1] * finite[i + 1][1] <= 0.0:
            pair = (finite[i], finite[i + 1])
            break

    note = None
    if pair is None:
        lam_best, v_best = min(finite, key=lambda p: abs(p[1]))
        if abs(v_best) > 0.02:
            lo_v = min(v for _, v in finite)
            hi_v = max(v for _, v in finite)
            raise ConvergenceError(
                "morgenstern_price",
                "no lambda satisfies both closures over "
                f"[{lam_range[0]:g}, {lam_range[1]:g}] (Fm-Ff in [{lo_v:.4f}, {hi_v:.4f}]); "
                "try another f(x) shape",
            )
        lam = lam_best
        note = f"no exact lambda bracket; accepted lambda={lam:.3f} with |Fm-Ff|={abs(v_best):.4f}"
    else:
        lam = _false_position(delta, pair[0][0], pair[0][1], pair[1][0], pair[1][1], tol)

    fs = _fs_from(moment_imbalance, slices, f_values, lam)
    n_force = normal_forces(slices, fs, lam, f_values)
    shear = shear_forces(slices, fs, n_force)
    residuals = equilibrium_residuals(slices, n_force, shear)

    info: dict[str, object] = {
        "lambda": float(lam),
        "f_function": f_name,
        "force_moment_gap": abs(float(delta(lam))),
    }
    if note:
        info["note"] = note
    negative = int(np.sum(n_force < 0.0))
    if negative:
        info["negative_normal_slices"] = negative

    return LEMResult(
        method="morgenstern_price",
        fs=check_factor_of_safety(fs, "morgenstern_price"),
        slices=slices,
        normal_forces=n_force,
        shear_forces=shear,
        residuals=residuals,
        converged=note is None,
        iterations=0,
        info=info,
    )
