"""Spencer's method (1967).

Assumes all interslice force resultants are parallel (``X/E = tan(theta)``,
i.e. ``f(x) = 1``) and solves two simultaneous equations - horizontal force
equilibrium and moment equilibrium - for the two unknowns ``(F, theta)``.
It is one of the classic "rigorous" checks: it satisfies both equilibrium
conditions exactly while using a single free parameter.
"""

from __future__ import annotations

import numpy as np

from ..exceptions import ConvergenceError, GeometryError
from ..geometry.slices import SliceSet
from .base import LEMResult, check_factor_of_safety, equilibrium_residuals
from .gle import constant_f, force_imbalance, moment_imbalance, normal_forces, shear_forces


def spencer(
    slices: SliceSet,
    fs0: float | None = None,
    lam0: float = 0.3,
    tol: float = 1e-6,
    max_iter: int = 80,
) -> LEMResult:
    """Analyse a prepared slice set with Spencer's method."""
    if slices.surface_kind != "circle":
        raise GeometryError(
            "spencer: moment equilibrium about a fixed centre is only defined for "
            "circular surfaces in v0.1; use bishop/janbu for a prescribed non-circular "
            "surface (rigorous non-circular GLE is on the v0.2 roadmap)"
        )
    n = len(slices)
    f = constant_f(n)
    scale = max(float(np.sum(slices.w)), 1e-9)
    span = max(float(slices.x_right[-1] - slices.x_left[0]), 1e-6)

    def residual(fs: float, lam: float) -> np.ndarray:
        return np.array(
            [
                force_imbalance(slices, fs, lam, f) / scale,
                moment_imbalance(slices, fs, lam, f) / (scale * span),
            ]
        )

    if fs0 is None:
        from .bishop import bishop

        try:
            fs0 = bishop(slices).fs
        except Exception:  # noqa: BLE001 - fall back to a neutral guess
            fs0 = 1.3
    x = np.array([float(fs0), float(lam0)])
    if not np.isfinite(x[0]) or x[0] <= 0.0:
        x[0] = 1.3

    converged = False
    it = 0
    for step in range(1, max_iter + 1):
        it = step
        r = residual(x[0], x[1])
        if np.max(np.abs(r)) < tol:
            converged = True
            break
        jac = np.empty((2, 2))
        dx_fs = 1e-4 * max(x[0], 1e-3)
        dx_lam = 1e-4
        r_fs = residual(x[0] + dx_fs, x[1])
        r_lam = residual(x[0], x[1] + dx_lam)
        jac[:, 0] = (r_fs - r) / dx_fs
        jac[:, 1] = (r_lam - r) / dx_lam
        try:
            step = np.linalg.solve(jac, -r)
        except np.linalg.LinAlgError as exc:
            raise ConvergenceError(
                "spencer", "singular Jacobian; the surface may be degenerate", x[0]
            ) from exc

        alpha = 1.0
        accepted = False
        for _ in range(8):
            trial = x + alpha * step
            if trial[0] <= 0.0 or not np.all(np.isfinite(trial)):
                alpha *= 0.5
                continue
            r_trial = residual(trial[0], trial[1])
            if np.max(np.abs(r_trial)) < np.max(np.abs(r)):
                x = trial
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            x = x + 0.05 * step
            if x[0] <= 0.0 or not np.all(np.isfinite(x)):
                raise ConvergenceError("spencer", "Newton iteration stagnated", x[0])

    if not converged:
        raise ConvergenceError(
            "spencer",
            f"no convergence within {max_iter} iterations (tol={tol:g})",
            float(x[0]),
        )

    fs, lam = float(x[0]), float(x[1])
    n_force = normal_forces(slices, fs, lam, f)
    shear = shear_forces(slices, fs, n_force)
    residuals = equilibrium_residuals(slices, n_force, shear)

    info: dict[str, object] = {"lambda": lam, "theta_deg": float(np.degrees(np.arctan(lam)))}
    negative = int(np.sum(n_force < 0.0))
    if negative:
        info["negative_normal_slices"] = negative
    return LEMResult(
        method="spencer",
        fs=check_factor_of_safety(fs, "spencer"),
        slices=slices,
        normal_forces=n_force,
        shear_forces=shear,
        residuals=residuals,
        converged=converged,
        iterations=it,
        info=info,
    )
