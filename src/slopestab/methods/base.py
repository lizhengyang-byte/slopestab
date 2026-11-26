"""Shared limit-equilibrium machinery: results and equilibrium auditing."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ..geometry.slices import SliceSet


@dataclass
class LEMResult:
    """Outcome of one limit-equilibrium analysis.

    Attributes
    ----------
    method:
        Method name (``"bishop"``, ``"spencer"``, ...).
    fs:
        Factor of safety.
    slices:
        The :class:`~slopestab.geometry.slices.SliceSet` that was solved.
    normal_forces, shear_forces:
        Per-slice effective normal force ``N`` and mobilised shear force ``S``
        (kN/m of slope width).
    residuals:
        Global force and moment imbalance, expressed as a percentage of
        ``sum(W)``.  Small values mean the method's equilibrium assumptions are
        satisfied; Bishop/Fellenius are only in moment equilibrium, so their
        force residuals are expected to be non-zero.
    """

    method: str
    fs: float
    slices: SliceSet
    normal_forces: np.ndarray
    shear_forces: np.ndarray
    residuals: dict[str, float] = field(default_factory=dict)
    converged: bool = True
    iterations: int = 0
    info: dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------ summary
    @property
    def resisting_force(self) -> float:
        return float(np.sum(self.shear_forces))

    @property
    def driving_force(self) -> float:
        return float(self.resisting_force / self.fs) if self.fs else float("nan")

    def summary(self) -> dict[str, Any]:
        """Compact, JSON-friendly summary of the analysis."""
        out: dict[str, Any] = {
            "method": self.method,
            "fs": round(float(self.fs), 6),
            "n_slices": len(self.slices),
            "total_weight_kN_per_m": round(self.slices.total_weight, 3),
            "residual_percent": {k: round(v, 4) for k, v in self.residuals.items()},
            "converged": self.converged,
            "iterations": self.iterations,
        }
        out.update(self.info)
        if self.slices.notes:
            out["notes"] = list(self.slices.notes)
        return out

    def to_dict(self) -> dict[str, Any]:
        return {"summary": self.summary(), "slices": self.slices.as_table()}

    def to_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
        return path

    def to_csv(self, path: str | Path) -> Path:
        rows = self.slices.as_table()
        path = Path(path)
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        return path

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"LEMResult({self.method}, Fs={self.fs:.4f}, slices={len(self.slices)})"


def equilibrium_residuals(
    slices: SliceSet, n_force: np.ndarray, s_force: np.ndarray
) -> dict[str, float]:
    """Global force/moment imbalance of a solution, as percentages of ``sum(W)``.

    ``Fx`` and ``Fy`` are computed with the base reaction ``N`` along the inward
    normal and the mobilised shear ``S`` opposed to the sliding direction.
    """
    n_force = np.asarray(n_force, dtype=float)
    s_force = np.asarray(s_force, dtype=float)
    n_total = n_force + slices.u * slices.l  # effective + pore-water resultant

    fx = float(np.sum(n_total * slices.nx - s_force * slices.ux) + np.sum(slices.kh_force))
    fy = float(np.sum(n_total * slices.ny - s_force * slices.uy) - np.sum(slices.w_eff))

    base_dy = n_total * slices.ny - s_force * slices.uy
    base_dx = n_total * slices.nx - s_force * slices.ux
    m_base = float(
        np.sum((slices.x_mid - slices.x_ref) * base_dy - (slices.y_base - slices.y_ref) * base_dx)
    )
    m_load = float(
        np.sum(
            (slices.x_mid - slices.x_ref) * (-slices.w_eff)
            - (slices.y_c - slices.y_ref) * slices.kh_force
        )
    )

    normaliser = max(float(np.sum(slices.w)), 1e-9)
    return {
        "Fx": 100.0 * fx / normaliser,
        "Fy": 100.0 * fy / normaliser,
        "M": 100.0 * (m_base + m_load) / normaliser,
        "Fx_abs": abs(fx),
        "Fy_abs": abs(fy),
        "M_abs": abs(m_base + m_load),
    }


def check_factor_of_safety(fs: float, method: str) -> float:
    """Reject non-physical factors of safety before they reach a report."""
    if not np.isfinite(fs):
        raise ArithmeticError(f"{method}: factor of safety is not finite")
    if fs <= 0.0:
        raise ArithmeticError(
            f"{method}: non-positive factor of safety ({fs:.4f}); "
            "the surface may not be a kinematically admissible failure mechanism"
        )
    return float(fs)
