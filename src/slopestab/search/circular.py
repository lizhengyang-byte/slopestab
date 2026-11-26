"""Circular slip-surface search: grid scan plus local refinement."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

from ..analysis import run_on_slices
from ..geometry.profile import Slope
from ..geometry.slices import build_slices
from ..geometry.soil import SoilModel
from ..geometry.surface import SlipCircle
from ..loads import Seismic, SurchargeSet, WaterCondition

Candidate = tuple[float, float, float, float]  # cx, cy, r, fs


@dataclass
class SearchGrid:
    """Ranges and resolutions of the three-dimensional centre/radius grid."""

    cx_range: tuple[float, float]
    cy_range: tuple[float, float]
    r_range: tuple[float, float]
    n_cx: int = 24
    n_cy: int = 14
    n_r: int = 12

    def __post_init__(self) -> None:
        for name in ("cx_range", "cy_range", "r_range"):
            lo, hi = getattr(self, name)
            if hi <= lo:
                raise ValueError(f"{name} must be increasing")
        if min(self.n_cx, self.n_cy, self.n_r) < 3:
            raise ValueError("grid resolutions must be at least 3")

    @classmethod
    def default(cls, slope: Slope) -> SearchGrid:
        """A sensible grid derived from the slope height and extent."""
        h = max(slope.height, 1e-3)
        return cls(
            cx_range=(slope.x_min - 0.1 * h, slope.x_min + 2.6 * h),
            cy_range=(slope.y_max + 0.05 * h, slope.y_max + 2.8 * h),
            r_range=(0.5 * h, 3.2 * h),
        )

    @property
    def size(self) -> int:
        return self.n_cx * self.n_cy * self.n_r


@dataclass
class SearchResult:
    """Outcome of a slip-surface search."""

    critical: Any  # LEMResult
    candidates: list[Candidate]
    grid: SearchGrid
    method: str
    elapsed: float = 0.0
    n_evaluated: int = 0
    slope: Any = None
    soil: Any = None
    water: Any = None
    surcharges: Any = None
    seismic: Any = None
    notes: list[str] = field(default_factory=list)

    @property
    def fs(self) -> float:
        return float(self.critical.fs)

    def summary(self) -> dict[str, Any]:
        circle = getattr(self.critical, "surface", None)
        out: dict[str, Any] = {
            "method": self.method,
            "critical_fs": round(self.fs, 6),
            "n_evaluated": self.n_evaluated,
            "elapsed_s": round(self.elapsed, 3),
        }
        if isinstance(circle, SlipCircle):
            out["critical_circle"] = {
                "cx": round(circle.cx, 4),
                "cy": round(circle.cy, 4),
                "r": round(circle.r, 4),
            }
        if self.notes:
            out["notes"] = list(self.notes)
        return out

    def to_dataframe(self):
        try:
            import pandas as pd
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError("pandas is required for to_dataframe()") from exc
        return pd.DataFrame(
            [{"cx": c[0], "cy": c[1], "r": c[2], "fs": c[3]} for c in self.candidates]
        )


def search_circular(
    slope: Slope,
    soil: SoilModel,
    method: str = "bishop",
    grid: SearchGrid | None = None,
    n_slices: int = 20,
    refine: bool = True,
    refine_candidates: int = 3,
    refine_n_slices: int = 30,
    max_refine_iter: int = 60,
    water: WaterCondition | None = None,
    surcharges: SurchargeSet | None = None,
    seismic: Seismic | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> SearchResult:
    """Search for the critical circular slip surface.

    The scan covers a three-dimensional grid of centres ``(cx, cy)`` and radii
    ``r``.  Candidates are filtered for geometric admissibility (two ground
    intersections, non-empty sliding mass) before any mechanics is computed,
    and the best few are then polished with a Nelder-Mead simplex.

    Parameters
    ----------
    method:
        Limit-equilibrium method used as the objective.
    grid:
        Search domain; defaults to :meth:`SearchGrid.default`.
    n_slices:
        Slices used during the coarse scan (refinement uses more).
    refine:
        Polish the best candidates locally.

    Returns
    -------
    SearchResult
        ``.critical`` is the full :class:`~slopestab.methods.base.LEMResult` of
        the governing surface, ``.candidates`` the ranked scan results.
    """
    grid = grid or SearchGrid.default(slope)

    t0 = time.perf_counter()
    cx_values = np.linspace(*grid.cx_range, grid.n_cx)
    cy_values = np.linspace(*grid.cy_range, grid.n_cy)
    r_values = np.linspace(*grid.r_range, grid.n_r)

    candidates: list[Candidate] = []
    evaluated = 0
    skipped = 0
    total = grid.size
    counter = 0

    for cx in cx_values:
        for cy in cy_values:
            for r in r_values:
                counter += 1
                if progress is not None and counter % max(1, total // 20) == 0:
                    progress(counter, total)
                circle = SlipCircle(float(cx), float(cy), float(r))
                if not _admissible(slope, circle):
                    skipped += 1
                    continue
                try:
                    slices = build_slices(
                        slope,
                        soil,
                        circle,
                        water=water,
                        surcharges=surcharges,
                        seismic=seismic,
                        n_slices=n_slices,
                    )
                    result = run_on_slices(method, slices)
                except Exception:  # noqa: BLE001 - invalid trial surfaces are expected
                    skipped += 1
                    continue
                evaluated += 1
                candidates.append((float(cx), float(cy), float(r), float(result.fs)))

    if not candidates:
        raise ValueError(
            "no admissible circular slip surface was found inside the search grid; "
            "widen the ranges or check the geometry"
        )

    candidates.sort(key=lambda c: c[3])
    notes = [f"coarse scan: {evaluated} surfaces analysed, {skipped} rejected as inadmissible"]

    best = candidates[0]
    critical: Any = None
    if refine:
        best, critical = _refine(
            best,
            slope,
            soil,
            method,
            water,
            surcharges,
            seismic,
            refine_candidates,
            candidates,
            refine_n_slices,
            max_refine_iter,
        )
        notes.append(
            f"refined {min(refine_candidates, len(candidates))} candidate(s) with Nelder-Mead"
        )
    else:
        critical = _analyse_circle(
            slope, soil, (best[0], best[1], best[2]), method, water, surcharges, seismic, n_slices
        )

    elapsed = time.perf_counter() - t0
    return SearchResult(
        critical=critical,
        candidates=candidates,
        grid=grid,
        method=method,
        elapsed=elapsed,
        n_evaluated=evaluated,
        slope=slope,
        soil=soil,
        water=water,
        surcharges=surcharges,
        seismic=seismic,
        notes=notes,
    )


def _admissible(slope: Slope, circle: SlipCircle) -> bool:
    """Cheap geometric filter applied before any slice work is done."""
    if circle.r <= 0.0:
        return False
    if not (slope.x_min <= circle.cx <= slope.x_max):
        return False
    lowest = circle.cy - circle.r
    if lowest >= slope.y_max:  # entirely above the ground
        return False
    if circle.base_y(circle.cx) >= slope.y_at(circle.cx):  # no soil above the arc
        return False
    ints = circle.intersections(slope)
    if ints.size < 2:
        return False
    span = float(ints[-1] - ints[0])
    return span > 0.05 * max(slope.height, 1e-3)


def _analyse_circle(
    slope: Slope,
    soil: SoilModel,
    circle: tuple[float, float, float],
    method: str,
    water: WaterCondition | None,
    surcharges: SurchargeSet | None,
    seismic: Seismic | None,
    n_slices: int,
):
    slices = build_slices(
        slope,
        soil,
        SlipCircle(*circle),
        water=water,
        surcharges=surcharges,
        seismic=seismic,
        n_slices=n_slices,
    )
    result = run_on_slices(method, slices)
    result.slope = slope
    result.soil = soil
    result.surface = SlipCircle(*circle)
    return result


def _objective(
    params: np.ndarray,
    slope: Slope,
    soil: SoilModel,
    method: str,
    water: WaterCondition | None,
    surcharges: SurchargeSet | None,
    seismic: Seismic | None,
    n_slices: int,
) -> float:
    cx, cy, r = float(params[0]), float(params[1]), float(params[2])
    circle = SlipCircle(cx, cy, r)
    if not _admissible(slope, circle):
        return float("inf")
    try:
        return float(
            _analyse_circle(
                slope, soil, (cx, cy, r), method, water, surcharges, seismic, n_slices
            ).fs
        )
    except Exception:  # noqa: BLE001
        return float("inf")


def _refine(
    best: Candidate,
    slope: Slope,
    soil: SoilModel,
    method: str,
    water: WaterCondition | None,
    surcharges: SurchargeSet | None,
    seismic: Seismic | None,
    refine_candidates: int,
    candidates: list[Candidate],
    n_slices: int,
    max_iter: int,
) -> tuple[Candidate, Any]:
    h = max(slope.height, 1e-3)
    starts = candidates[: max(1, refine_candidates)]
    best_params = np.array([best[0], best[1], best[2]], dtype=float)
    best_fs = best[3]

    for start in starts:
        x0 = np.array([start[0], start[1], start[2]], dtype=float)
        steps = np.array([0.15 * h, 0.15 * h, 0.1 * h], dtype=float)

        def fun(p: np.ndarray) -> float:
            return _objective(p, slope, soil, method, water, surcharges, seismic, n_slices)

        p, fv = _nelder_mead(fun, x0, steps, max_iter=max_iter, tol=1e-6)
        if fv < best_fs:
            best_fs = float(fv)
            best_params = np.asarray(p, dtype=float)

    critical = _analyse_circle(
        slope,
        soil,
        (float(best_params[0]), float(best_params[1]), float(best_params[2])),
        method,
        water,
        surcharges,
        seismic,
        n_slices,
    )
    return (float(best_params[0]), float(best_params[1]), float(best_params[2]), best_fs), critical


def _nelder_mead(
    fun: Callable[[np.ndarray], float],
    x0: np.ndarray,
    steps: np.ndarray,
    max_iter: int = 80,
    tol: float = 1e-6,
) -> tuple[np.ndarray, float]:
    """Minimal Nelder-Mead simplex (keeps the core free of SciPy)."""
    n = x0.size
    simplex = [np.asarray(x0, dtype=float)]
    for i in range(n):
        p = np.array(x0, dtype=float)
        p[i] += steps[i]
        simplex.append(p)
    values = [float(fun(p)) for p in simplex]

    for _ in range(max_iter):
        order = np.argsort(values)
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        if not np.all(np.isfinite(values)):
            break
        if abs(values[-1] - values[0]) <= tol:
            break

        centroid = np.mean(simplex[:-1], axis=0)
        x_worst = simplex[-1]
        x_r = centroid + (centroid - x_worst)
        f_r = float(fun(x_r))

        if f_r < values[0]:
            x_e = centroid + 2.0 * (centroid - x_worst)
            f_e = float(fun(x_e))
            if f_e < f_r:
                simplex[-1], values[-1] = x_e, f_e
            else:
                simplex[-1], values[-1] = x_r, f_r
        elif f_r < values[-2]:
            simplex[-1], values[-1] = x_r, f_r
        else:
            x_c = centroid + 0.5 * (x_worst - centroid)
            f_c = float(fun(x_c))
            if f_c < values[-1]:
                simplex[-1], values[-1] = x_c, f_c
            else:
                for i in range(1, n + 1):
                    simplex[i] = simplex[0] + 0.5 * (simplex[i] - simplex[0])
                    values[i] = float(fun(simplex[i]))

    i_best = int(np.argmin(values))
    return simplex[i_best], values[i_best]
