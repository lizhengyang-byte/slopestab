"""Slice discretisation.

The slice builder is the part that decides whether the numbers are trustworthy,
so it works a little harder than the textbook formulas:

* cuts are inserted wherever the slip surface crosses a layer boundary, the
  water table or a surcharge edge, so every slice base lies inside exactly one
  soil layer with a single pore-pressure rule;
* the slice weight is integrated piecewise with the correct unit weight above
  and below the water table, instead of averaging the whole slice;
* the base inclination is expressed through the direction cosines
  ``(sin_a, cos_a)`` derived from the actual geometry, which removes the
  sign ambiguity between "down-slope" and "up-slope" slices.

Conventions used throughout the package
---------------------------------------
``cos_a > 0`` always (the base carries compression);
``sin_a > 0`` means the base descends in the direction the mass is moving, so
``W * sin_a`` is a driving component.  ``(u_x, u_y)`` is the unit tangent
pointing in the direction of movement and ``(n_x, n_y)`` the unit normal
pointing from the base into the sliding mass.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from ..loads import NoWater, Seismic, SurchargeSet, WaterCondition
from .profile import Slope
from .soil import SoilModel
from .surface import SlipSurface


@dataclass
class SliceSet:
    """Vectorised slice geometry and loads for one slip surface."""

    x_left: np.ndarray
    x_right: np.ndarray
    b: np.ndarray
    x_mid: np.ndarray
    y_base: np.ndarray
    y_ground: np.ndarray
    h: np.ndarray
    l: np.ndarray
    c: np.ndarray
    phi: np.ndarray
    tan_phi: np.ndarray
    w: np.ndarray
    w_eff: np.ndarray
    kh_force: np.ndarray
    u: np.ndarray
    sin_a: np.ndarray
    cos_a: np.ndarray
    ux: np.ndarray
    uy: np.ndarray
    nx: np.ndarray
    ny: np.ndarray
    y_c: np.ndarray
    area: np.ndarray
    slide_dir: float
    x_ref: float
    y_ref: float
    layer_index: np.ndarray
    water_table: np.ndarray | None = None
    notes: list[str] | None = None
    surface_kind: str = "circle"

    # ------------------------------------------------------------------ helpers
    def __len__(self) -> int:
        return int(self.b.size)

    @property
    def total_weight(self) -> float:
        return float(self.w.sum())

    @property
    def span(self) -> tuple[float, float]:
        return (float(self.x_left[0]), float(self.x_right[-1]))

    def as_table(self) -> list[dict[str, Any]]:
        """Slice data as a list of plain dictionaries (CSV/JSON friendly)."""
        rows = []
        for i in range(len(self)):
            rows.append(
                {
                    "slice": i + 1,
                    "x_left": round(float(self.x_left[i]), 6),
                    "x_right": round(float(self.x_right[i]), 6),
                    "b": round(float(self.b[i]), 6),
                    "x_mid": round(float(self.x_mid[i]), 6),
                    "y_base": round(float(self.y_base[i]), 6),
                    "h": round(float(self.h[i]), 6),
                    "l": round(float(self.l[i]), 6),
                    "w": round(float(self.w[i]), 6),
                    "u": round(float(self.u[i]), 6),
                    "sin_a": round(float(self.sin_a[i]), 8),
                    "cos_a": round(float(self.cos_a[i]), 8),
                    "c": round(float(self.c[i]), 6),
                    "phi": round(float(self.phi[i]), 6),
                }
            )
        return rows

    def to_dataframe(self):
        """Slice data as a :class:`pandas.DataFrame` (pandas is optional)."""
        try:
            import pandas as pd
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError(
                "pandas is required for to_dataframe(); use .as_table() instead"
            ) from exc
        return pd.DataFrame(self.as_table())


def build_slices(
    slope: Slope,
    soil: SoilModel,
    surface: SlipSurface,
    water: WaterCondition | None = None,
    surcharges: SurchargeSet | None = None,
    seismic: Seismic | None = None,
    n_slices: int = 30,
) -> SliceSet:
    """Discretise the sliding mass into vertical slices.

    Parameters
    ----------
    slope, soil, surface:
        Geometry, materials and the slip surface to analyse.
    water, surcharges, seismic:
        Optional loads; defaults to dry, unloaded, static conditions.
    n_slices:
        Number of equal-width slices *before* extra cuts at layer
        boundaries, the water table and surcharge edges are inserted.
        The final slice count is therefore >= ``n_slices``.

    Raises
    ------
    ValueError
        If the surface does not intersect the ground surface twice, if the
        resulting sliding mass is empty, or if the pore pressure anywhere
        exceeds the slice weight (water table drawn above the ground).
    """
    water = water or NoWater()
    surcharges = surcharges or SurchargeSet()
    seismic = seismic or Seismic()
    if n_slices < 5:
        raise ValueError("n_slices should be at least 5")

    ints = surface.intersections(slope)
    if ints.size < 2:
        raise ValueError(
            "the slip surface must cross the ground surface at least twice "
            f"(found {ints.size} intersection(s))"
        )
    x_l, x_r = float(ints[0]), float(ints[-1])
    if x_r - x_l < 1e-6:
        raise ValueError("the slip surface intersection span is degenerate")

    # ------------------------------------------------------------------ cut points
    cuts = [x_l, x_r]
    cuts.extend(np.linspace(x_l, x_r, n_slices + 1).tolist())
    cuts.extend(float(v) for v in slope.vertices_between(x_l, x_r))
    cuts.extend(_layer_crossings(surface, soil, x_l, x_r))
    cuts.extend(float(v) for v in water.crossing_x(surface) if x_l < float(v) < x_r)
    cuts.extend(x for x in surcharges.edges() if x_l < x < x_r)
    xs = _unique_sorted(cuts, tol=1e-6 * max(1.0, x_r - x_l))

    # --------------------------------------------------------------- slice geometry
    y_ground_ref = float(slope.y_at(x_l)), float(slope.y_at(x_r))
    slide_dir = -1.0 if y_ground_ref[1] >= y_ground_ref[0] else 1.0

    x_left = xs[:-1]
    x_right = xs[1:]
    b = x_right - x_left
    x_mid = 0.5 * (x_left + x_right)
    y_base = np.asarray(surface.base_y(x_mid), dtype=float)
    y_ground = np.asarray(slope.y_at(x_mid), dtype=float)
    h = y_ground - y_base

    keep = (b > 1e-6) & (h > 1e-4)
    x_left, x_right, b, x_mid = x_left[keep], x_right[keep], b[keep], x_mid[keep]
    y_base, y_ground, h = y_base[keep], y_ground[keep], h[keep]
    if b.size == 0:
        raise ValueError("the slip surface encloses no soil above it")

    slope_ratio = np.asarray(surface.base_slope(x_mid), dtype=float)
    cos_a = 1.0 / np.sqrt(1.0 + slope_ratio**2)
    sin_a = -slide_dir * slope_ratio * cos_a
    ux = slide_dir * cos_a
    uy = -sin_a
    nx = slide_dir * sin_a
    ny = cos_a

    l = np.array([surface.length(a, c) for a, c in zip(x_left, x_right)], dtype=float)
    l = np.maximum(l, b * cos_a * 0.5 + 1e-9)

    # ------------------------------------------------------------------ slice loads
    n = b.size
    w = np.zeros(n)
    moment_y = np.zeros(n)
    area = np.zeros(n)
    u_press = np.zeros(n)
    layer_idx = np.zeros(n, dtype=int)
    water_table = np.full(n, np.nan)

    for i in range(n):
        xm, yb, yg = float(x_mid[i]), float(y_base[i]), float(y_ground[i])
        wy = water.water_y(xm)
        wy = None if wy is None else float(np.asarray(wy))
        if wy is not None:
            water_table[i] = wy

        ys = [yb]
        if wy is not None and yb < wy < yg:
            ys.append(wy)
        ys.extend(soil.boundaries_between(xm, yb, yg))
        ys.append(yg)
        ys = _unique_sorted(ys, tol=1e-9)

        layer_idx[i] = soil.layer_index(xm, yb)
        for ya, yb2 in zip(ys[:-1], ys[1:]):
            if yb2 <= ya:
                continue
            mid = 0.5 * (ya + yb2)
            layer = soil.layer_at(xm, mid)
            below = wy is not None and mid < wy
            gamma = layer.unit_weight(below_water=below)
            thick = yb2 - ya
            a_sub = b[i] * thick
            w[i] += gamma * a_sub
            moment_y[i] += gamma * a_sub * mid
            area[i] += a_sub

        q = surcharges.q_at(xm)
        if q > 0.0:
            w_q = q * b[i]
            w[i] += w_q
            moment_y[i] += w_q * yg

        u_press[i] = water.pressure(xm, yb, float(w[i]), float(b[i]))

    c_arr = np.zeros(n)
    phi_arr = np.zeros(n)
    for i in range(n):
        layer = soil.layer_at(float(x_mid[i]), float(y_base[i]))
        c_arr[i] = layer.c
        phi_arr[i] = layer.phi

    tan_phi = np.tan(np.radians(phi_arr))
    w_eff = w * (1.0 + seismic.kv)
    kh_force = seismic.kh * w * slide_dir

    # A water table above the ground surface drives the pore pressure past the
    # slice weight, which makes the effective normal force negative and turns the
    # factor of safety into a meaningless number.  Refuse instead of reporting it.
    if np.any(u_press * b > w + 1e-9):
        worst = int(np.argmax(u_press * b - w))
        raise ValueError(
            "pore pressure exceeds the slice weight at x="
            f"{float(x_mid[worst]):.2f} m (u*b={float(u_press[worst] * b[worst]):.2f} kN/m vs "
            f"W={float(w[worst]):.2f} kN/m): the water table is above the ground surface. "
            "Lower the phreatic line to or below the ground surface, or model the artesian "
            "pressure explicitly with r_u / a pressure grid."
        )

    y_c = np.where(w > 0.0, moment_y / np.where(w > 0, w, 1.0), 0.5 * (y_base + y_ground))

    x_ref, y_ref = surface.moment_center()

    notes: list[str] = []
    if np.any(l <= 0.0):
        notes.append("zero-length slice base detected")

    return SliceSet(
        x_left=x_left,
        x_right=x_right,
        b=b,
        x_mid=x_mid,
        y_base=y_base,
        y_ground=y_ground,
        h=h,
        l=l,
        c=c_arr,
        phi=phi_arr,
        tan_phi=tan_phi,
        w=w,
        w_eff=w_eff,
        kh_force=kh_force,
        u=u_press,
        sin_a=sin_a,
        cos_a=cos_a,
        ux=ux,
        uy=uy,
        nx=nx,
        ny=ny,
        y_c=y_c,
        area=area,
        slide_dir=slide_dir,
        x_ref=float(x_ref),
        y_ref=float(y_ref),
        layer_index=layer_idx,
        water_table=None if np.all(np.isnan(water_table)) else water_table,
        notes=notes,
        surface_kind=getattr(surface, "kind", "abstract"),
    )


def _layer_crossings(surface: SlipSurface, soil: SoilModel, x_l: float, x_r: float) -> list[float]:
    """Abscissae where the slip surface crosses an internal layer boundary."""
    from .profile import Slope

    out: list[float] = []
    for i in range(len(soil) - 1):
        layer = soil.layers[i]
        if layer.boundary is None:
            continue
        boundary = Slope(layer.boundary)
        for x in surface.intersections(boundary):
            if x_l < float(x) < x_r:
                out.append(float(x))
    return out


def _unique_sorted(values: Sequence[float], tol: float) -> np.ndarray:
    out: list[float] = []
    for v in sorted(float(x) for x in values):
        if not out or v - out[-1] > tol:
            out.append(v)
    return np.asarray(out, dtype=float)
