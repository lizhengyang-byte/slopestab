"""External loads and hydraulic conditions applied to a slope.

Three pore-pressure models are provided:

* :class:`NoWater` - dry / total-stress analysis, ``u = 0``.
* :class:`RuCoefficient` - pore-pressure ratio ``r_u``, ``u = r_u * gamma * h``.
* :class:`PhreaticLine` - hydrostatic pressure below a water table polyline,
  ``u = gamma_w * (y_water - y_base)``.

All pressures are in kPa (kN/m2) and all forces per unit width of slope in kN/m.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

GAMMA_W = 9.81  # kN/m3

Point = tuple[float, float]


class WaterCondition:
    """Base class for pore-pressure models."""

    kind = "none"

    def pressure(self, x: float, y_base: float, weight: float, width: float) -> float:
        """Pore pressure (kPa) acting on the slice base."""
        return 0.0

    def water_y(self, x) -> float | None:
        """Water-table elevation at ``x`` (``None`` when there is no table)."""
        return None

    def crossing_x(self, surface) -> list[float]:
        """Abscissae where the water table crosses the slip surface."""
        return []


class NoWater(WaterCondition):
    """Dry slope (or a total-stress analysis with no pore pressures)."""

    kind = "none"

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return "NoWater()"


@dataclass
class RuCoefficient(WaterCondition):
    """Pore-pressure ratio model.

    ``u = r_u * sigma_v`` with ``sigma_v`` taken as the slice weight divided by
    its base width, i.e. ``u * b = r_u * W``.  This is the standard ``r_u``
    shortcut used for rapid drawdown and undrained loading screening.
    """

    ru: float

    kind = "ru"

    def __post_init__(self) -> None:
        if not 0.0 <= self.ru < 1.0:
            raise ValueError("r_u must be in [0, 1)")

    def pressure(self, x: float, y_base: float, weight: float, width: float) -> float:
        return float(self.ru * weight / width)


class PhreaticLine(WaterCondition):
    """Hydrostatic pore pressure below a phreatic line.

    Parameters
    ----------
    points:
        Water-table polyline ``[(x, y), ...]`` with increasing ``x``.  It is
        linearly extrapolated beyond its own range.
    gamma_w:
        Unit weight of water (kN/m3), 9.81 by default.
    """

    kind = "phreatic"

    def __init__(self, points: Sequence[Point], gamma_w: float = GAMMA_W) -> None:
        pts = tuple((float(x), float(y)) for x, y in points)
        if len(pts) < 2:
            raise ValueError("phreatic line needs at least 2 points")
        xs = np.array([p[0] for p in pts])
        if np.any(np.diff(xs) <= 0):
            raise ValueError("phreatic line x coordinates must be increasing")
        self.points = pts
        self.gamma_w = float(gamma_w)
        self._x = xs
        self._y = np.array([p[1] for p in pts], dtype=float)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"PhreaticLine({len(self.points)} points)"

    def water_y(self, x):
        return np.interp(x, self._x, self._y)

    def pressure(self, x: float, y_base: float, weight: float, width: float) -> float:
        depth = float(np.interp(x, self._x, self._y)) - y_base
        return float(self.gamma_w * depth) if depth > 0.0 else 0.0

    def crossing_x(self, surface) -> list[float]:
        # Re-use the surface/ground intersection machinery by treating the
        # water table as a one-segment "ground profile".
        from .geometry.profile import Slope

        return [float(v) for v in surface.intersections(Slope(self.points))]


@dataclass
class Seismic:
    """Pseudo-static seismic coefficients.

    ``kh`` acts horizontally in the direction of sliding; ``kv`` is positive
    downwards (i.e. a destabilising vertical acceleration, as used by USACE).
    """

    kh: float = 0.0
    kv: float = 0.0

    def __post_init__(self) -> None:
        if abs(self.kh) > 1.0 or abs(self.kv) > 1.0:
            raise ValueError("seismic coefficients are expected to be within [-1, 1]")

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"Seismic(kh={self.kh:g}, kv={self.kv:g})"


@dataclass
class UniformSurcharge:
    """Uniform surcharge ``q`` (kPa) applied between ``x_start`` and ``x_end``."""

    q: float
    x_start: float = float("-inf")
    x_end: float = float("inf")

    def __post_init__(self) -> None:
        if self.q < 0:
            raise ValueError("surcharge intensity must be non-negative")
        if self.x_end <= self.x_start:
            raise ValueError("x_end must be greater than x_start")

    def applies_at(self, x: float) -> bool:
        return self.x_start <= x <= self.x_end


class SurchargeSet:
    """Collection of surcharges, evaluated at a given ``x``."""

    def __init__(self, surcharges: Sequence[UniformSurcharge] | None = None) -> None:
        self.items: tuple[UniformSurcharge, ...] = tuple(surcharges or ())

    def q_at(self, x: float) -> float:
        return float(sum(s.q for s in self.items if s.applies_at(x)))

    def edges(self) -> list[float]:
        out: list[float] = []
        for s in self.items:
            if np.isfinite(s.x_start):
                out.append(float(s.x_start))
            if np.isfinite(s.x_end):
                out.append(float(s.x_end))
        return out

    def __len__(self) -> int:
        return len(self.items)
