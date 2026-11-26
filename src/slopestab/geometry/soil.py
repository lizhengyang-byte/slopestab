"""Layered soil model: strength and unit weight as functions of position."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

Point = tuple[float, float]


@dataclass
class SoilLayer:
    """One soil layer with Mohr-Coulomb parameters.

    Parameters
    ----------
    name:
        Label used in plots and reports.
    c, phi:
        Effective cohesion (kPa) and friction angle (degrees).
    gamma:
        Bulk unit weight above the water table (kN/m3).
    gamma_sat:
        Saturated unit weight below the water table (kN/m3).  When omitted,
        ``gamma`` is used everywhere (total-stress / dry analysis).
    boundary:
        Optional polyline defining the *bottom* of the layer.  Layers are
        listed from top to bottom, so ``boundary`` must be given for every
        layer except the last one, which extends downwards indefinitely.
    """

    name: str
    c: float
    phi: float
    gamma: float
    gamma_sat: float | None = None
    boundary: Sequence[Point] | None = None

    def __post_init__(self) -> None:
        if self.c < 0:
            raise ValueError(f"layer {self.name!r}: c must be >= 0")
        if not 0.0 <= self.phi < 90.0:
            raise ValueError(f"layer {self.name!r}: phi must be in [0, 90) degrees")
        if self.gamma <= 0:
            raise ValueError(f"layer {self.name!r}: gamma must be > 0")
        if self.gamma_sat is not None and self.gamma_sat <= 0:
            raise ValueError(f"layer {self.name!r}: gamma_sat must be > 0")
        if self.boundary is not None:
            pts = tuple((float(x), float(y)) for x, y in self.boundary)
            if len(pts) < 2:
                raise ValueError(f"layer {self.name!r}: boundary needs >= 2 points")
            xs = np.array([p[0] for p in pts])
            if np.any(np.diff(xs) <= 0.0):
                raise ValueError(f"layer {self.name!r}: boundary x must be increasing")
            self.boundary = pts

    @property
    def tan_phi(self) -> float:
        return float(np.tan(np.radians(self.phi)))

    def unit_weight(self, below_water: bool) -> float:
        """Bulk (or saturated) unit weight."""
        if below_water and self.gamma_sat is not None:
            return float(self.gamma_sat)
        return float(self.gamma)


class SoilModel:
    """Ordered stack of :class:`SoilLayer` objects.

    The first layer is the topmost one.  ``layer_index(x, y)`` answers
    "which layer contains this point", using the layer bottom polylines
    (horizontally extrapolated beyond their own x range).
    """

    def __init__(self, layers: Sequence[SoilLayer]) -> None:
        if not layers:
            raise ValueError("at least one soil layer is required")
        self.layers: tuple[SoilLayer, ...] = tuple(layers)
        for layer in self.layers[:-1]:
            if layer.boundary is None:
                raise ValueError(
                    f"layer {layer.name!r} is not the last layer and needs a 'boundary' polyline"
                )

    # ------------------------------------------------------------------ queries
    def __len__(self) -> int:
        return len(self.layers)

    def bottom_y(self, index: int, x: float) -> float:
        """Elevation of the bottom of layer ``index`` at ``x``."""
        layer = self.layers[index]
        if layer.boundary is None:
            return float("-inf")
        bx = np.array([p[0] for p in layer.boundary], dtype=float)
        by = np.array([p[1] for p in layer.boundary], dtype=float)
        return float(np.interp(x, bx, by))

    def layer_index(self, x: float, y: float) -> int:
        """Index of the layer containing the point ``(x, y)``."""
        for i, layer in enumerate(self.layers):
            if layer.boundary is None:
                return i
            if y >= self.bottom_y(i, x):
                return i
        return len(self.layers) - 1

    def layer_at(self, x: float, y: float) -> SoilLayer:
        return self.layers[self.layer_index(x, y)]

    def params_at(self, x: float, y: float) -> tuple[float, float]:
        """Effective ``(c, phi)`` in kPa / degrees at ``(x, y)``."""
        layer = self.layer_at(x, y)
        return float(layer.c), float(layer.phi)

    def boundaries_between(self, x: float, y_lo: float, y_hi: float) -> list[float]:
        """Elevations of layer boundaries strictly inside ``(y_lo, y_hi)`` at ``x``."""
        cuts = []
        for i in range(len(self.layers) - 1):
            b = self.bottom_y(i, x)
            if y_lo < b < y_hi:
                cuts.append(b)
        return cuts

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        names = ", ".join(layer.name for layer in self.layers)
        return f"SoilModel({len(self.layers)} layers: {names})"
