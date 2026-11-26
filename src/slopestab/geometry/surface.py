"""Slip surfaces: circular arcs and user-defined polylines.

A slip surface exposes the few queries the slice builder needs:

* ``intersections(slope)`` - where it crosses the ground surface,
* ``base_y(x)`` - the elevation of the surface at ``x``,
* ``base_slope(x)`` - ``dy/dx`` of the surface at ``x``,
* ``length(x1, x2)`` - arc/chord length between two abscissae,
* ``moment_center()`` - a sensible reference point for moment equilibrium.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

Point = tuple[float, float]


class SlipSurface:
    """Base class for slip surfaces."""

    kind = "abstract"

    def intersections(self, slope) -> np.ndarray:  # pragma: no cover - abstract
        raise NotImplementedError

    def base_y(self, x):  # pragma: no cover - abstract
        raise NotImplementedError

    def base_slope(self, x):  # pragma: no cover - abstract
        raise NotImplementedError

    def length(self, x1: float, x2: float) -> float:  # pragma: no cover - abstract
        raise NotImplementedError

    def moment_center(self) -> Point:  # pragma: no cover - abstract
        raise NotImplementedError

    def sample(
        self, n: int = 200, span: tuple[float, float] | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Points along the surface, for plotting.

        ``span`` restricts the sampling to the part of the surface that lies
        between the ground intersections.
        """
        x0, x1 = span if span is not None else self.x_span()
        xs = np.linspace(x0, x1, n)
        return xs, np.asarray(self.base_y(xs), dtype=float)

    def x_span(self) -> tuple[float, float]:  # pragma: no cover - abstract
        raise NotImplementedError

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"{type(self).__name__}()"


class SlipCircle(SlipSurface):
    """Circular slip surface defined by a centre and a radius."""

    kind = "circle"

    def __init__(self, cx: float, cy: float, r: float) -> None:
        if r <= 0:
            raise ValueError("radius must be positive")
        self.cx = float(cx)
        self.cy = float(cy)
        self.r = float(r)

    # ------------------------------------------------------------------ geometry
    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"SlipCircle(cx={self.cx:.3f}, cy={self.cy:.3f}, r={self.r:.3f})"

    def base_y(self, x):
        x = np.asarray(x, dtype=float)
        u = np.clip(x - self.cx, -self.r, self.r)
        return self.cy - np.sqrt(np.maximum(self.r**2 - u**2, 0.0))

    def base_slope(self, x):
        x = np.asarray(x, dtype=float)
        u = np.clip(x - self.cx, -self.r + 1e-12, self.r - 1e-12)
        v = np.sqrt(np.maximum(self.r**2 - u**2, 1e-12))
        return u / v

    def length(self, x1: float, x2: float) -> float:
        t1 = np.arcsin(np.clip((x1 - self.cx) / self.r, -1.0, 1.0))
        t2 = np.arcsin(np.clip((x2 - self.cx) / self.r, -1.0, 1.0))
        return float(abs(t2 - t1) * self.r)

    def moment_center(self) -> Point:
        return (self.cx, self.cy)

    def x_span(self) -> tuple[float, float]:
        return (self.cx - self.r, self.cx + self.r)

    def lowest_point(self) -> Point:
        return (self.cx, self.cy - self.r)

    def arc_angles(self, x1: float, x2: float) -> tuple[float, float]:
        return (
            float(np.arcsin(np.clip((x1 - self.cx) / self.r, -1.0, 1.0))),
            float(np.arcsin(np.clip((x2 - self.cx) / self.r, -1.0, 1.0))),
        )

    def sample(
        self, n: int = 200, span: tuple[float, float] | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        x0, x1 = span if span is not None else self.x_span()
        t0, t1 = self.arc_angles(x0, x1)
        t = np.linspace(t0, t1, n)
        return self.cx + self.r * np.sin(t), self.cy - self.r * np.cos(t)

    # -------------------------------------------------------------- intersections
    def intersections(self, slope) -> np.ndarray:
        xs: list[float] = []
        for (x1, y1), (x2, y2) in slope.segments():
            dx, dy = x2 - x1, y2 - y1
            a = dx * dx + dy * dy
            b = 2.0 * ((x1 - self.cx) * dx + (y1 - self.cy) * dy)
            c = (x1 - self.cx) ** 2 + (y1 - self.cy) ** 2 - self.r**2
            disc = b * b - 4.0 * a * c
            if disc < 0.0 or a == 0.0:
                continue
            root = np.sqrt(disc)
            for t in ((-b - root) / (2 * a), (-b + root) / (2 * a)):
                if -1e-12 <= t <= 1.0 + 1e-12:
                    xs.append(x1 + t * dx)
        return _unique(sorted(xs), tol=1e-9 * max(1.0, abs(self.cx) + self.r))


class PolylineSurface(SlipSurface):
    """User-defined non-circular slip surface (x-monotone polyline)."""

    kind = "polyline"

    def __init__(self, points: Sequence[Point]) -> None:
        pts = tuple((float(x), float(y)) for x, y in points)
        if len(pts) < 2:
            raise ValueError("a polyline slip surface needs at least 2 points")
        xs = np.array([p[0] for p in pts])
        if np.any(np.diff(xs) <= 0):
            raise ValueError("polyline slip surface x coordinates must be increasing")
        self.points = pts
        self._x = xs
        self._y = np.array([p[1] for p in pts], dtype=float)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"PolylineSurface({len(self.points)} points)"

    def x_span(self) -> tuple[float, float]:
        return (float(self._x[0]), float(self._x[-1]))

    def base_y(self, x):
        return np.interp(x, self._x, self._y)

    def base_slope(self, x):
        x = np.asarray(x, dtype=float)
        idx = np.clip(np.searchsorted(self._x, x, side="right") - 1, 0, len(self._x) - 2)
        dx = self._x[idx + 1] - self._x[idx]
        return (self._y[idx + 1] - self._y[idx]) / dx

    def length(self, x1: float, x2: float) -> float:
        lo, hi = min(x1, x2), max(x1, x2)
        total = 0.0
        for i in range(len(self._x) - 1):
            seg_lo = max(lo, self._x[i])
            seg_hi = min(hi, self._x[i + 1])
            if seg_hi > seg_lo:
                s = (self._y[i + 1] - self._y[i]) / (self._x[i + 1] - self._x[i])
                total += (seg_hi - seg_lo) * np.sqrt(1.0 + s * s)
        return float(total)

    def moment_center(self) -> Point:
        return (float(self._x.mean()), float(self._y.mean()))

    def intersections(self, slope) -> np.ndarray:
        xs: list[float] = []
        for (sx1, sy1), (sx2, sy2) in slope.segments():
            for i in range(len(self._x) - 1):
                px1, py1 = self._x[i], self._y[i]
                px2, py2 = self._x[i + 1], self._y[i + 1]
                lo, hi = max(sx1, px1), min(sx2, px2)
                if hi <= lo:
                    continue
                f1 = _y_on_segment(sx1, sy1, sx2, sy2, lo) - _y_on_segment(px1, py1, px2, py2, lo)
                f2 = _y_on_segment(sx1, sy1, sx2, sy2, hi) - _y_on_segment(px1, py1, px2, py2, hi)
                if f1 == 0.0:
                    xs.append(lo)
                elif f1 * f2 < 0.0:
                    xs.append(lo + (hi - lo) * f1 / (f1 - f2))
                if f2 == 0.0:
                    xs.append(hi)
        return _unique(sorted(xs), tol=1e-9 * max(1.0, abs(self._x[-1]) - abs(self._x[0]) + 1.0))


def _y_on_segment(x1: float, y1: float, x2: float, y2: float, x: float) -> float:
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def _unique(values: Sequence[float], tol: float) -> np.ndarray:
    out: list[float] = []
    for v in values:
        if not out or abs(v - out[-1]) > tol:
            out.append(float(v))
    return np.asarray(out, dtype=float)
