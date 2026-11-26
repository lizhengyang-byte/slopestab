"""Ground surface geometry: a piecewise-linear slope profile.

All lengths are in metres and all coordinates follow the usual geotechnical
convention: ``x`` points to the right, ``y`` points upwards.

A :class:`Slope` is nothing but an x-monotone polyline.  Anything more
complicated (benches, berms, layered ground) is expressed by adding vertices,
which keeps the geometry exact and the intersection code simple.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np

Point = tuple[float, float]


@dataclass
class Slope:
    """Piecewise-linear ground surface.

    Parameters
    ----------
    profile:
        Vertices ``[(x0, y0), (x1, y1), ...]`` with strictly increasing ``x``.
        The profile must cover the whole zone of interest (the search domain
        and every slip surface are clipped to ``[x_min, x_max]``).

    Examples
    --------
    A 10 m high 2H:1V slope with a 10 m wide crest and a 5 m wide toe:

    >>> Slope.simple(height=10.0, ratio=2.0, crest_length=10.0, toe_length=5.0).profile
    ((0.0, 0.0), (5.0, 0.0), (25.0, 10.0), (35.0, 10.0))
    """

    profile: Sequence[Point]

    def __post_init__(self) -> None:
        pts = tuple((float(x), float(y)) for x, y in self.profile)
        if len(pts) < 2:
            raise ValueError("profile needs at least two vertices")
        xs = np.array([p[0] for p in pts], dtype=float)
        ys = np.array([p[1] for p in pts], dtype=float)
        if np.any(np.diff(xs) <= 0.0):
            raise ValueError("profile x coordinates must be strictly increasing")
        if not (np.all(np.isfinite(xs)) and np.all(np.isfinite(ys))):
            raise ValueError("profile coordinates must be finite")
        self.profile = pts
        self._x = xs
        self._y = ys

    # ------------------------------------------------------------------ geometry
    @classmethod
    def simple(
        cls,
        height: float = 10.0,
        ratio: float = 2.0,
        crest_length: float = 10.0,
        toe_length: float = 5.0,
    ) -> Slope:
        """Build the classic textbook slope: flat toe, straight face, flat crest.

        Parameters
        ----------
        height:
            Slope height ``H`` (m).
        ratio:
            Horizontal : vertical ratio of the face, e.g. ``2.0`` for 2H:1V.
        crest_length, toe_length:
            Horizontal extent of the flat crest and toe (m).
        """
        if height <= 0 or ratio <= 0:
            raise ValueError("height and ratio must be positive")
        if crest_length < 0 or toe_length < 0:
            raise ValueError("crest_length and toe_length must be non-negative")
        return cls(
            (
                (0.0, 0.0),
                (toe_length, 0.0),
                (toe_length + ratio * height, height),
                (toe_length + ratio * height + crest_length, height),
            )
        )

    @property
    def x_min(self) -> float:
        return float(self._x[0])

    @property
    def x_max(self) -> float:
        return float(self._x[-1])

    @property
    def y_min(self) -> float:
        return float(self._y.min())

    @property
    def y_max(self) -> float:
        return float(self._y.max())

    @property
    def height(self) -> float:
        """Total relief of the profile (crest elevation minus lowest point)."""
        return self.y_max - self.y_min

    @property
    def toe(self) -> Point:
        """Lowest vertex of the profile (the toe of the slope)."""
        i = int(np.argmin(self._y))
        return (float(self._x[i]), float(self._y[i]))

    @property
    def crest(self) -> Point:
        """Highest vertex of the profile."""
        i = int(np.argmax(self._y))
        return (float(self._x[i]), float(self._y[i]))

    @property
    def face_angle_deg(self) -> float:
        """Steepest inclination of the profile, in degrees (0 = horizontal)."""
        dx = np.diff(self._x)
        dy = np.diff(self._y)
        return float(np.degrees(np.arctan(np.max(np.abs(dy / dx))))) if len(dx) else 0.0

    def y_at(self, x: float | np.ndarray) -> float | np.ndarray:
        """Ground elevation at ``x`` (linear interpolation, clamped at the ends)."""
        return np.interp(x, self._x, self._y)

    def vertices_between(self, x0: float, x1: float) -> np.ndarray:
        """Profile vertices whose x lies strictly inside ``(x0, x1)``."""
        lo, hi = min(x0, x1), max(x0, x1)
        return self._x[(self._x > lo) & (self._x < hi)]

    def segments(self) -> Iterable[tuple[Point, Point]]:
        """Iterate over ``((x1, y1), (x2, y2))`` segments of the profile."""
        for i in range(len(self.profile) - 1):
            yield self.profile[i], self.profile[i + 1]

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"Slope(x=[{self.x_min:.2f}, {self.x_max:.2f}], height={self.height:.2f}, "
            f"face={self.face_angle_deg:.1f} deg)"
        )
