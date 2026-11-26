"""Plane (wedge) sliding with a tension crack - the classic planar failure check.

Used for slopes where a continuous planar discontinuity (bedding plane, fault,
interface) controls stability.  The wedge is bounded by the ground surface, a
vertical tension crack and the sliding plane; water may be present in the
crack and a pseudo-static seismic coefficient may be applied.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from .geometry.profile import Slope
from .geometry.soil import SoilModel
from .loads import Seismic


@dataclass
class WedgeResult:
    """Outcome of a planar wedge analysis."""

    fs: float
    weight: float
    plane_length: float
    plane_angle_deg: float
    water_force_plane: float
    water_force_crack: float
    area: float
    cohesion_force: float
    friction_force: float
    driving_force: float

    def summary(self) -> dict[str, Any]:
        return {
            "method": "planar_wedge",
            "fs": round(self.fs, 6),
            "weight_kN_per_m": round(self.weight, 3),
            "area_m2": round(self.area, 3),
            "plane_length_m": round(self.plane_length, 3),
            "plane_angle_deg": round(self.plane_angle_deg, 3),
            "water_force_plane_kN_per_m": round(self.water_force_plane, 3),
            "water_force_crack_kN_per_m": round(self.water_force_crack, 3),
        }


def planar_wedge(
    slope: Slope,
    soil: SoilModel,
    x_toe: float,
    x_top: float,
    crack_depth: float = 0.0,
    water_depth: float = 0.0,
    seismic: Seismic | None = None,
    gamma_w: float = 9.81,
) -> WedgeResult:
    """Analyse a planar wedge defined by its toe and crest abscissae.

    Parameters
    ----------
    x_toe, x_top:
        Abscissae where the sliding plane meets the ground at the toe and where
        the tension crack is located (``x_top > x_toe``).
    crack_depth:
        Depth of the vertical tension crack measured down from the ground
        surface at ``x_top`` (m).  ``0`` means the plane exits at the surface.
    water_depth:
        Water depth inside the crack (m), measured from the crack bottom;
        must not exceed ``crack_depth``.
    seismic:
        Pseudo-static coefficients (``kh`` acts in the sliding direction).
    """
    if x_top <= x_toe:
        raise ValueError("x_top must be greater than x_toe")
    if crack_depth < 0.0:
        raise ValueError("crack_depth must be non-negative")
    if not 0.0 <= water_depth <= crack_depth + 1e-12:
        raise ValueError("water_depth must be between 0 and crack_depth")
    seismic = seismic or Seismic()

    y_toe = float(slope.y_at(x_toe))
    y_top = float(slope.y_at(x_top))
    p_toe = np.array([float(x_toe), y_toe])
    p_crack = np.array([float(x_top), y_top - float(crack_depth)])

    # ------------------------------------------------------------- wedge geometry
    poly: list[tuple[float, float]] = [(float(x_toe), y_toe)]
    for x in slope.vertices_between(x_toe, x_top):
        poly.append((float(x), float(slope.y_at(float(x)))))
    poly.append((float(x_top), y_top))
    if crack_depth > 0.0:
        poly.append((float(x_top), float(p_crack[1])))
    poly.append((float(x_toe), y_toe))

    area = _polygon_area(poly)
    if area <= 0.0:
        raise ValueError("the wedge has no area; check x_toe / x_top")

    centroid = _polygon_centroid(poly)
    layer = soil.layer_at(centroid[0], centroid[1])
    weight = area * layer.gamma

    plane_vec = p_toe - p_crack
    plane_length = float(np.linalg.norm(plane_vec))
    if plane_length <= 1e-9:
        raise ValueError("the sliding plane is degenerate")
    u_x, u_y = plane_vec / plane_length
    if abs(u_x) < 1e-12:
        raise ValueError("the sliding plane is vertical")
    n_x, n_y = -u_y, u_x
    if n_y < 0.0:
        n_x, n_y = -n_x, -n_y
    slide_dir = float(np.sign(u_x))

    # ------------------------------------------------------------------ forces
    applied = np.array([seismic.kh * weight * slide_dir, -(1.0 + seismic.kv) * weight])
    crack_force = np.array([0.0, 0.0])
    if water_depth > 0.0:
        sin_psi = abs(u_y)
        u_plane = 0.5 * gamma_w * water_depth**2 / max(sin_psi, 1e-9)
        v_crack = 0.5 * gamma_w * water_depth**2
        crack_force = np.array([v_crack * slide_dir, 0.0])
    else:
        u_plane = 0.0
        v_crack = 0.0

    total = applied + crack_force
    n_eff = -(float(total @ np.array([n_x, n_y]))) - u_plane
    driving = float(total @ np.array([u_x, u_y]))
    if driving <= 0.0:
        raise ArithmeticError(
            "the wedge is not kinematically admissible (driving force is not positive); "
            "check that the plane dips towards the free face"
        )

    c, phi = layer.c, layer.phi
    cohesion_force = c * plane_length
    friction_force = max(n_eff, 0.0) * math.tan(math.radians(phi))
    fs = (cohesion_force + friction_force) / driving

    return WedgeResult(
        fs=float(fs),
        weight=float(weight),
        plane_length=plane_length,
        plane_angle_deg=float(math.degrees(math.atan2(-u_y, abs(u_x)))),
        water_force_plane=float(u_plane),
        water_force_crack=float(v_crack),
        area=float(area),
        cohesion_force=float(cohesion_force),
        friction_force=float(friction_force),
        driving_force=float(driving),
    )


def _polygon_area(poly: list[tuple[float, float]]) -> float:
    if len(poly) < 3:
        return 0.0
    x = np.array([p[0] for p in poly], dtype=float)
    y = np.array([p[1] for p in poly], dtype=float)
    return 0.5 * float(np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _polygon_centroid(poly: list[tuple[float, float]]) -> tuple[float, float]:
    x = np.array([p[0] for p in poly], dtype=float)
    y = np.array([p[1] for p in poly], dtype=float)
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    cross = x * y2 - x2 * y
    a = 0.5 * float(np.sum(cross))
    if abs(a) < 1e-12:
        return float(x.mean()), float(y.mean())
    cx = float(np.sum((x + x2) * cross) / (6.0 * a))
    cy = float(np.sum((y + y2) * cross) / (6.0 * a))
    return cx, cy
