"""Geometry package: ground profile, soil model, slip surfaces and slices."""

from __future__ import annotations

from .profile import Slope
from .slices import SliceSet, build_slices
from .soil import SoilLayer, SoilModel
from .surface import PolylineSurface, SlipCircle, SlipSurface

__all__ = [
    "Slope",
    "SoilLayer",
    "SoilModel",
    "SlipCircle",
    "SlipSurface",
    "PolylineSurface",
    "SliceSet",
    "build_slices",
]
