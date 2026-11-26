"""slopestab - 2D slope stability analysis by limit equilibrium methods.

Quick start
-----------
>>> import slopestab as sst
>>> slope = sst.Slope.simple(height=10, ratio=2.0)
>>> soil = sst.SoilModel([sst.SoilLayer("clay", c=3, phi=19.6, gamma=20)])
>>> res = sst.analyze("bishop", slope, soil, sst.SlipCircle(12, 22, 15))
>>> round(res.fs, 3)
1.0...

Seven analysis tools are exposed: the five limit-equilibrium methods
(``fellenius``, ``bishop``, ``janbu``, ``spencer``, ``morgenstern_price``),
the infinite-slope closed form and the planar-wedge check, plus a circular
slip-surface search and reporting helpers.
"""

from __future__ import annotations

__version__ = "0.1.0"

from .analysis import METHOD_NAMES, analyze, compare, run_on_slices
from .exceptions import ConvergenceError, GeometryError, SlopestabError
from .geometry.profile import Slope
from .geometry.slices import SliceSet, build_slices
from .geometry.soil import SoilLayer, SoilModel
from .geometry.surface import PolylineSurface, SlipCircle, SlipSurface
from .infinite import critical_depth, infinite_slope_fs
from .io.cases import build_problem, dump_case, load_case
from .loads import (
    NoWater,
    PhreaticLine,
    RuCoefficient,
    Seismic,
    SurchargeSet,
    UniformSurcharge,
    WaterCondition,
)
from .methods.base import LEMResult
from .methods.bishop import bishop
from .methods.fellenius import fellenius
from .methods.janbu import janbu
from .methods.morgenstern_price import morgenstern_price
from .methods.spencer import spencer
from .search.circular import SearchGrid, SearchResult, search_circular
from .viz.plot import plot
from .viz.report import report_html
from .wedge import WedgeResult, planar_wedge

__all__ = [
    "__version__",
    "METHOD_NAMES",
    # geometry
    "Slope",
    "SoilLayer",
    "SoilModel",
    "SlipCircle",
    "SlipSurface",
    "PolylineSurface",
    "SliceSet",
    "build_slices",
    # loads
    "NoWater",
    "RuCoefficient",
    "PhreaticLine",
    "Seismic",
    "UniformSurcharge",
    "SurchargeSet",
    "WaterCondition",
    # analysis
    "analyze",
    "compare",
    "run_on_slices",
    "LEMResult",
    "fellenius",
    "bishop",
    "janbu",
    "spencer",
    "morgenstern_price",
    # search
    "search_circular",
    "SearchGrid",
    "SearchResult",
    # special cases
    "infinite_slope_fs",
    "critical_depth",
    "planar_wedge",
    "WedgeResult",
    # reporting
    "plot",
    "report_html",
    # io
    "load_case",
    "dump_case",
    "build_problem",
    # exceptions
    "SlopestabError",
    "GeometryError",
    "ConvergenceError",
]
