"""Public entry points: ``analyze`` and ``compare``."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from .geometry.profile import Slope
from .geometry.slices import SliceSet, build_slices
from .geometry.soil import SoilModel
from .geometry.surface import SlipSurface
from .loads import Seismic, SurchargeSet, UniformSurcharge, WaterCondition
from .methods.base import LEMResult
from .methods.bishop import bishop
from .methods.fellenius import fellenius
from .methods.janbu import janbu
from .methods.morgenstern_price import morgenstern_price
from .methods.spencer import spencer

_METHODS = {
    "fellenius": fellenius,
    "ordinary": fellenius,
    "swedish": fellenius,
    "bishop": bishop,
    "bishop_simplified": bishop,
    "janbu": janbu,
    "janbu_simplified": janbu,
    "spencer": spencer,
    "morgenstern_price": morgenstern_price,
    "morgenstern-price": morgenstern_price,
    "mp": morgenstern_price,
    "gle": morgenstern_price,
}

METHOD_NAMES = ("fellenius", "bishop", "janbu", "spencer", "morgenstern_price")


def _as_surcharges(
    surcharge: UniformSurcharge | Iterable[UniformSurcharge] | None,
    surcharges: SurchargeSet | None,
) -> SurchargeSet | None:
    if surcharges is not None and surcharge is not None:
        raise ValueError("pass either 'surcharge' or 'surcharges', not both")
    if surcharges is not None:
        return surcharges
    if surcharge is None:
        return None
    if isinstance(surcharge, UniformSurcharge):
        return SurchargeSet([surcharge])
    return SurchargeSet(list(surcharge))


def analyze(
    method: str = "bishop",
    slope: Slope | None = None,
    soil: SoilModel | None = None,
    surface: SlipSurface | None = None,
    water: WaterCondition | None = None,
    surcharge: UniformSurcharge | Iterable[UniformSurcharge] | None = None,
    surcharges: SurchargeSet | None = None,
    seismic: Seismic | None = None,
    n_slices: int = 30,
    **kwargs: Any,
) -> LEMResult:
    """Analyse one slip surface with a chosen limit-equilibrium method.

    Parameters
    ----------
    method:
        ``"fellenius"`` (``"ordinary"``, ``"swedish"``), ``"bishop"``,
        ``"janbu"``, ``"spencer"`` or ``"morgenstern_price"`` (``"mp"``,
        ``"gle"``).
    slope, soil, surface:
        Geometry, materials and the slip surface.  Required.
    water:
        Pore-pressure model; dry when omitted.
    surcharge, surcharges:
        One or more uniform surcharges.
    seismic:
        Pseudo-static coefficients.
    n_slices:
        Number of slices before extra cuts are inserted.
    **kwargs:
        Forwarded to the method (e.g. ``f0`` for Janbu, ``f`` for
        Morgenstern-Price).

    Examples
    --------
    >>> res = analyze("bishop", slope, soil, SlipCircle(12, 22, 15))
    >>> res.fs
    1.42...
    """
    if slope is None or soil is None or surface is None:
        raise ValueError("slope, soil and surface are required")
    key = method.lower().replace(" ", "_")
    if key not in _METHODS:
        raise ValueError(f"unknown method {method!r}; choose from {sorted(_METHODS)}")

    slices = build_slices(
        slope,
        soil,
        surface,
        water=water,
        surcharges=_as_surcharges(surcharge, surcharges),
        seismic=seismic,
        n_slices=n_slices,
    )
    return run_on_slices(key, slices, **kwargs)


def run_on_slices(method: str, slices: SliceSet, **kwargs: Any) -> LEMResult:
    """Run a method on an already built slice set (avoids re-discretising)."""
    key = method.lower().replace(" ", "_")
    if key not in _METHODS:
        raise ValueError(f"unknown method {method!r}; choose from {sorted(_METHODS)}")
    return _METHODS[key](slices, **kwargs)


def compare(
    methods: Iterable[str] = METHOD_NAMES,
    slope: Slope | None = None,
    soil: SoilModel | None = None,
    surface: SlipSurface | None = None,
    **kwargs: Any,
) -> dict[str, LEMResult]:
    """Run several methods on the *same* slice set and return their results.

    Discretising once and sharing the slices is what makes a five-method
    cross-check meaningful: any difference between the factors of safety comes
    from the equilibrium assumptions, not from a different slice layout.
    """
    if slope is None or soil is None or surface is None:
        raise ValueError("slope, soil and surface are required")
    n_slices = int(kwargs.pop("n_slices", 30))
    water = kwargs.pop("water", None)
    seismic = kwargs.pop("seismic", None)
    surcharge = kwargs.pop("surcharge", None)
    surcharges = kwargs.pop("surcharges", None)
    slices = build_slices(
        slope,
        soil,
        surface,
        water=water,
        surcharges=_as_surcharges(surcharge, surcharges),
        seismic=seismic,
        n_slices=n_slices,
    )
    out: dict[str, LEMResult] = {}
    for name in methods:
        key = name.lower().replace(" ", "_")
        out[key] = _METHODS[key](slices, **kwargs)
    return out
