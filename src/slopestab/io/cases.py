"""JSON/YAML case files: build a full problem from a plain dictionary."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from ..geometry.profile import Slope
from ..geometry.soil import SoilLayer, SoilModel
from ..geometry.surface import PolylineSurface, SlipCircle
from ..loads import NoWater, PhreaticLine, RuCoefficient, Seismic, SurchargeSet, UniformSurcharge

CASE_TEMPLATE: dict[str, Any] = {
    "method": "bishop",
    "n_slices": 30,
    "slope": {"profile": [[0, 0], [5, 0], [25, 10], [35, 10]]},
    "soil": {
        "layers": [
            {"name": "填土", "c": 10.0, "phi": 25.0, "gamma": 18.5, "gamma_sat": 19.8},
        ]
    },
    "surface": {"type": "circle", "cx": 12.0, "cy": 22.0, "r": 15.0},
    "water": None,
    "seismic": {"kh": 0.0, "kv": 0.0},
    "surcharge": [],
}


def load_case(path: str | Path) -> dict[str, Any]:
    """Read a ``.json`` or ``.yaml`` case file into a dictionary."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError(
                "PyYAML is required for YAML cases; install 'slopestab[yaml]' or use JSON"
            ) from exc
        return dict(yaml.safe_load(text))
    return dict(json.loads(text))


def dump_case(case: dict[str, Any], path: str | Path) -> Path:
    """Write a case dictionary to disk (JSON, or YAML when PyYAML is available)."""
    path = Path(path)
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise ImportError("PyYAML is required to write YAML cases") from exc
        path.write_text(yaml.safe_dump(case, allow_unicode=True, sort_keys=False), encoding="utf-8")
    else:
        path.write_text(json.dumps(case, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


def build_problem(case: dict[str, Any]) -> dict[str, Any]:
    """Turn a case dictionary into ready-to-use slopestab objects.

    Returns a dict with keys ``slope, soil, surface, water, surcharges,
    seismic, method, n_slices``.
    """
    if "slope" not in case or "soil" not in case or "surface" not in case:
        raise ValueError("a case needs 'slope', 'soil' and 'surface' sections")

    profile = [(float(x), float(y)) for x, y in case["slope"]["profile"]]
    slope = Slope(profile)

    layers: list[SoilLayer] = []
    for raw in case["soil"]["layers"]:
        boundary = raw.get("boundary")
        layers.append(
            SoilLayer(
                name=str(raw.get("name", f"layer{len(layers) + 1}")),
                c=float(raw["c"]),
                phi=float(raw["phi"]),
                gamma=float(raw["gamma"]),
                gamma_sat=None if raw.get("gamma_sat") is None else float(raw["gamma_sat"]),
                boundary=None if boundary is None else [(float(x), float(y)) for x, y in boundary],
            )
        )
    soil = SoilModel(layers)

    surf = case["surface"]
    kind = str(surf.get("type", "circle")).lower()
    if kind == "circle":
        surface = SlipCircle(float(surf["cx"]), float(surf["cy"]), float(surf["r"]))
    elif kind in ("polyline", "noncircular", "general"):
        surface = PolylineSurface([(float(x), float(y)) for x, y in surf["points"]])
    else:
        raise ValueError(f"unknown surface type {kind!r}")

    water_raw = case.get("water")
    water = NoWater()
    if water_raw:
        wkind = str(water_raw.get("type", "phreatic")).lower()
        if wkind == "phreatic":
            water = PhreaticLine(
                [(float(x), float(y)) for x, y in water_raw["points"]],
                gamma_w=float(water_raw.get("gamma_w", 9.81)),
            )
        elif wkind == "ru":
            water = RuCoefficient(float(water_raw["ru"]))
        elif wkind in ("none", "dry"):
            water = NoWater()
        else:
            raise ValueError(f"unknown water type {wkind!r}")

    seismic_raw = case.get("seismic") or {}
    seismic = Seismic(kh=float(seismic_raw.get("kh", 0.0)), kv=float(seismic_raw.get("kv", 0.0)))

    sur_raw: Sequence[dict[str, Any]] = case.get("surcharge") or []
    surcharges = SurchargeSet(
        [
            UniformSurcharge(
                q=float(s["q"]),
                x_start=float(s.get("x_start", float("-inf"))),
                x_end=float(s.get("x_end", float("inf"))),
            )
            for s in sur_raw
        ]
    )

    return {
        "slope": slope,
        "soil": soil,
        "surface": surface,
        "water": water,
        "surcharges": surcharges,
        "seismic": seismic,
        "method": str(case.get("method", "bishop")),
        "n_slices": int(case.get("n_slices", 30)),
        "search": case.get("search") or {},
    }
