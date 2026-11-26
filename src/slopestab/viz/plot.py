"""Plotting helpers.

``matplotlib`` is an optional dependency: it is imported inside the functions,
so the numerical core keeps working in environments where only numpy is
available.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np


def plot(
    result: Any,
    path: str | Path | None = None,
    *,
    slope: Any = None,
    soil: Any = None,
    water: Any = None,
    figsize: tuple[float, float] = (9.0, 5.5),
    dpi: int = 150,
    show_slices: bool = True,
    title: str | None = None,
    ax: Any = None,
):
    """Draw the slope section, the slip surface and (optionally) the slices.

    Parameters
    ----------
    result:
        A :class:`~slopestab.methods.base.LEMResult` or
        :class:`~slopestab.search.circular.SearchResult`.
    path:
        If given, the figure is written to this file (PNG/PDF/SVG by suffix).
    slope, soil, water:
        Override the context stored on the result.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Polygon

    lem, search = _unpack(result)
    slope = slope or getattr(search, "slope", None) or getattr(lem, "slope", None)
    soil = soil or getattr(search, "soil", None) or getattr(lem, "soil", None)
    water = water or getattr(search, "water", None) or getattr(lem, "water", None)
    if slope is None:
        raise ValueError("the slope geometry is unknown; pass slope=...")
    surface = getattr(lem, "surface", None)
    slices = lem.slices

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure

    h = max(slope.height, 1e-3)
    px, py = zip(*slope.profile)
    y_floor = min(np.min(py), float(np.min(slices.y_base))) - 0.15 * h

    # soil body
    body = Polygon(
        list(zip(px, py)) + [(px[-1], y_floor), (px[0], y_floor)],
        closed=True,
        facecolor="#d9c9a3",
        edgecolor="none",
        alpha=0.9,
        zorder=0,
    )
    ax.add_patch(body)
    ax.plot(px, py, color="#4a4032", linewidth=2.0, zorder=5)

    # layer boundaries
    if soil is not None:
        for layer in soil.layers:
            boundary = getattr(layer, "boundary", None)
            if boundary is None:
                continue
            bx = np.array([p[0] for p in boundary], dtype=float)
            by = np.array([p[1] for p in boundary], dtype=float)
            inside = (bx >= slope.x_min) & (bx <= slope.x_max)
            if inside.any():
                ax.plot(bx, by, color="#8a7b62", linewidth=1.0, linestyle="--", zorder=3)

    # water table
    wpoints = getattr(water, "points", None)
    if wpoints:
        wx = np.array([p[0] for p in wpoints], dtype=float)
        wy = np.array([p[1] for p in wpoints], dtype=float)
        ax.plot(wx, wy, color="#2a7fb8", linewidth=1.4, zorder=4)
        ax.fill_between(wx, wy, np.interp(wx, px, py), color="#2a7fb8", alpha=0.18, zorder=1)

    # slices
    if show_slices:
        segs = [
            [(float(xl), float(slope.y_at(xl))), (float(xl), float(slices.y_base[i]))]
            for i, xl in enumerate(slices.x_left)
        ]
        segs.append(
            [
                (float(slices.x_right[-1]), float(slope.y_at(slices.x_right[-1]))),
                (float(slices.x_right[-1]), float(slices.y_base[-1])),
            ]
        )
        ax.add_collection(
            LineCollection(segs, colors="#6b8fb5", linewidths=0.8, alpha=0.75, zorder=2)
        )

    # slip surface
    if surface is not None:
        x0, x1 = slices.span
        sx, sy = surface.sample(240, span=(x0, x1))
        ax.plot(sx, sy, color="#c62828", linewidth=2.4, zorder=6)
        ax.plot(
            [x0, x1], [slope.y_at(x0), slope.y_at(x1)], "o", color="#c62828", markersize=5, zorder=7
        )
        gx = np.linspace(x0, x1, 120)
        gy = np.array([slope.y_at(float(x)) for x in gx])
        ax.add_patch(
            Polygon(
                np.concatenate([np.column_stack([sx, sy]), np.column_stack([gx[::-1], gy[::-1]])]),
                closed=True,
                facecolor="#e8743b",
                edgecolor="none",
                alpha=0.25,
                zorder=1,
            )
        )

    label = f"{lem.method}: Fs = {lem.fs:.3f}\n{len(slices)} slices"
    ax.annotate(
        label,
        xy=(0.02, 0.97),
        xycoords="axes fraction",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85, edgecolor="#888"),
    )

    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.03)
    ax.grid(True, alpha=0.25, linestyle=":")
    if title:
        ax.set_title(title)

    if path is not None:
        fig.savefig(Path(path), bbox_inches="tight")
        plt.close(fig)
        return Path(path)
    return fig


def _unpack(result: Any) -> tuple[Any, Any]:
    """Return ``(LEMResult, SearchResult | None)`` for either accepted type."""
    from ..search.circular import SearchResult

    if isinstance(result, SearchResult):
        return result.critical, result
    if hasattr(result, "critical"):
        return result.critical, result
    if hasattr(result, "fs") and hasattr(result, "slices"):
        return result, None
    raise TypeError("plot() expects an LEMResult or a SearchResult")


__all__: Sequence[str] = ("plot",)
