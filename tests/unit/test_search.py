"""Slip-surface search tests."""

import numpy as np
import pytest

import slopestab as sst
from slopestab import SearchGrid, SlipCircle, Slope, SoilLayer, SoilModel

SLOPE = Slope.simple(height=10, ratio=2.0, crest_length=20, toe_length=10)
SOIL = SoilModel([SoilLayer("s", c=10.0, phi=20.0, gamma=18.0)])
SMALL_GRID = SearchGrid(
    cx_range=(0.0, 25.0),
    cy_range=(11.0, 26.0),
    r_range=(10.0, 26.0),
    n_cx=8,
    n_cy=6,
    n_r=6,
)


def test_search_finds_a_surface_and_ranks_it_first():
    best = sst.search_circular(SLOPE, SOIL, method="bishop", grid=SMALL_GRID, n_slices=20)
    assert best.fs > 0.0
    assert best.candidates
    assert best.fs <= best.candidates[0][3] + 1e-6  # refinement may improve it
    assert best.n_evaluated > 0


def test_critical_surface_beats_a_hand_picked_circle():
    picked = sst.analyze("bishop", SLOPE, SOIL, SlipCircle(10.03, 27.70, 27.69)).fs
    best = sst.search_circular(SLOPE, SOIL, method="bishop", grid=SMALL_GRID, n_slices=20)
    assert best.fs <= picked + 1e-6


def test_refinement_never_makes_the_answer_worse():
    coarse = sst.search_circular(
        SLOPE, SOIL, method="bishop", grid=SMALL_GRID, n_slices=20, refine=False
    )
    refined = sst.search_circular(
        SLOPE, SOIL, method="bishop", grid=SMALL_GRID, n_slices=20, refine=True
    )
    assert refined.fs <= coarse.fs + 1e-6


def test_search_result_exports_candidates():
    best = sst.search_circular(SLOPE, SOIL, method="bishop", grid=SMALL_GRID, n_slices=20)
    df = best.to_dataframe()
    assert {"cx", "cy", "r", "fs"} <= set(df.columns)
    assert len(df) == len(best.candidates)


def test_unreachable_grid_raises_a_helpful_error():
    grid = SearchGrid(
        cx_range=(200.0, 220.0),
        cy_range=(200.0, 220.0),
        r_range=(10.0, 12.0),
        n_cx=4,
        n_cy=4,
        n_r=4,
    )
    with pytest.raises(ValueError, match="no admissible"):
        sst.search_circular(SLOPE, SOIL, grid=grid)


def test_search_is_deterministic():
    a = sst.search_circular(SLOPE, SOIL, method="bishop", grid=SMALL_GRID, n_slices=15)
    b = sst.search_circular(SLOPE, SOIL, method="bishop", grid=SMALL_GRID, n_slices=15)
    assert a.fs == pytest.approx(b.fs, rel=1e-12)
    assert np.allclose(a.critical.slices.x_mid, b.critical.slices.x_mid)
