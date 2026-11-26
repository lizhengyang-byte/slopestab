"""Limit-equilibrium method tests: benchmarks, identities, cross-checks."""

import numpy as np
import pytest

import slopestab as sst
from slopestab import SlipCircle, Slope, SoilLayer, SoilModel, build_slices

# ACADS 1(a): the canonical verification problem (Giam & Donald, 1989).
# Geometry: 10 m high 2H:1V slope; homogeneous soil; published Fs = 1.00 and
# Bishop = 0.987 / Spencer = 0.986 (Slide2 verification manual).
ACADS_SLOPE = Slope([(0, 0), (20, 10), (30, 10)])
ACADS_SOIL = SoilModel([SoilLayer("homogeneous", c=3.0, phi=19.6, gamma=20.0)])
ACADS_CIRCLE = SlipCircle(0.03, 27.70, 27.69)

# A generic homogeneous slope used for cross-method comparisons.
PLAIN_SLOPE = Slope.simple(height=10, ratio=2.0, crest_length=20, toe_length=10)
PLAIN_SOIL = SoilModel([SoilLayer("s", c=10.0, phi=20.0, gamma=18.0, gamma_sat=20.0)])
PLAIN_CIRCLE = SlipCircle(10.03, 27.70, 27.69)


# --------------------------------------------------------------------- benchmarks
def test_acads_1a_bishop_matches_published_value():
    res = sst.analyze("bishop", ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=30)
    assert res.fs == pytest.approx(0.987, abs=0.02)
    assert res.converged


def test_acads_1a_spencer_matches_published_value():
    res = sst.analyze("spencer", ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=30)
    assert res.fs == pytest.approx(0.986, abs=0.02)


def test_acads_1a_morgenstern_price_matches_published_value():
    res = sst.analyze("morgenstern_price", ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=30)
    assert res.fs == pytest.approx(0.986, abs=0.02)


def test_acads_1a_janbu_correction_reaches_the_expected_range():
    raw = sst.analyze("janbu", ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=30)
    corrected = sst.analyze("janbu", ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=30, f0=1.05)
    assert raw.fs < corrected.fs
    assert corrected.fs == pytest.approx(0.990, abs=0.03)


def test_acads_1a_all_methods_stay_within_a_sane_band():
    results = sst.compare(
        ["fellenius", "bishop", "janbu", "spencer", "morgenstern_price"],
        ACADS_SLOPE,
        ACADS_SOIL,
        ACADS_CIRCLE,
        n_slices=30,
    )
    for name, res in results.items():
        assert 0.90 < res.fs < 1.06, f"{name} drifted to {res.fs:.3f}"


def test_ordinary_method_is_the_lowest_bound():
    results = sst.compare(["fellenius", "bishop", "spencer"], PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    assert results["fellenius"].fs < results["bishop"].fs
    assert results["bishop"].fs <= results["spencer"].fs * 1.05


# --------------------------------------------------------------------- identities
def test_undrained_slices_make_bishop_equal_to_fellenius():
    soil = SoilModel([SoilLayer("undrained", c=40.0, phi=0.0, gamma=18.0)])
    bishop = sst.analyze("bishop", ACADS_SLOPE, soil, ACADS_CIRCLE)
    fellenius = sst.analyze("fellenius", ACADS_SLOPE, soil, ACADS_CIRCLE)
    assert bishop.fs == pytest.approx(fellenius.fs, rel=1e-3)


def test_undrained_factor_of_safety_equals_the_moment_identity():
    """phi = 0: Fs = c*L*R / sum(W*(x - x0)) is the exact moment solution."""
    soil = SoilModel([SoilLayer("undrained", c=40.0, phi=0.0, gamma=18.0)])
    slices = build_slices(ACADS_SLOPE, soil, ACADS_CIRCLE, n_slices=30)
    length = ACADS_CIRCLE.length(*slices.span)
    expected = (
        40.0 * length * ACADS_CIRCLE.r / float(np.sum(slices.w * (slices.x_mid - ACADS_CIRCLE.cx)))
    )
    res = sst.analyze("fellenius", ACADS_SLOPE, soil, ACADS_CIRCLE, n_slices=30)
    assert res.fs == pytest.approx(expected, rel=2e-3)


def test_bishop_and_spencer_agree_on_homogeneous_slopes():
    bishop = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    spencer = sst.analyze("spencer", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    mp = sst.analyze("morgenstern_price", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    assert abs(bishop.fs - spencer.fs) / spencer.fs < 0.02
    assert abs(mp.fs - spencer.fs) / spencer.fs < 0.02


def test_mp_with_constant_function_reproduces_spencer():
    spencer = sst.analyze("spencer", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    mp = sst.analyze("morgenstern_price", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE, f="constant")
    assert mp.fs == pytest.approx(spencer.fs, rel=1e-3)
    assert mp.info["lambda"] == pytest.approx(spencer.info["lambda"], abs=1e-3)


def test_rigorous_methods_satisfy_equilibrium():
    for method in ("spencer", "morgenstern_price"):
        res = sst.analyze(method, PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
        assert abs(res.residuals["Fx"]) < 0.05
        assert abs(res.residuals["M"]) < 0.05


def test_slice_count_does_not_change_the_answer_materially():
    coarse = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE, n_slices=12)
    fine = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE, n_slices=60)
    assert abs(coarse.fs - fine.fs) / fine.fs < 0.02


# ------------------------------------------------------------------ loads effects
def test_water_reduces_stability():
    dry = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    wet = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE, water=sst.RuCoefficient(0.3))
    assert wet.fs < dry.fs


def test_horizontal_seismic_reduces_stability():
    static = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    seismic = sst.analyze(
        "bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE, seismic=sst.Seismic(kh=0.15)
    )
    assert seismic.fs < static.fs


def test_downward_vertical_seismic_reduces_stability():
    static = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    seismic = sst.analyze(
        "bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE, seismic=sst.Seismic(kh=0.0, kv=0.2)
    )
    assert seismic.fs < static.fs


def test_surcharge_reduces_stability():
    static = sst.analyze("bishop", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)
    loaded = sst.analyze(
        "bishop",
        PLAIN_SLOPE,
        PLAIN_SOIL,
        PLAIN_CIRCLE,
        surcharge=sst.UniformSurcharge(q=30.0, x_start=25.0, x_end=60.0),
    )
    assert loaded.fs < static.fs


def test_layered_profile_still_converges():
    soil = SoilModel(
        [
            SoilLayer("fill", c=5, phi=25, gamma=18, boundary=[(0, 3), (60, 3)]),
            SoilLayer("clay", c=15, phi=18, gamma=19),
        ]
    )
    res = sst.analyze("morgenstern_price", PLAIN_SLOPE, soil, PLAIN_CIRCLE)
    assert 0.5 < res.fs < 10.0


# ---------------------------------------------------------------------- failures
def test_unknown_method_is_rejected():
    with pytest.raises(ValueError, match="unknown method"):
        sst.analyze("taylor", PLAIN_SLOPE, PLAIN_SOIL, PLAIN_CIRCLE)


def test_non_circular_surface_is_rejected_by_rigorous_methods():
    surface = sst.PolylineSurface([(10, 0), (18, 2), (26, 5), (34, 10)])
    with pytest.raises(sst.GeometryError):
        sst.analyze("spencer", PLAIN_SLOPE, PLAIN_SOIL, surface)


def test_compare_returns_every_method():
    results = sst.compare(
        ["fellenius", "bishop", "janbu", "spencer", "morgenstern_price"],
        PLAIN_SLOPE,
        PLAIN_SOIL,
        PLAIN_CIRCLE,
    )
    assert set(results) == {"fellenius", "bishop", "janbu", "spencer", "morgenstern_price"}
