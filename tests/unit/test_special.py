"""Closed-form checks: infinite slope and planar wedge."""

import math

import pytest

import slopestab as sst
from slopestab import Slope, SoilLayer, SoilModel


def test_infinite_slope_matches_hand_calculation():
    beta, z, c, phi, gamma = 30.0, 2.0, 5.0, 25.0, 18.0
    b = math.radians(beta)
    expected = (c + gamma * z * math.cos(b) ** 2 * math.tan(math.radians(phi))) / (
        gamma * z * math.sin(b) * math.cos(b)
    )
    assert sst.infinite_slope_fs(beta, z, c=c, phi=phi, gamma=gamma) == pytest.approx(expected)


def test_infinite_slope_with_seepage_is_lower():
    dry = sst.infinite_slope_fs(30.0, 2.0, c=5.0, phi=25.0, gamma=18.0)
    wet = sst.infinite_slope_fs(
        30.0, 2.0, c=5.0, phi=25.0, gamma=18.0, gamma_sat=20.0, seepage=True
    )
    assert wet < dry


def test_critical_depth_reaches_unit_factor_of_safety():
    beta, c, phi, gamma = 30.0, 5.0, 25.0, 18.0
    depth = sst.critical_depth(beta, c, phi, gamma)
    assert depth is not None
    assert sst.infinite_slope_fs(beta, depth, c=c, phi=phi, gamma=gamma) == pytest.approx(
        1.0, rel=1e-9
    )


def test_critical_depth_is_none_for_strong_soils():
    assert sst.critical_depth(30.0, c=0.0, phi=45.0, gamma=18.0) is None


def test_infinite_slope_rejects_bad_angles():
    with pytest.raises(ValueError):
        sst.infinite_slope_fs(90.0, 2.0)


# ------------------------------------------------------------------------- wedge
SLOPE = Slope([(0, 0), (10, 0), (30, 10), (50, 10)])
SOIL = SoilModel([SoilLayer("s", c=15.0, phi=30.0, gamma=18.0)])


def test_wedge_returns_physical_values():
    res = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0)
    assert res.fs > 1.0
    assert res.weight > 0.0
    assert res.plane_length > 0.0
    assert 0.0 < res.plane_angle_deg < 90.0
    assert res.cohesion_force + res.friction_force == pytest.approx(res.fs * res.driving_force)


def test_wedge_with_water_in_crack_is_less_stable():
    dry = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0, crack_depth=4.0)
    wet = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0, crack_depth=4.0, water_depth=2.0)
    assert wet.fs < dry.fs
    assert wet.water_force_crack == pytest.approx(0.5 * 9.81 * 2.0**2)


def test_deeper_tension_crack_enlarges_the_wedge():
    """A deeper crack lowers the plane, so more soil sits above it."""
    shallow = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0, crack_depth=1.0)
    deep = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0, crack_depth=4.0)
    assert deep.weight > shallow.weight
    assert deep.plane_length > 0.0


def test_wedge_seismic_reduces_stability():
    static = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0)
    seismic = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0, seismic=sst.Seismic(kh=0.2))
    assert seismic.fs < static.fs


def test_wedge_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        sst.planar_wedge(SLOPE, SOIL, x_toe=30.0, x_top=10.0)
    with pytest.raises(ValueError):
        sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0, crack_depth=2.0, water_depth=3.0)


def test_wedge_matches_the_hand_formula():
    """Triangle toe(10,0)-crest(30,10)-crack(40,10): A=50, W=900, L=31.623, psi=18.43."""
    res = sst.planar_wedge(SLOPE, SOIL, x_toe=10.0, x_top=40.0)
    assert res.area == pytest.approx(50.0, rel=1e-9)
    assert res.weight == pytest.approx(900.0, rel=1e-9)
    assert res.plane_length == pytest.approx(31.6228, rel=1e-6)
    assert res.plane_angle_deg == pytest.approx(18.4349, abs=1e-4)
    expected = (
        15.0 * res.plane_length
        + 900.0 * math.cos(math.radians(18.4349)) * math.tan(math.radians(30.0))
    ) / (900.0 * math.sin(math.radians(18.4349)))
    assert res.fs == pytest.approx(expected, rel=1e-3)
