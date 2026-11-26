"""Geometry tests: profile, surfaces, intersections."""

import numpy as np
import pytest

from slopestab import PolylineSurface, SlipCircle, Slope


def test_slope_simple_geometry():
    slope = Slope.simple(height=10.0, ratio=2.0, crest_length=10.0, toe_length=5.0)
    assert slope.profile == ((0.0, 0.0), (5.0, 0.0), (25.0, 10.0), (35.0, 10.0))
    assert slope.height == pytest.approx(10.0)
    assert slope.toe == (0.0, 0.0)
    assert slope.face_angle_deg == pytest.approx(26.565, abs=1e-3)


def test_slope_rejects_non_monotone_profile():
    with pytest.raises(ValueError):
        Slope([(0, 0), (10, 5), (5, 8)])


def test_y_at_is_clamped_outside_the_profile():
    slope = Slope([(0, 0), (10, 0), (20, 10)])
    assert slope.y_at(-5) == pytest.approx(0.0)
    assert slope.y_at(25) == pytest.approx(10.0)
    assert slope.y_at(15) == pytest.approx(5.0)


def test_circle_intersects_ground_twice():
    slope = Slope([(0, 0), (20, 10), (30, 10)])
    circle = SlipCircle(0.03, 27.70, 27.69)
    ints = circle.intersections(slope)
    assert ints.size == 2
    assert ints[0] == pytest.approx(0.0, abs=0.05)
    assert ints[1] == pytest.approx(21.32, abs=0.05)


def test_circle_base_geometry_consistency():
    circle = SlipCircle(12.0, 22.0, 15.0)
    xs = np.linspace(1.0, 20.0, 7)
    ys = circle.base_y(xs)
    # every sampled point must sit on the circle
    assert np.allclose(np.hypot(xs - circle.cx, ys - circle.cy), circle.r, atol=1e-9)
    # the analytic slope must match a finite difference of base_y
    h = 1e-6
    fd = (circle.base_y(xs + h) - circle.base_y(xs - h)) / (2 * h)
    assert np.allclose(circle.base_slope(xs), fd, atol=1e-4)


def test_arc_length_matches_analytic_value():
    circle = SlipCircle(0.0, 0.0, 10.0)
    x1, x2 = -6.0, 6.0
    expected = 10.0 * (np.arcsin(0.6) - np.arcsin(-0.6))
    assert circle.length(x1, x2) == pytest.approx(expected, rel=1e-9)


def test_polyline_surface_intersections():
    slope = Slope([(0, 0), (10, 0), (30, 10), (50, 10)])
    surface = PolylineSurface([(10, 0), (18, 2.0), (26, 5.0), (34, 10)])
    ints = surface.intersections(slope)
    assert ints.size == 2
    assert ints[0] == pytest.approx(10.0)
    assert ints[-1] == pytest.approx(34.0)
    # both ends lie on the ground surface
    for x in ints:
        assert surface.base_y(x) == pytest.approx(slope.y_at(x), abs=1e-9)


def test_polyline_rejects_non_monotone_points():
    with pytest.raises(ValueError):
        PolylineSurface([(0, 0), (5, 2), (4, 3)])
