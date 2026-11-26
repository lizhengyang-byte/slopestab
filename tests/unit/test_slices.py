"""Slice builder tests: masses, layer cuts, water and seismic loads."""

import numpy as np
import pytest

from slopestab import (
    PhreaticLine,
    RuCoefficient,
    Seismic,
    SlipCircle,
    Slope,
    SoilLayer,
    SoilModel,
    SurchargeSet,
    UniformSurcharge,
    build_slices,
)

SLOPE = Slope([(0, 0), (20, 10), (30, 10)])
HOMOGENEOUS = SoilModel([SoilLayer("hom", c=3.0, phi=19.6, gamma=20.0)])
ACADS_CIRCLE = SlipCircle(0.03, 27.70, 27.69)


def test_slice_mass_matches_published_area():
    """ACADS 1(a): the published sliding-mass area is 47.17 m2 (reconstruction 48.18)."""
    slices = build_slices(SLOPE, HOMOGENEOUS, ACADS_CIRCLE, n_slices=40)
    assert slices.area.sum() == pytest.approx(48.18, rel=0.02)
    assert slices.total_weight == pytest.approx(20.0 * slices.area.sum(), rel=1e-9)


def test_slices_are_contiguous_and_ordered():
    slices = build_slices(SLOPE, HOMOGENEOUS, ACADS_CIRCLE, n_slices=20)
    assert np.allclose(slices.x_right[:-1], slices.x_left[1:])
    assert np.allclose(slices.b, slices.x_right - slices.x_left)
    assert (slices.b > 0).all()


def test_base_direction_cosines_are_normalised():
    slices = build_slices(SLOPE, HOMOGENEOUS, ACADS_CIRCLE)
    assert np.allclose(slices.cos_a**2 + slices.sin_a**2, 1.0)
    assert (slices.cos_a > 0).all()
    assert np.allclose(slices.ux**2 + slices.uy**2, 1.0)
    assert np.allclose(slices.nx**2 + slices.ny**2, 1.0)
    assert (slices.ny > 0).all()
    # tangent and normal are orthogonal
    assert np.allclose(slices.ux * slices.nx + slices.uy * slices.ny, 0.0)


def test_layer_boundaries_are_honoured_by_cuts():
    soil = SoilModel(
        [
            SoilLayer("fill", c=5, phi=25, gamma=18, boundary=[(0, 3), (40, 3)]),
            SoilLayer("clay", c=15, phi=18, gamma=19),
        ]
    )
    circle = SlipCircle(0.03, 27.70, 27.69)
    slices = build_slices(SLOPE, soil, circle, n_slices=20)
    # every base must sit in the layer recorded for it
    for i in range(len(slices)):
        layer = soil.layer_at(slices.x_mid[i], slices.y_base[i])
        assert layer.c == pytest.approx(slices.c[i])
        assert layer.phi == pytest.approx(slices.phi[i])


def test_phreatic_line_lowers_factor_of_safety():
    from slopestab import analyze

    soil = SoilModel([SoilLayer("s", c=10, phi=20, gamma=18, gamma_sat=20)])
    slope = Slope.simple(height=10, ratio=2.0, crest_length=20, toe_length=10)
    circle = SlipCircle(10.03, 27.70, 27.69)
    dry = analyze("bishop", slope, soil, circle).fs
    wet = analyze("bishop", slope, soil, circle, water=PhreaticLine([(10, 0), (30, 9), (50, 9)])).fs
    assert wet < dry


def test_ru_model_reduces_normal_force():
    soil = SoilModel([SoilLayer("s", c=10, phi=20, gamma=18)])
    slope = Slope.simple(height=10, ratio=2.0)
    circle = SlipCircle(10.03, 27.70, 27.69)
    dry = build_slices(slope, soil, circle)
    wet = build_slices(slope, soil, circle, water=RuCoefficient(0.4))
    assert np.allclose(wet.u * wet.b, 0.4 * wet.w)
    assert np.all(dry.u == 0.0)


def test_seismic_arrays_are_signed_with_slide_direction():
    soil = SoilModel([SoilLayer("s", c=10, phi=25, gamma=18)])
    slope = Slope.simple(height=10, ratio=2.0)
    circle = SlipCircle(10.03, 27.70, 27.69)
    slices = build_slices(slope, soil, circle, seismic=Seismic(kh=0.2, kv=0.1))
    assert slices.slide_dir == -1.0  # crest on the right, mass slides left
    assert np.all(slices.kh_force <= 0)
    assert np.allclose(slices.w_eff, slices.w * 1.1)


def test_surcharge_adds_weight_inside_its_range():
    base = build_slices(SLOPE, HOMOGENEOUS, ACADS_CIRCLE)
    loaded = build_slices(
        SLOPE,
        HOMOGENEOUS,
        ACADS_CIRCLE,
        surcharges=SurchargeSet([UniformSurcharge(q=20.0, x_start=10.0, x_end=21.0)]),
    )
    overlap = max(0.0, min(21.0, loaded.span[1]) - 10.0)
    extra = float(loaded.total_weight) - float(base.total_weight)
    assert extra == pytest.approx(20.0 * overlap, rel=0.02)


def test_surface_without_two_intersections_is_rejected():
    soil = SoilModel([SoilLayer("s", c=10, phi=25, gamma=18)])
    small_slope = Slope([(0, 0), (10, 10)])
    with pytest.raises(ValueError):
        build_slices(small_slope, soil, SlipCircle(100, 100, 5))
