"""Infinite slope: closed-form solutions used as validation anchors.

Two classical cases are covered:

* dry slope, failure plane at depth ``z`` parallel to the surface;
* steady seepage flowing parallel to the slope, water table at the surface.

Both reduce to hand-checkable formulas, which is why they are part of the
test suite: if the closed form and the method-of-slices machinery disagree,
something deeper is wrong than a rounding error.
"""

from __future__ import annotations

import math

GAMMA_W = 9.81


def infinite_slope_fs(
    beta_deg: float,
    z: float,
    c: float = 0.0,
    phi: float = 30.0,
    gamma: float = 18.0,
    gamma_sat: float | None = None,
    gamma_w: float = GAMMA_W,
    seepage: bool = False,
) -> float:
    """Factor of safety of an infinite slope with a plane failure surface.

    Parameters
    ----------
    beta_deg:
        Slope angle (degrees).
    z:
        Vertical depth of the failure plane below the surface (m).
    c, phi:
        Effective cohesion (kPa) and friction angle (degrees).
    gamma, gamma_sat:
        Bulk and saturated unit weights (kN/m3).
    seepage:
        When true, steady seepage parallel to the slope with the water table at
        the ground surface is assumed.

    Notes
    -----
    Dry case::

        F = [c' + gamma*z*cos^2(beta)*tan(phi')] / [gamma*z*sin(beta)*cos(beta)]

    Seepage case (flow parallel to the surface)::

        F = [c' + (gamma_sat - gamma_w)*z*cos^2(beta)*tan(phi')] /
            [gamma_sat*z*sin(beta)*cos(beta)]
    """
    if not 0.0 < beta_deg < 90.0:
        raise ValueError("beta_deg must be in (0, 90)")
    if z <= 0.0:
        raise ValueError("z must be positive")
    beta = math.radians(beta_deg)
    gamma_sat = gamma if gamma_sat is None else gamma_sat
    tan_phi = math.tan(math.radians(phi))

    if seepage:
        normal_gamma = gamma_sat - gamma_w
        driving_gamma = gamma_sat
    else:
        normal_gamma = gamma
        driving_gamma = gamma

    numerator = c + normal_gamma * z * math.cos(beta) ** 2 * tan_phi
    denominator = driving_gamma * z * math.sin(beta) * math.cos(beta)
    return numerator / denominator


def critical_depth(
    beta_deg: float,
    c: float,
    phi: float,
    gamma: float,
    gamma_sat: float | None = None,
    gamma_w: float = GAMMA_W,
    seepage: bool = False,
) -> float | None:
    """Depth at which an infinite slope reaches ``F = 1``.

    Returns ``None`` when the slope cannot fail at any depth (the friction
    angle is steep enough relative to the slope angle).
    """
    beta = math.radians(beta_deg)
    gamma_sat = gamma if gamma_sat is None else gamma_sat
    tan_phi = math.tan(math.radians(phi))
    normal_gamma = (gamma_sat - gamma_w) if seepage else gamma
    driving_gamma = gamma_sat if seepage else gamma
    denom = (
        driving_gamma
        * math.cos(beta)
        * (math.sin(beta) - (normal_gamma / driving_gamma) * math.cos(beta) * tan_phi)
    )
    if denom <= 0.0:
        return None
    return c / denom
