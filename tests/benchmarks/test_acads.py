"""Benchmarks against published verification values.

The reference numbers below come from:

* ACADS 1(a) - Giam & Donald (1989), the survey answer is Fs = 1.00;
* Slide2 verification manual (Rocscience), ACADS 1(a) case:
  Bishop 0.987, Spencer 0.986, GLE 0.986, Janbu corrected 0.990;
* Zhao (2014) reporting the same problem: Bishop 0.992, M-P 0.991.

The prescribed circle is the one published with the ACADS solution
(centre 0.03, 27.70; radius 27.69; area 47.17 m2).
"""

import pytest

import slopestab as sst
from slopestab import SlipCircle, Slope, SoilLayer, SoilModel

ACADS_SLOPE = Slope([(0, 0), (20, 10), (30, 10)])
ACADS_SOIL = SoilModel([SoilLayer("homogeneous", c=3.0, phi=19.6, gamma=20.0)])
ACADS_CIRCLE = SlipCircle(0.03, 27.70, 27.69)

PUBLISHED = {
    "bishop": 0.987,
    "spencer": 0.986,
    "morgenstern_price": 0.986,
}


@pytest.mark.slow
@pytest.mark.parametrize("method,expected", list(PUBLISHED.items()))
def test_acads_1a_against_the_published_value(method, expected):
    res = sst.analyze(method, ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=30)
    assert res.fs == pytest.approx(expected, abs=0.02), f"{method}: got {res.fs:.4f}"


@pytest.mark.slow
def test_acads_1a_survey_answer_is_reproduced_within_five_percent():
    res = sst.analyze("morgenstern_price", ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=30)
    assert abs(res.fs - 1.00) / 1.00 < 0.05


@pytest.mark.slow
def test_mesh_refinement_converges_to_a_stable_value():
    values = [
        sst.analyze("bishop", ACADS_SLOPE, ACADS_SOIL, ACADS_CIRCLE, n_slices=n).fs
        for n in (10, 20, 40, 80)
    ]
    assert max(values) - min(values) < 0.01
