"""I/O, reporting and CLI tests."""

import json

import pytest

import slopestab as sst
from slopestab.cli import main

CASE = {
    "method": "bishop",
    "n_slices": 24,
    "slope": {"profile": [[0, 0], [10, 0], [30, 10], [50, 10]]},
    "soil": {
        "layers": [
            {"name": "fill", "c": 5.0, "phi": 25.0, "gamma": 18.0, "gamma_sat": 19.5},
        ]
    },
    "surface": {"type": "circle", "cx": 10.03, "cy": 27.70, "r": 27.69},
    "water": {"type": "phreatic", "points": [[10, 0], [30, 9], [50, 9]]},
    "seismic": {"kh": 0.05, "kv": 0.0},
    "surcharge": [{"q": 10.0, "x_start": 40.0, "x_end": 50.0}],
}


def test_case_roundtrip(tmp_path):
    path = tmp_path / "case.json"
    sst.dump_case(CASE, path)
    loaded = sst.load_case(path)
    assert loaded == json.loads(json.dumps(CASE))
    problem = sst.build_problem(loaded)
    res = sst.analyze(
        problem["method"],
        problem["slope"],
        problem["soil"],
        problem["surface"],
        water=problem["water"],
        surcharges=problem["surcharges"],
        seismic=problem["seismic"],
        n_slices=problem["n_slices"],
    )
    assert 0.0 < res.fs < 20.0


def test_case_with_ru_and_surcharge():
    case = dict(CASE)
    case["water"] = {"type": "ru", "ru": 0.3}
    problem = sst.build_problem(case)
    assert isinstance(problem["water"], sst.RuCoefficient)
    assert len(problem["surcharges"]) == 1


def test_result_exports(tmp_path):
    case = dict(CASE)
    case["water"] = None
    problem = sst.build_problem(case)
    res = sst.analyze("bishop", problem["slope"], problem["soil"], problem["surface"])
    res.slope = problem["slope"]
    res.soil = problem["soil"]
    res.surface = problem["surface"]

    csv_path = res.to_csv(tmp_path / "slices.csv")
    assert csv_path.exists() and "w" in csv_path.read_text(encoding="utf-8").splitlines()[0]

    json_path = res.to_json(tmp_path / "result.json")
    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["summary"]["method"] == "bishop"
    assert len(data["slices"]) == len(res.slices)


def test_html_report_is_self_contained(tmp_path):
    case = dict(CASE)
    case["water"] = None
    problem = sst.build_problem(case)
    res = sst.analyze("bishop", problem["slope"], problem["soil"], problem["surface"])
    res.slope, res.soil, res.surface = problem["slope"], problem["soil"], problem["surface"]

    html_path = tmp_path / "report.html"
    sst.report_html(res, html_path)
    text = html_path.read_text(encoding="utf-8")
    assert "边坡稳定分析计算书" in text
    assert "Fs" in text
    assert "data:image/png;base64" in text  # figure embedded, no external assets


def test_plot_writes_a_png(tmp_path):
    case = dict(CASE)
    case["water"] = None
    problem = sst.build_problem(case)
    res = sst.analyze("spencer", problem["slope"], problem["soil"], problem["surface"])
    res.slope, res.soil, res.surface = problem["slope"], problem["soil"], problem["surface"]

    png = tmp_path / "section.png"
    sst.plot(res, png)
    assert png.exists() and png.stat().st_size > 10_000


def test_cli_run_and_search(tmp_path, capsys):
    case_path = tmp_path / "case.json"
    sst.dump_case(CASE, case_path)

    out_json = tmp_path / "out.json"
    code = main(["run", str(case_path), "--json", str(out_json)])
    assert code == 0
    assert json.loads(out_json.read_text(encoding="utf-8"))["summary"]["fs"] > 0

    out_csv = tmp_path / "search.csv"
    code = main(["search", str(case_path), "--method", "bishop", "--csv", str(out_csv)])
    assert code == 0
    assert "cx" in out_csv.read_text(encoding="utf-8").splitlines()[0]


def test_bad_case_file_returns_error_code(tmp_path, capsys):
    bad = tmp_path / "bad.json"
    bad.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError):
        main(["run", str(bad)])
