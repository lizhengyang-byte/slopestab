"""Command line interface: ``slopestab run`` and ``slopestab search``."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .analysis import METHOD_NAMES, analyze
from .exceptions import SlopestabError
from .io.cases import build_problem, load_case
from .search.circular import SearchGrid, search_circular
from .viz.report import report_html


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("case", help="case file (.json or .yaml)")
    parser.add_argument("--method", choices=sorted(METHOD_NAMES), help="override the method")
    parser.add_argument("--plot", metavar="PNG", help="write a section drawing")
    parser.add_argument("--report", metavar="HTML", help="write an HTML calculation report")
    parser.add_argument("--json", dest="json_out", metavar="JSON", help="write the result as JSON")
    parser.add_argument("--csv", metavar="CSV", help="write the slice table as CSV")


def _write_outputs(args, result) -> None:
    from .viz.plot import plot

    if args.plot:
        plot(result, args.plot)
    if args.report:
        report_html(result, args.report)

    if hasattr(result, "critical"):  # search result: export the candidate table
        rows = [{"cx": c[0], "cy": c[1], "r": c[2], "fs": c[3]} for c in result.candidates]
        if args.json_out:
            Path(args.json_out).write_text(
                json.dumps(
                    {"summary": result.summary(), "candidates": rows}, ensure_ascii=False, indent=2
                ),
                encoding="utf-8",
            )
        if args.csv:
            import csv as _csv

            with Path(args.csv).open("w", newline="", encoding="utf-8") as fh:
                writer = _csv.DictWriter(fh, fieldnames=["cx", "cy", "r", "fs"])
                writer.writeheader()
                writer.writerows(rows)
        return

    if args.json_out:
        result.to_json(args.json_out)
    if args.csv:
        result.to_csv(args.csv)


def cmd_run(args) -> int:
    case = load_case(args.case)
    problem = build_problem(case)
    method = args.method or problem["method"]
    result = analyze(
        method,
        problem["slope"],
        problem["soil"],
        problem["surface"],
        water=problem["water"],
        surcharges=problem["surcharges"],
        seismic=problem["seismic"],
        n_slices=problem["n_slices"],
    )
    result.slope = problem["slope"]
    result.soil = problem["soil"]
    result.surface = problem["surface"]
    result.water = problem["water"]

    print(json.dumps(result.summary(), ensure_ascii=False, indent=2))
    _write_outputs(args, result)
    return 0


def cmd_search(args) -> int:
    case = load_case(args.case)
    problem = build_problem(case)
    method = args.method or problem["method"]
    raw = problem["search"]
    grid = None
    if raw:
        grid = SearchGrid(
            cx_range=tuple(raw.get("cx_range", (0.0, 1.0))),
            cy_range=tuple(raw.get("cy_range", (0.0, 1.0))),
            r_range=tuple(raw.get("r_range", (0.0, 1.0))),
            n_cx=int(raw.get("n_cx", 24)),
            n_cy=int(raw.get("n_cy", 14)),
            n_r=int(raw.get("n_r", 12)),
        )
    result = search_circular(
        problem["slope"],
        problem["soil"],
        method=method,
        grid=grid,
        n_slices=problem["n_slices"],
        water=problem["water"],
        surcharges=problem["surcharges"],
        seismic=problem["seismic"],
    )
    print(json.dumps(result.summary(), ensure_ascii=False, indent=2))
    _write_outputs(args, result)
    return 0


def _cmd_methods(_args) -> int:
    print("\n".join(sorted(METHOD_NAMES)))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="slopestab",
        description="2D slope stability analysis by limit equilibrium methods",
    )
    parser.add_argument("--version", action="version", version=f"slopestab {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="analyse one prescribed slip surface")
    _add_common(p_run)
    p_run.set_defaults(func=cmd_run)

    p_search = sub.add_parser("search", help="search for the critical circular slip surface")
    _add_common(p_search)
    p_search.set_defaults(func=cmd_search)

    p_methods = sub.add_parser("methods", help="list available methods")
    p_methods.set_defaults(func=_cmd_methods)

    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except SlopestabError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
