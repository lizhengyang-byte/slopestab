# Changelog

All notable changes to this project are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-11-26

### Added

- Geometry kernel: piecewise-linear slope profiles, multi-layer soil models with
  polyline layer boundaries, circular and polyline slip surfaces.
- Slice builder with forced cuts at layer boundaries, water-table crossings and
  surcharge edges; piecewise weight integration above/below the water table.
- Limit-equilibrium methods: Fellenius (ordinary), Bishop simplified, Janbu
  simplified (with user-supplied `f0` correction), Spencer, Morgenstern-Price
  (constant / half-sine / sine / trapezoid interslice functions).
- Global equilibrium audit for every result: force and moment imbalance reported
  as a percentage of the total slice weight.
- Pore-pressure models: dry, `r_u` coefficient, phreatic line (hydrostatic).
- Pseudo-static seismic loading (`kh`, downward-positive `kv`) and uniform
  surcharges.
- Circular slip-surface search: three-dimensional centre/radius grid with
  geometric pre-filters and Nelder-Mead polishing.
- Special cases: infinite-slope closed form (dry and seepage parallel to the
  slope) and planar wedge with tension crack and crack water.
- Outputs: section drawing (matplotlib optional), self-contained HTML report,
  CSV/JSON slice tables, YAML/JSON case files, `slopestab run` / `slopestab search`
  CLI.
- Verification suite: ACADS 1(a) benchmark, undrained moment identity,
  closed-form infinite slope, cross-method agreement, CI on Python 3.9/3.11/3.13.
- A referee test that recomputes Bishop's simplified method from scratch (own
  circle/ground intersections, own equal-width slices, own `m_alpha` iteration)
  and must agree with the library to 0.5%; plus hand-formula checks for the
  planar wedge and a mass-conservation check against a shoelace polygon.
- Bilingual documentation: `README.md` (Chinese) and `README.en.md` (English),
  kept structurally parallel.

### Fixed

- Pore water was missing from the generalised (Spencer / Morgenstern-Price)
  per-slice equilibrium, which inflated the rigorous factor of safety by up to
  22% on saturated profiles (layered sand/clay example: Spencer 3.105 -> 2.554).
  The pore resultant now enters both the normal force and the force/moment
  closure; `normal_forces(..., lambda=0)` reproduces Bishop's closed form exactly.
- Morgenstern-Price no longer crashes when a trial `lambda` leaves the solvable
  domain: the scan drops undefined points, the root finder is a secant/bisection
  hybrid, and a `lambda` with a residual force/moment gap <= 0.02 is accepted
  with `converged=False` and an explanatory `info["note"]` instead of a spurious
  error or a silently wrong number.
- Section drawing: the slide mass is now filled as a polygon between the slip
  surface and the ground, and the legend no longer repeats the title.
- A phreatic line above the ground surface used to run silently: pore pressure
  exceeded the slice weight, effective normal forces went negative and Bishop
  returned a meaningless `Fs = 0.46` with `converged=True`.  The slice builder
  now rejects that input with a `ValueError` naming the offending abscissa and
  quoting `u*b` versus `W`.

### Known limitations

- Rigorous GLE (Spencer, Morgenstern-Price) is restricted to circular surfaces;
  non-circular rigorous analysis is planned for v0.2.
- Janbu's correction factor `f0` is not auto-computed from the chart; the
  uncorrected value is returned unless the user supplies `f0`.
- No finite-element strength reduction, no 3D wedge analysis, no seepage flow
  solver.
