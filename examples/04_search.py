"""Example 04 - search for the critical circular slip surface.

Run:  python examples/04_search.py
"""

import pathlib

import slopestab as sst

out = pathlib.Path(__file__).parent / "output"
out.mkdir(exist_ok=True)

slope = sst.Slope.simple(height=10, ratio=2.0, crest_length=20, toe_length=10)
soil = sst.SoilModel([sst.SoilLayer("silty clay", c=10.0, phi=20.0, gamma=18.0, gamma_sat=20.0)])

best = sst.search_circular(slope, soil, method="bishop", n_slices=20)
print("critical Fs :", round(best.fs, 4))
print("circle      :", best.critical.surface)
print("evaluated   :", best.n_evaluated, "surfaces in", round(best.elapsed, 2), "s")

print("\ntop 5 circles:")
for cx, cy, r, fs in best.candidates[:5]:
    print(f"  cx={cx:7.2f} cy={cy:7.2f} r={r:6.2f}  Fs={fs:.4f}")

# export the whole scan, then draw and report the governing surface
best.to_dataframe().to_csv(out / "04_candidates.csv", index=False)
sst.plot(best, out / "04_critical_surface.png")
sst.report_html(best, out / "04_search_report.html", title="最危险滑面搜索 · Bishop")
print(f"\nwrote {out}/04_candidates.csv, 04_critical_surface.png, 04_search_report.html")
