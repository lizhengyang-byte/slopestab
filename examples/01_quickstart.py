"""Example 01 - analyse the ACADS 1(a) benchmark with all five methods.

Run:  python examples/01_quickstart.py
"""

import slopestab as sst

slope = sst.Slope([(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)])
soil = sst.SoilModel([sst.SoilLayer("homogeneous", c=3.0, phi=19.6, gamma=20.0)])
circle = sst.SlipCircle(0.03, 27.70, 27.69)

METHODS = ("fellenius", "bishop", "janbu", "spencer", "morgenstern_price")

print(f"{'method':<20}{'Fs':>8}{'iterations':>12}   residuals (%)")
print("-" * 60)

results = sst.compare(METHODS, slope, soil, circle, n_slices=30)
for name, res in results.items():
    usable = {k: v for k, v in res.residuals.items() if not k.endswith("_abs")}
    residuals = "  ".join(f"{k}={v:+.2f}" for k, v in usable.items())
    print(f"{name:<20}{res.fs:>8.4f}{res.iterations:>12}   {residuals}")

print("\nACADS reference answer: Fs = 1.00; Bishop / Spencer published: 0.987 / 0.986")

best = results["bishop"]
print(f"\nbishop slice table (first 5 rows): {best.slices.as_table()[:5]}")
