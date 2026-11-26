"""Example 02 - layered slope with a phreatic line, cross-checked by two methods.

Run:  python examples/02_layered_water.py
"""

import pathlib

import slopestab as sst

out = pathlib.Path(__file__).parent / "output"
out.mkdir(exist_ok=True)

slope = sst.Slope([(0, 0), (10, 0), (30, 10), (50, 10)])
soil = sst.SoilModel(
    [
        sst.SoilLayer(
            "fill", c=5.0, phi=25.0, gamma=18.0, gamma_sat=19.5, boundary=[(0, 3), (50, 3)]
        ),
        sst.SoilLayer(
            "clay", c=15.0, phi=18.0, gamma=19.0, gamma_sat=20.0, boundary=[(0, 7), (50, 7)]
        ),
        sst.SoilLayer("bedrock", c=50.0, phi=30.0, gamma=21.0),
    ]
)
# water table outcrops at the toe and rises under the slope
water = sst.PhreaticLine([(10, 0.0), (30, 9.0), (50, 9.0)])
circle = sst.SlipCircle(10.03, 27.70, 27.69)

for method in ("bishop", "spencer", "morgenstern_price"):
    res = sst.analyze(method, slope, soil, circle, water=water, n_slices=30)
    print(f"{method:<20} Fs = {res.fs:.4f}   lambda = {res.info.get('lambda', '-')}")

main = sst.analyze("bishop", slope, soil, circle, water=water, n_slices=30)
main.slope, main.soil, main.water, main.surface = slope, soil, water, circle

sst.plot(main, out / "02_layered_water.png")
sst.report_html(main, out / "02_layered_water.html", title="分层土坡 + 浸润线 · Bishop")
main.to_csv(out / "02_slices.csv")
print(f"\nwrote {out}/02_layered_water.png, .html and 02_slices.csv")
