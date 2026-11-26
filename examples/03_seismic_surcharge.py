"""Example 03 - seismic and surcharge screening for the same section.

Run:  python examples/03_seismic_surcharge.py
"""

import slopestab as sst

slope = sst.Slope.simple(height=10, ratio=2.0, crest_length=20, toe_length=10)
soil = sst.SoilModel([sst.SoilLayer("silty sand", c=8.0, phi=26.0, gamma=18.5, gamma_sat=20.0)])
circle = sst.SlipCircle(10.03, 27.70, 27.69)

cases = {
    "static": {},
    "kh = 0.10": {"seismic": sst.Seismic(kh=0.10)},
    "kh = 0.20": {"seismic": sst.Seismic(kh=0.20)},
    "kh = 0.20, kv = 0.10": {"seismic": sst.Seismic(kh=0.20, kv=0.10)},
    "q = 30 kPa on crest": {"surcharge": sst.UniformSurcharge(q=30.0, x_start=30.0, x_end=60.0)},
    "kh = 0.15 + q = 30 kPa": {
        "seismic": sst.Seismic(kh=0.15),
        "surcharge": sst.UniformSurcharge(q=30.0, x_start=30.0, x_end=60.0),
    },
}

print(f"{'case':<28}{'Fs':>8}")
print("-" * 38)
for label, kwargs in cases.items():
    res = sst.analyze("spencer", slope, soil, circle, **kwargs)
    print(f"{label:<28}{res.fs:>8.4f}")

print("\nNote: kh acts in the sliding direction; kv is positive downwards (destabilising).")
