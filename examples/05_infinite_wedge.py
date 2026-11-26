"""Example 05 - infinite slope and planar wedge closed forms.

Run:  python examples/05_infinite_wedge.py
"""

import slopestab as sst

print("infinite slope, beta=26.57 deg (2H:1V), c=5 kPa, phi=25 deg, gamma=18 kN/m3")
for depth in (1.0, 2.0, 3.0, 5.0):
    dry = sst.infinite_slope_fs(26.57, depth, c=5.0, phi=25.0, gamma=18.0)
    wet = sst.infinite_slope_fs(
        26.57, depth, c=5.0, phi=25.0, gamma=18.0, gamma_sat=20.0, seepage=True
    )
    print(f"  z={depth:>4.1f} m   dry={dry:6.3f}   seepage={wet:6.3f}")

critical = sst.critical_depth(26.57, 5.0, 25.0, 18.0)
print(f"critical depth (F=1, dry): {critical:.2f} m" if critical else "never fails")

print("\nplanar wedge, plane from the toe (x=10) to the crest (x=40)")
slope = sst.Slope([(0, 0), (10, 0), (30, 10), (50, 10)])
soil = sst.SoilModel([sst.SoilLayer("weathered shale", c=15.0, phi=30.0, gamma=18.0)])
for crack, water in ((0.0, 0.0), (3.0, 0.0), (3.0, 1.5), (6.0, 3.0)):
    res = sst.planar_wedge(
        slope, soil, x_toe=10.0, x_top=40.0, crack_depth=crack, water_depth=water
    )
    print(
        f"  crack={crack:>4.1f} m  water={water:>3.1f} m  "
        f"W={res.weight:7.1f} kN/m  Fs={res.fs:6.3f}"
    )

print(
    "\nwedge with seismic kh=0.2:",
    round(
        sst.planar_wedge(
            slope, soil, x_toe=10.0, x_top=40.0, crack_depth=3.0, seismic=sst.Seismic(kh=0.2)
        ).fs,
        3,
    ),
)
