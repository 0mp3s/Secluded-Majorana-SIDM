#!/usr/bin/env python3
"""Quick MAP compatibility check against 13 observations."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))
from v22_raw_scan import sigma_T_vpm

mc, mp, alpha = 94.07, 11.10e-3, 5.734e-3

obs = [
    ("Draco dSph",          12,   0.6,  0.1,  2.0),
    ("Fornax dSph",         12,   0.8,  0.2,  3.0),
    ("NGC 2976",            60,   2.0,  0.5,  5.0),
    ("NGC 1560",            55,   3.0,  1.0,  8.0),
    ("IC 2574",             50,   1.5,  0.3,  5.0),
    ("NGC 720 (group)",    250,   0.5,  0.1,  1.5),
    ("NGC 1332 (group)",   280,   0.3,  0.05, 1.0),
    ("Abell 611",         1200,   0.1,  0.02, 0.3),
    ("Abell 2537",        1100,   0.15, 0.03, 0.4),
    ("Diverse RC band",     40,   3.0,  0.5,  10.0),
    ("Bullet Cluster",    4700,   0.7,  0.0,  1.25),
    ("72 cluster mergers", 1000,  0.2,  0.0,  0.47),
    ("TBTF dwarfs",         30,   1.0,  0.5,  5.0),
]

print("MAP compatibility (94.07, 11.10 MeV, 5.734e-3):")
print("-" * 80)
n = 0
for name, v, sm, lo, hi in obs:
    t = sigma_T_vpm(mc, mp, alpha, float(v))
    ok = lo <= t <= hi
    if ok:
        n += 1
        s = "OK"
    elif t < lo:
        s = f"BELOW ({t:.3f} < {lo:.2f})"
    else:
        s = f"ABOVE ({t:.3f} > {hi:.2f})"
    print(f"  {name:<25} v={v:5.0f}  theory={t:.4f}  [{lo:.2f}, {hi:.2f}]  {s}")

print(f"\nMAP compatible: {n}/{len(obs)}")

# Key velocities summary
print("\nKey sigma/m values:")
for v in [12, 30, 40, 50, 100, 220, 1000, 1200]:
    s = sigma_T_vpm(mc, mp, alpha, float(v))
    print(f"  sigma/m({v}) = {s:.4f} cm^2/g")
