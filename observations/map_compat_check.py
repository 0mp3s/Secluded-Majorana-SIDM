#!/usr/bin/env python3
"""Quick MAP compatibility check against 13 observations."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))
from global_config import GC
from v22_raw_scan import sigma_T_vpm

_map_bp = GC.benchmark("MAP")
mc, mp, alpha = _map_bp["m_chi_GeV"], _map_bp["m_phi_MeV"] * 1e-3, _map_bp["alpha"]

# Read observations from global config
obs = [(o[0], o[1], o[2], o[3], o[4]) for o in GC.observations_as_tuples()]

print(f"MAP compatibility ({mc}, {mp*1e3:.2f} MeV, {alpha}):")
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
