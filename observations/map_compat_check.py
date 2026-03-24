#!/usr/bin/env python3
"""Quick MAP compatibility check against 13 observations."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))
from config_loader import load_config
from v22_raw_scan import sigma_T_vpm

_CFG = load_config(__file__)
_benchmarks = {b[0]: b[1:] for b in _CFG.get("benchmarks", [])}
_map = _benchmarks.get("MAP", [94.07, 11.10e-3, 5.734e-3])
mc, mp, alpha = _map[0], _map[1], _map[2]

# Read observations from config (same source as chi2_fit)
_raw_obs = _CFG.get("observations", [])
obs = [(o[0], o[1], o[2], o[3], o[4]) for o in _raw_obs]

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
