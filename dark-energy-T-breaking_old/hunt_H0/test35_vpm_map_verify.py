#!/usr/bin/env python3
"""
Test 35: Independent VPM verification of MAP point
===================================================
MAP = (m_χ = 94.07 GeV, m_φ = 11.10 MeV, α = 5.734×10⁻³)

Verifies σ_T/m at all SIDM velocity windows using the
full VPM solver (identical Majorana fermions, spin-statistical weights).

NOTE: old/_vpm_map_check.py had a UNITS BUG: used m_chi = 94.07e-3 (MeV)
instead of 94.07 (GeV). This script uses correct units.
"""
import sys, os, math, time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from consistency_check_sidm import sigma_T_vpm

# MAP parameters — ALL IN GeV
M_CHI = 94.07        # GeV  (NOT 94.07e-3!)
M_PHI = 11.10e-3     # GeV  (= 11.10 MeV)
ALPHA = 5.734e-3     # dimensionless

# Derived scales
LAM = ALPHA * M_CHI / M_PHI          # λ = α m_χ / m_φ
V_MED = (M_PHI / M_CHI) * 299792.458  # km/s  (mediator velocity scale)

# Full SIDM constraint set (13 windows from Paper 1)
CONSTRAINTS = [
    # (label, v [km/s], σ_m_lo, σ_m_hi)
    ("Fornax dSph",        12,   0.5,  100.0),
    ("Draco/Sculptor",     30,   0.5,   50.0),
    ("LSB galaxies",       60,   0.5,   10.0),
    ("MW substructure",   100,   0.0,   35.0),
    ("MW satellites",     200,   0.0,   35.0),
    ("Milky Way",         220,   0.0,   10.0),
    ("Galaxy groups",     500,   0.1,    3.0),
    ("Galaxy clusters",  1000,   0.1,    1.0),
    ("Abell 3827",       1500,   0.0,    1.5),
    ("MACS J0025",       2000,   0.0,    3.0),
    ("Bullet Cluster",   3000,   0.0,    1.25),
    ("El Gordo",         4000,   0.0,    2.0),
    ("Musket Ball",      4500,   0.0,    7.0),
]

def main():
    print("=" * 90)
    print("  TEST 35: Independent VPM Verification of MAP Point")
    print("=" * 90)
    print()
    print(f"  m_χ   = {M_CHI} GeV")
    print(f"  m_φ   = {M_PHI*1e3} MeV")
    print(f"  α     = {ALPHA}")
    print(f"  λ     = α m_χ/m_φ = {LAM:.3f}")
    print(f"  v_med = m_φ/m_χ × c = {V_MED:.1f} km/s")
    print()

    # JIT warmup
    print("  [JIT warmup]...", end="", flush=True)
    _ = sigma_T_vpm(20.0, 0.01, 1e-3, 30.0)
    print(" done.")
    print()

    # Run VPM at all velocity windows
    print(f"  {'Constraint':>20}  {'v [km/s]':>10}  {'σ/m [cm²/g]':>14}  "
          f"{'lo':>6}  {'hi':>6}  {'Result':>8}")
    print("  " + "-" * 75)

    t0 = time.time()
    passed = 0
    results = []
    for label, v, lo, hi in CONSTRAINTS:
        sm = sigma_T_vpm(M_CHI, M_PHI, ALPHA, float(v))
        ok = lo <= sm <= hi
        if ok:
            passed += 1
        tag = "PASS" if ok else "FAIL"
        results.append((label, v, sm, lo, hi, ok))
        print(f"  {label:>20}  {v:>10}  {sm:>14.4f}  "
              f"{lo:>6.2f}  {hi:>6.2f}  [{tag:>4}]")

    dt = time.time() - t0
    print()
    print(f"  Time: {dt:.1f}s")
    print(f"  Result: {passed}/{len(CONSTRAINTS)} constraints passed")
    print()

    # Also compute relic-matched MAP (α=4.523e-3)
    ALPHA_RELIC = 4.523e-3
    print("  --- MAP_relic (α = 4.523e-3, relic-matched) ---")
    print(f"  {'Constraint':>20}  {'v [km/s]':>10}  {'σ/m [cm²/g]':>14}  {'Result':>8}")
    print("  " + "-" * 60)
    passed_r = 0
    for label, v, lo, hi in CONSTRAINTS:
        sm = sigma_T_vpm(M_CHI, M_PHI, ALPHA_RELIC, float(v))
        ok = lo <= sm <= hi
        if ok:
            passed_r += 1
        tag = "PASS" if ok else "FAIL"
        print(f"  {label:>20}  {v:>10}  {sm:>14.4f}  [{tag:>4}]")

    print()
    print(f"  MAP_relic: {passed_r}/{len(CONSTRAINTS)} constraints passed")
    print()

    # Conclusion
    if passed == len(CONSTRAINTS):
        print("  ✅ CONCLUSION: MAP point is FULLY SIDM-consistent (all 13 windows).")
    elif passed >= 10:
        print(f"  ⚠️  CONCLUSION: MAP passes {passed}/13 — mostly consistent.")
    else:
        print(f"  ❌ CONCLUSION: MAP fails {len(CONSTRAINTS)-passed}/13 constraints.")
    print()

    # Bug note
    print("  HISTORICAL NOTE: old/_vpm_map_check.py used m_chi = 94.07e-3 (MeV)")
    print("  instead of 94.07 (GeV) — a factor-1000 mass error that gave σ/m ~ 36000.")
    print("  This script uses the correct m_chi = 94.07 GeV.")

if __name__ == "__main__":
    main()
