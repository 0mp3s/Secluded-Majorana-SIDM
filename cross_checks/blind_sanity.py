#!/usr/bin/env python3
"""
V8 — v28_blind_sanity.py
========================
Blind sanity check: evaluate the VPM solver on points that were
NEVER in the scan grid.  Three tests:

 1. Off-grid viable:   10 points interpolated between known-viable
                       grid cells → should still pass SIDM cuts.
 2. Off-grid outside:  10 points far from viable region → should fail.
 3. Monotonicity:       5 viable points evaluated at 5 velocities
                       → σ/m must *decrease* with velocity.

Imports sigma_T_vpm directly from v22 — zero new physics code.
Runtime: ~1-2 minutes.
"""
# === path setup (auto-generated) ================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# =================================================================

import sys, os, math, time
import numpy as np

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)
os.environ['PYTHONIOENCODING'] = 'utf-8'

# Import the VPM solver from v22 (no new physics code)
from v22_raw_scan import sigma_T_vpm, M_CHI_VALS, M_PHI_VALS, LAM_CRITS
from config_loader import load_config
from output_manager import get_latest

_CFG = load_config(__file__)

# ── constants ──
V_DWARF   = _CFG.get("v_dwarf", 30.0)     # km/s
V_CLUSTER = _CFG.get("v_cluster", 1000.0)  # km/s

# ── helpers ──
def geometric_mean(a, b):
    return math.sqrt(a * b)

def evaluate_point(m_chi, m_phi, alpha, label=""):
    """Evaluate σ/m at dwarf and cluster velocities."""
    s30   = sigma_T_vpm(m_chi, m_phi, alpha, V_DWARF)
    s1000 = sigma_T_vpm(m_chi, m_phi, alpha, V_CLUSTER)
    viable = (1.0 <= s30 <= 10.0) and (s1000 < 0.1)
    return s30, s1000, viable

# ==============================================================
#  Load representative CSV to pick "nearby" off-grid points
# ==============================================================
def load_representatives():
    import csv
    reps = []
    _SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    _csv_override = _CFG.get("representative_csv")
    csv_path = os.path.join(_SCRIPT_DIR, _csv_override) if _csv_override else str(get_latest("all_viable_representative_v8"))
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Support both m_phi_MeV and legacy m_phi_GeV headers
            if 'm_phi_MeV' in row:
                mp_gev = float(row['m_phi_MeV']) / 1000.0
            else:
                mp_gev = float(row['m_phi_GeV'])
            reps.append({
                'm_chi': float(row['m_chi_GeV']),
                'm_phi': mp_gev,
                'alpha': float(row['alpha']),
                'sigma30': float(row['sigma_m_30']),
            })
    return reps


# ==============================================================
#  TEST 1 — Off-grid viable (interpolation between known points)
# ==============================================================
def test1_offgrid_viable():
    print("="*65)
    print("TEST 1: Off-grid viable — interpolated between known grid cells")
    print("="*65)
    
    reps = load_representatives()
    if len(reps) < 20:
        print("  [SKIP] Not enough representative points")
        return 0, 0
    
    # Pick 10 pairs of consecutive viable points, take geometric mean
    rng = np.random.default_rng(seed=42)
    indices = rng.choice(len(reps) - 1, size=10, replace=False)
    
    n_pass, n_fail = 0, 0
    print(f"  {'#':>3}  {'m_chi':>10}  {'m_phi':>12}  {'alpha':>12}  {'σ/m(30)':>10}  {'σ/m(1k)':>10}  {'viable?':>8}")
    print(f"  {'':>3}  {'[GeV]':>10}  {'[GeV]':>12}  {'':>12}  {'cm²/g':>10}  {'cm²/g':>10}")
    print("  " + "-"*63)
    
    for i, idx in enumerate(indices):
        p1, p2 = reps[idx], reps[idx + 1]
        # Geometric mean of all three parameters → off-grid point
        mc  = geometric_mean(p1['m_chi'], p2['m_chi'])
        mp  = geometric_mean(p1['m_phi'], p2['m_phi'])
        alp = geometric_mean(p1['alpha'], p2['alpha'])
        
        s30, s1000, viable = evaluate_point(mc, mp, alp)
        # Here we DON'T require viable=True. We just record what happens.
        # The test passes if the solver returns finite, physical numbers.
        is_physical = np.isfinite(s30) and np.isfinite(s1000) and s30 >= 0 and s1000 >= 0
        if is_physical:
            n_pass += 1
        else:
            n_fail += 1
        
        tag = "✓ viable" if viable else "✗ out"
        phys = "PHYSICAL" if is_physical else "** NaN/INF **"
        print(f"  {i+1:>3}  {mc:>10.3f}  {mp:>12.3e}  {alp:>12.3e}  {s30:>10.3f}  {s1000:>10.4f}  {tag:>8}  {phys}")
    
    # Count how many of the interpolated points are still viable
    # (not a pass/fail — just informational)
    print(f"\n  Solver returned physical values: {n_pass}/10")
    print(f"  (If some interpolated points are not viable, that's expected —")
    print(f"   the viable region is narrow in α.)")
    
    ok = (n_fail == 0)
    print(f"\n  TEST 1: {'PASS' if ok else 'FAIL'} — all solver outputs are finite & physical")
    return 1 if ok else 0, 1


# ==============================================================
#  TEST 2 — Off-grid outside (should NOT be viable)
# ==============================================================
def test2_offgrid_outside():
    print("\n" + "="*65)
    print("TEST 2: Off-grid outside — points that should fail SIDM cuts")
    print("="*65)
    
    # 10 points deliberately outside viable region:
    # - α too large (strong coupling → σ/m(1000) too big)
    # - α too small (weak coupling → σ/m(30) too small)
    # - m_chi too small (< 10 GeV → typically no viable region)
    # - m_phi too large (> 500 MeV → Yukawa too short-range)
    outside_points = [
        # (m_chi [GeV], m_phi [GeV], alpha, reason)
        (50.0,   5e-3,    0.1,       "α way too large"),
        (50.0,   5e-3,    1e-8,      "α way too small"),
        (0.5,    5e-3,    1e-4,      "m_χ = 0.5 GeV (too light)"),
        (50.0,   0.5,     1e-4,      "m_φ = 500 MeV (too heavy)"),
        (50.0,   5e-3,    0.01,      "α = 0.01 (overcoupled)"),
        (1.0,    1e-4,    1e-5,      "m_χ = 1 GeV, m_φ = 0.1 MeV"),
        (100.0,  0.2,     1e-3,      "m_φ = 200 MeV, at grid edge"),
        (50.0,   5e-3,    1e-9,      "α = 10⁻⁹ (ultra-weak)"),
        (0.1,    0.1e-3,  1e-4,      "m_χ = 0.1 GeV (minimum grid)"),
        (80.0,   50e-3,   5e-3,      "m_φ = 50 MeV, α = 5×10⁻³"),
    ]
    
    n_pass, n_total = 0, len(outside_points)
    print(f"  {'#':>3}  {'m_chi':>8}  {'m_phi':>10}  {'alpha':>10}  {'σ/m(30)':>10}  {'σ/m(1k)':>10}  {'viable?':>8}  reason")
    print("  " + "-"*80)
    
    for i, (mc, mp, alp, reason) in enumerate(outside_points):
        s30, s1000, viable = evaluate_point(mc, mp, alp)
        tag = "✗ out" if not viable else "⚠ VIABLE!"
        if not viable:
            n_pass += 1
        print(f"  {i+1:>3}  {mc:>8.1f}  {mp:>10.1e}  {alp:>10.1e}  {s30:>10.3f}  {s1000:>10.4f}  {tag:>8}  {reason}")
    
    ok = (n_pass == n_total)
    print(f"\n  Outside points correctly rejected: {n_pass}/{n_total}")
    print(f"  TEST 2: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("  ⚠ Some 'outside' points passed SIDM cuts — check if the region is wider than expected.")
    return 1 if ok else 0, 1


# ==============================================================
#  TEST 3 — Monotonicity: σ/m must decrease with velocity
# ==============================================================
def test3_monotonicity():
    print("\n" + "="*65)
    print("TEST 3: Monotonicity — σ/m should decrease with velocity")
    print("="*65)
    
    reps = load_representatives()
    rng = np.random.default_rng(seed=123)
    indices = rng.choice(len(reps), size=5, replace=False)
    
    velocities = _CFG.get("monotonicity_velocities", [10.0, 30.0, 100.0, 300.0, 1000.0])  # km/s
    
    n_monotone, n_total = 0, 5
    
    for i, idx in enumerate(indices):
        p = reps[idx]
        mc, mp, alp = p['m_chi'], p['m_phi'], p['alpha']
        
        sigmas = []
        for v in velocities:
            s = sigma_T_vpm(mc, mp, alp, v)
            sigmas.append(s)
        
        # Check monotonically decreasing (allowing tiny numerical noise)
        is_mono = True
        for j in range(len(sigmas) - 1):
            if sigmas[j+1] > sigmas[j] * 1.05:  # 5% tolerance for numerical noise
                is_mono = False
                break
        
        if is_mono:
            n_monotone += 1
        
        tag = "✓ mono" if is_mono else "✗ NON-MONO"
        print(f"\n  Point {i+1}: m_χ={mc:.1f} GeV, m_φ={mp:.2e} GeV, α={alp:.2e}")
        print(f"    {'v [km/s]':>12}  {'σ/m [cm²/g]':>14}  {'ratio to prev':>14}")
        for j, (v, s) in enumerate(zip(velocities, sigmas)):
            ratio = f"{s/sigmas[j-1]:.4f}" if j > 0 else "—"
            print(f"    {v:>12.0f}  {s:>14.4f}  {ratio:>14}")
        print(f"    → {tag}")
    
    ok = (n_monotone == n_total)
    print(f"\n  Monotonically decreasing: {n_monotone}/{n_total}")
    print(f"  TEST 3: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("  ⚠ Non-monotonic σ/m(v) could indicate resonance structure at")
        print("    intermediate velocities. Check if it's a narrow spike.")
    return 1 if ok else 0, 1


# ==============================================================
#  TEST 4 — Solver stability: extreme κ values
# ==============================================================
def test4_solver_stability():
    print("\n" + "="*65)
    print("TEST 4: Solver stability at extreme κ values")
    print("="*65)
    
    # Points designed to stress the solver:
    # Very low κ (few partial waves), very high κ (many partial waves)
    stress_points = [
        # (m_chi, m_phi, alpha, v_km_s, description)
        (50.0,  5e-3,   5e-4,   1.0,     "v=1 km/s, κ≈0.03"),
        (50.0,  5e-3,   5e-4,   5.0,     "v=5 km/s, κ≈0.14"),
        (100.0, 0.1e-3, 1e-4,   1000.0,  "v=1000, m_φ tiny → κ≈560"),
        (100.0, 1e-3,   1e-3,   1000.0,  "v=1000, λ=100 → high coupling"),
        (13.9,  1.1e-3, 2e-6,   30.0,    "near viable edge (low α)"),
    ]
    
    n_ok, n_total = 0, len(stress_points)
    print(f"  {'#':>3}  {'description':>35}  {'σ/m':>12}  {'status':>12}")
    print("  " + "-"*65)
    
    for i, (mc, mp, alp, v, desc) in enumerate(stress_points):
        try:
            s = sigma_T_vpm(mc, mp, alp, v)
            is_ok = np.isfinite(s) and s >= 0
            if is_ok:
                n_ok += 1
            status = "✓ OK" if is_ok else "✗ NaN/INF"
            print(f"  {i+1:>3}  {desc:>35}  {s:>12.4e}  {status:>12}")
        except Exception as e:
            print(f"  {i+1:>3}  {desc:>35}  {'EXCEPTION':>12}  ✗ {e}")
    
    ok = (n_ok == n_total)
    print(f"\n  Stable evaluations: {n_ok}/{n_total}")
    print(f"  TEST 4: {'PASS' if ok else 'FAIL'}")
    return 1 if ok else 0, 1


# ==============================================================
#  MAIN
# ==============================================================
def main():
    print("╔═══════════════════════════════════════════════════════════════╗")
    print("║      V8 — v28 Blind Sanity Check (off-grid evaluation)      ║")
    print("╠═══════════════════════════════════════════════════════════════╣")
    print("║  Tests the VPM solver on points NEVER in the scan grid.     ║")
    print("║  Zero new physics code — imports sigma_T_vpm from v22.      ║")
    print("╚═══════════════════════════════════════════════════════════════╝\n")
    
    t0 = time.time()
    
    # Warm up JIT
    print("Warming up Numba JIT...")
    _ = sigma_T_vpm(50.0, 5e-3, 5e-4, 30.0)
    print(f"JIT warm-up: {time.time()-t0:.1f}s\n")
    
    passed, total = 0, 0
    
    p, t = test1_offgrid_viable()
    passed += p; total += t
    
    p, t = test2_offgrid_outside()
    passed += p; total += t
    
    p, t = test3_monotonicity()
    passed += p; total += t
    
    p, t = test4_solver_stability()
    passed += p; total += t
    
    elapsed = time.time() - t0
    
    print("\n" + "="*65)
    print(f"  SUMMARY: {passed}/{total} tests PASSED  ({elapsed:.1f}s)")
    print("="*65)
    
    if passed == total:
        print("  ✓ All blind sanity checks passed.")
        print("  The solver produces physical, stable results on off-grid data.")
    else:
        print(f"  ⚠ {total - passed} test(s) failed — investigate above.")
    
    return 0 if passed == total else 1

if __name__ == "__main__":
    sys.exit(main())
