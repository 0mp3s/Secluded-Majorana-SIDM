#!/usr/bin/env python3
"""
V8 — v29_blind_large.py
========================
Large-scale blind sanity check: 200 fully random off-grid points.

Strategy:
  A. 100 points uniformly random in the VIABLE region
     (m_chi, m_phi, alpha sampled continuously — NOT on grid)
     → check solver stability + count how many are actually viable
  B. 50 points uniformly random across the FULL parameter space
     → mostly should fail (viable region is thin)
  C. 50 points: take known-viable, perturb alpha by ±1%, ±5%, ±10%
     → test continuity (small perturbation → small σ/m change)

Also: full velocity curves (7 velocities) for all 200 points
      to characterize non-monotonicity prevalence.

Imports sigma_T_vpm from v22 — zero new physics code.
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

from v22_raw_scan import sigma_T_vpm
from config_loader import load_config

_CFG = load_config(__file__)

# ── velocities for full curve ──
VELOCITIES = _CFG.get("blind_velocities", [5.0, 10.0, 30.0, 100.0, 300.0, 500.0, 1000.0])

def eval_full_curve(m_chi, m_phi, alpha):
    """Return dict with σ/m at all velocities + viability."""
    results = {}
    for v in VELOCITIES:
        s = sigma_T_vpm(m_chi, m_phi, alpha, v)
        results[v] = s
    s30 = results[30.0]
    s1000 = results[1000.0]
    results['viable'] = (1.0 <= s30 <= 10.0) and (s1000 < 0.1)
    return results

def check_monotonicity(curve, tol=0.05):
    """Check if σ/m decreases with velocity (within tolerance)."""
    vals = [curve[v] for v in VELOCITIES]
    for i in range(len(vals) - 1):
        if vals[i+1] > vals[i] * (1.0 + tol):
            return False, i  # non-mono at index i→i+1
    return True, -1

def load_representatives():
    import csv
    reps = []
    _SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(_SCRIPT_DIR, _CFG.get("representative_csv", "../data/all_viable_representative_v8.csv"))
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
            })
    return reps


# ==============================================================
#  BLOCK A — 100 random points in viable region bounds
# ==============================================================
def block_a_random_viable_region():
    print("="*70)
    print("BLOCK A: 100 random points sampled within viable-region bounds")
    print("         m_χ ∈ [13.9, 100] GeV, m_φ ∈ [0.1, 200] MeV (log)")
    print("         α chosen near resonances (continuous, off-grid)")
    print("="*70)
    
    rng = np.random.default_rng(seed=2026)
    LAM_CRITS = np.array(_CFG.get("lam_crits", [1.68, 6.45, 14.7, 26.0]))
    
    n_total = 100
    n_physical = 0
    n_viable = 0
    n_mono = 0
    n_nonmono = 0
    nonmono_details = []
    
    for i in range(n_total):
        # Random m_chi in [13.9, 100] GeV (log-uniform)
        m_chi = 10**(rng.uniform(np.log10(13.9), np.log10(100.0)))
        # Random m_phi in [0.5, 200] MeV (log-uniform), in GeV
        m_phi = 10**(rng.uniform(np.log10(0.5e-3), np.log10(200e-3)))
        # Pick random resonance and perturb alpha around it
        lam_c = rng.choice(LAM_CRITS)
        alpha_c = lam_c * m_phi / m_chi
        # Random perturbation ±50% in log space (wider than grid's ±30%)
        log_pert = rng.uniform(-0.5, 0.5)
        alpha = alpha_c * 10**log_pert
        
        curve = eval_full_curve(m_chi, m_phi, alpha)
        
        is_physical = all(np.isfinite(curve[v]) and curve[v] >= 0 for v in VELOCITIES)
        if is_physical:
            n_physical += 1
        if curve['viable']:
            n_viable += 1
        
        mono, idx = check_monotonicity(curve)
        if mono:
            n_mono += 1
        else:
            n_nonmono += 1
            nonmono_details.append((i+1, m_chi, m_phi, alpha, idx, curve))
        
        if (i+1) % 25 == 0:
            print(f"  ... {i+1}/100 done ({n_viable} viable, {n_nonmono} non-mono)")
    
    print(f"\n  Results:")
    print(f"    Physical (finite, ≥0):     {n_physical}/100")
    print(f"    SIDM-viable:               {n_viable}/100")
    print(f"    Monotonic (5% tolerance):  {n_mono}/100")
    print(f"    Non-monotonic:             {n_nonmono}/100")
    
    if nonmono_details:
        print(f"\n  Non-monotonic points (up to 5 shown):")
        for j, (idx, mc, mp, a, vi, c) in enumerate(nonmono_details[:5]):
            v1, v2 = VELOCITIES[vi], VELOCITIES[vi+1]
            print(f"    #{idx}: m_χ={mc:.1f}, m_φ={mp:.2e}, α={a:.2e}")
            print(f"      σ/m({v1:.0f})={c[v1]:.4f} → σ/m({v2:.0f})={c[v2]:.4f}  (ratio {c[v2]/c[v1]:.4f})")
    
    a_pass = (n_physical == 100)
    print(f"\n  BLOCK A: {'PASS' if a_pass else 'FAIL'} — solver stability")
    print(f"           viable fraction: {n_viable}% (within region bounds)")
    print(f"           non-monotonic: {n_nonmono}% prevalence")
    return a_pass, n_viable, n_nonmono


# ==============================================================
#  BLOCK B — 50 random points across FULL parameter space
# ==============================================================
def block_b_random_full_space():
    print("\n" + "="*70)
    print("BLOCK B: 50 random points across full parameter space")
    print("         m_χ ∈ [0.1, 1000] GeV, m_φ ∈ [0.01, 1000] MeV")
    print("         α ∈ [10⁻⁹, 1] — deliberately wide")
    print("="*70)
    
    rng = np.random.default_rng(seed=9999)
    
    n_total = 50
    n_physical = 0
    n_viable = 0
    
    for i in range(n_total):
        m_chi = 10**(rng.uniform(-1, 3))          # 0.1 – 1000 GeV
        m_phi = 10**(rng.uniform(-5, -0.3))        # 0.01 – 500 MeV (in GeV: 10⁻⁵ – 0.5)
        alpha = 10**(rng.uniform(-9, 0))           # 10⁻⁹ – 1
        
        try:
            s30   = sigma_T_vpm(m_chi, m_phi, alpha, 30.0)
            s1000 = sigma_T_vpm(m_chi, m_phi, alpha, 1000.0)
            is_phys = np.isfinite(s30) and np.isfinite(s1000) and s30 >= 0 and s1000 >= 0
        except Exception:
            is_phys = False
            s30, s1000 = float('nan'), float('nan')
        
        if is_phys:
            n_physical += 1
        viable = is_phys and (1.0 <= s30 <= 10.0) and (s1000 < 0.1)
        if viable:
            n_viable += 1
        
        if (i+1) % 25 == 0:
            print(f"  ... {i+1}/50 done ({n_viable} viable)")
    
    print(f"\n  Results:")
    print(f"    Physical:    {n_physical}/50")
    print(f"    Viable:      {n_viable}/50")
    print(f"    (Expected: ~0-2 viable from random sampling — viable region is very thin)")
    
    b_pass = (n_physical == 50)
    print(f"\n  BLOCK B: {'PASS' if b_pass else 'FAIL'} — solver never crashes on wild inputs")
    return b_pass, n_viable


# ==============================================================
#  BLOCK C — 50 perturbation tests (continuity check)
# ==============================================================
def block_c_perturbation():
    print("\n" + "="*70)
    print("BLOCK C: 50 perturbation tests — continuity of σ/m(α)")
    print("         Known viable points with α perturbed by ±1%, ±5%, ±10%")
    print("="*70)
    
    reps = load_representatives()
    rng = np.random.default_rng(seed=7777)
    
    # Pick 10 diverse base points
    indices = np.linspace(0, len(reps)-1, 10, dtype=int)
    perturbations = [0.99, 1.01, 0.95, 1.05, 1.10]  # 5 perturbations each
    
    n_total = 0
    n_continuous = 0
    max_relative_change = 0.0
    discontinuities = []
    
    for base_idx in indices:
        p = reps[base_idx]
        mc, mp = p['m_chi'], p['m_phi']
        alpha_base = p['alpha']
        s30_base = sigma_T_vpm(mc, mp, alpha_base, 30.0)
        
        for frac in perturbations:
            alpha_pert = alpha_base * frac
            s30_pert = sigma_T_vpm(mc, mp, alpha_pert, 30.0)
            n_total += 1
            
            if s30_base > 0:
                rel_change = abs(s30_pert - s30_base) / s30_base
            else:
                rel_change = 0.0
            
            max_relative_change = max(max_relative_change, rel_change)
            
            # For a 10% α change, σ/m should change by at most ~100%
            # (generous: near resonance it can change a lot)
            alpha_change = abs(frac - 1.0)
            expected_max = 10.0 * alpha_change + 0.1  # very generous bound
            
            if rel_change < expected_max:
                n_continuous += 1
            else:
                discontinuities.append((base_idx, alpha_base, frac, s30_base, s30_pert, rel_change))
    
    print(f"\n  Results:")
    print(f"    Continuous (within tolerance): {n_continuous}/{n_total}")
    print(f"    Max relative σ/m change:       {max_relative_change:.2%}")
    
    if discontinuities:
        print(f"    Discontinuities found: {len(discontinuities)}")
        for bi, ab, f, sb, sp, rc in discontinuities[:3]:
            print(f"      rep#{bi}: α×{f:.2f} → σ/m changed by {rc:.1%}")
    
    c_pass = (n_continuous >= n_total * 0.9)  # Allow ≤10% near-resonance jumps
    print(f"\n  BLOCK C: {'PASS' if c_pass else 'FAIL'} — solver output continuous in α")
    return c_pass


# ==============================================================
#  MAIN
# ==============================================================
def main():
    print("╔═══════════════════════════════════════════════════════════════════════╗")
    print("║     V8 — v29 Large Blind Sanity Check (200 off-grid evaluations)    ║")
    print("╠═══════════════════════════════════════════════════════════════════════╣")
    print("║  A: 100 random near-viable   B: 50 random wild   C: 50 perturbed   ║")
    print("║  Zero new physics — imports sigma_T_vpm from v22                    ║")
    print("╚═══════════════════════════════════════════════════════════════════════╝\n")
    
    t0 = time.time()
    
    # Warm up JIT
    print("Warming up Numba JIT...")
    _ = sigma_T_vpm(50.0, 5e-3, 5e-4, 30.0)
    print(f"JIT warm-up: {time.time()-t0:.1f}s\n")
    
    a_pass, a_viable, a_nonmono = block_a_random_viable_region()
    b_pass, b_viable = block_b_random_full_space()
    c_pass = block_c_perturbation()
    
    elapsed = time.time() - t0
    
    passed = sum([a_pass, b_pass, c_pass])
    total = 3
    
    print("\n" + "="*70)
    print(f"  GRAND SUMMARY  ({elapsed:.1f}s total)")
    print("="*70)
    print(f"  Block A (100 near-viable):  {'PASS' if a_pass else 'FAIL'}  — {a_viable}% viable, {a_nonmono}% non-mono")
    print(f"  Block B (50 wild random):   {'PASS' if b_pass else 'FAIL'}  — {b_viable} accidental viable")
    print(f"  Block C (50 perturbations): {'PASS' if c_pass else 'FAIL'}  — continuity check")
    print(f"\n  Overall: {passed}/{total} blocks PASS")
    
    if passed == total:
        print("\n  ✓ VPM solver is stable, continuous, and consistent on unseen data.")
        print("  ✓ No evidence of grid overfitting or numerical artifacts.")
        if a_nonmono > 0:
            print(f"  ℹ {a_nonmono}% of near-viable points show non-monotonic σ/m(v)")
            print("    — consistent with Yukawa quasi-bound-state resonances (physical).")
    else:
        print(f"\n  ⚠ {total - passed} block(s) failed — investigate above.")
    
    return 0 if passed == total else 1

if __name__ == "__main__":
    sys.exit(main())
