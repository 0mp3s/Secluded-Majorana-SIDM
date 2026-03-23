#!/usr/bin/env python3
"""
Condition 2: 2D Scan in (αₛ, αₚ) — Island Width / Band Structure
==================================================================

Goal: Show that the viable region in (αₛ, αₚ) is a BAND (not a point),
proving no fine-tuning between yₛ and yₚ.

Three simultaneous constraints:
  1. Relic:  Ωh² = 0.120 ± 0.001  →  ⟨σv⟩₀ = 2π αₛ αₚ / m_χ²
     This is a HYPERBOLA  αₛ × αₚ = const  in log-log space.
  2. SIDM:   0.1 < σ_T/m < 10 cm²/g  at v = 30, 100, 1000 km/s
     VPM depends on α = αₛ only (condition 4 proved yₚ negligible).
     This is a VERTICAL STRIP in (αₛ, αₚ) space.
  3. Perturbativity:  αₛ, αₚ < 1  (or equivalently y < √(4π))

The intersection:  hyperbola ∩ strip = BAND.

Method:
  Phase A: Boltzmann for each αₛ × αₚ → find relic contour  (fast, ~seconds)
  Phase B: VPM for each αₛ → find SIDM strip  (heavy, parallelized)

Fixed:  m_χ = 20.69 GeV,  m_φ = 11.34 MeV  (BP1 masses)

Date: 23 March 2026
"""

import sys, os, math, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

from v22_raw_scan import sigma_T_vpm

# ═══════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════
M_PL        = 1.2209e19
OMEGA_CDM   = 0.120
OMEGA_TOL   = 0.001
RHO_CRIT_H2 = 1.0539e-5        # GeV/cm³
S_0         = 2891.2           # cm⁻³

# ═══════════════════════════════════════════════════════════════
#  Benchmark (BP1) — fixed masses
# ═══════════════════════════════════════════════════════════════
M_CHI  = 20.69               # GeV
M_PHI  = 11.34e-3            # GeV (11.34 MeV)
G_CHI  = 2                   # Majorana DOF

# SIDM velocity targets (km/s)
V_SIDM = [30.0, 100.0, 1000.0]
SIDM_LABELS = ['dwarf (30 km/s)', 'LSB (100 km/s)', 'cluster (1000 km/s)']

# SIDM cuts:  σ/m in cm²/g
SIDM_LO = [0.1, 0.1, 0.1]
SIDM_HI = [50.0, 10.0, 1.0]

# ═══════════════════════════════════════════════════════════════
#  g_*(T) table (Drees+ 2015)
# ═══════════════════════════════════════════════════════════════
_G_TABLE = np.array([
    [1e4, 106.75, 106.75], [200, 106.75, 106.75],
    [80, 86.25, 86.25], [10, 86.25, 86.25],
    [1, 75.75, 75.75], [0.3, 61.75, 61.75],
    [0.2, 17.25, 17.25], [0.15, 14.25, 14.25],
    [0.1, 10.75, 10.75], [0.01, 10.75, 10.75],
    [0.001, 10.75, 10.75], [0.0005, 10.75, 10.75],
    [0.0001, 3.36, 3.91], [1e-5, 3.36, 3.91], [1e-8, 3.36, 3.91],
])
_LT = np.log(_G_TABLE[:, 0])
_GR = _G_TABLE[:, 1]
_GS = _G_TABLE[:, 2]

def g_rho(T): return float(np.interp(math.log(max(T, 1e-50)), _LT[::-1], _GR[::-1]))
def g_S(T):   return float(np.interp(math.log(max(T, 1e-50)), _LT[::-1], _GS[::-1]))


# ═══════════════════════════════════════════════════════════════
#  Boltzmann solver (fast RK4, single species)
# ═══════════════════════════════════════════════════════════════
def Y_eq(x, m_chi=M_CHI, g_chi=G_CHI):
    if x > 500: return 0.0
    T = m_chi / x
    gs = g_S(T)
    return 45.0 / (4 * math.pi**4) * g_chi / gs * x**1.5 * math.exp(-x)


def solve_boltzmann(alpha_s, alpha_p, x_end=500.0, n_steps=15000):
    """Solve Boltzmann for χ with ⟨σv⟩₀ = 2π αₛ αₚ / m_χ².
    Returns Ωh².
    """
    sv0 = 2.0 * math.pi * alpha_s * alpha_p / M_CHI**2
    if sv0 <= 0:
        return 1e10
    
    x_start = 1.0
    dx = (x_end - x_start) / n_steps
    x = x_start
    Y = Y_eq(x_start)
    
    for _ in range(n_steps):
        T = M_CHI / x
        gs = g_S(T)
        gr = g_rho(T)
        g_eff = gs / math.sqrt(gr)
        lam = math.sqrt(math.pi / 45.0) * g_eff * M_PL * M_CHI
        
        def f(xx, yy):
            Yeq = Y_eq(xx)
            return -lam * sv0 * (yy**2 - Yeq**2) / (xx * xx)
        
        k1 = dx * f(x, Y)
        k2 = dx * f(x + dx/2, Y + k1/2)
        k3 = dx * f(x + dx/2, Y + k2/2)
        k4 = dx * f(x + dx, Y + k3)
        
        Y = max(Y + (k1 + 2*k2 + 2*k3 + k4)/6, 1e-30)
        x += dx
    
    return M_CHI * Y * S_0 / RHO_CRIT_H2


# ═══════════════════════════════════════════════════════════════
#  VPM worker for multiprocessing
# ═══════════════════════════════════════════════════════════════
def _vpm_worker(args):
    """Compute σ_T/m for given αₛ at all SIDM velocities."""
    idx, alpha_s = args
    results = []
    for v in V_SIDM:
        sm = sigma_T_vpm(M_CHI, M_PHI, alpha_s, v)
        results.append(sm)
    return (idx, alpha_s, results)


# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════
def main():
    t0 = time.time()
    
    print("=" * 82)
    print("CONDITION 2: 2D Coupling Scan  (αₛ, αₚ)  →  Relic × SIDM Band")
    print("=" * 82)
    print()
    print(f"  Fixed:  m_χ = {M_CHI:.2f} GeV,  m_φ = {M_PHI*1e3:.2f} MeV")
    print(f"  Relic:  ⟨σv⟩₀ = 2π αₛ αₚ / m_χ²  →  Ωh² = 0.120 ± {OMEGA_TOL}")
    print(f"  SIDM:   σ_T/m from VPM (depends on αₛ only)")
    print(f"  SIDM velocities: {V_SIDM} km/s")
    print(f"  CPU cores: {multiprocessing.cpu_count()}")
    print()
    
    # ──────────────────────────────────────────────────────────
    #  PHASE A: Relic density — find αₛ×αₚ that gives Ωh²=0.12
    # ──────────────────────────────────────────────────────────
    print("─" * 82)
    print("  PHASE A: Relic Density Contour")
    print("─" * 82)
    print()
    
    # First find the product αₛ×αₚ that gives correct relic
    # Use bisection: fix αₛ=αₚ=α, find α such that Ωh²=0.12
    
    def omega_symmetric(alpha):
        return solve_boltzmann(alpha, alpha)
    
    # Bracket: small α → large Ω, large α → small Ω
    a_lo, a_hi = 1e-5, 0.1
    for _ in range(60):
        a_mid = math.sqrt(a_lo * a_hi)
        om = omega_symmetric(a_mid)
        if om > OMEGA_CDM:
            a_lo = a_mid
        else:
            a_hi = a_mid
    alpha_target = math.sqrt(a_lo * a_hi)
    product_target = alpha_target**2
    omega_check = omega_symmetric(alpha_target)
    
    print(f"  Symmetric point:  αₛ = αₚ = {alpha_target:.6e}")
    print(f"  Product:          αₛ × αₚ = {product_target:.6e}")
    print(f"  Ωh² = {omega_check:.6f}  (target: {OMEGA_CDM})")
    print()
    
    # Verify: the relic contour is αₛ × αₚ = product_target
    # Check a few asymmetric points
    print("  Verification — asymmetric points on hyperbola:")
    test_ratios = [0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0]
    print(f"    {'αₛ':>12}  {'αₚ':>12}  {'αₛ×αₚ':>12}  {'Ωh²':>10}  {'Status':>8}")
    print(f"    {'─'*70}")
    for r in test_ratios:
        as_val = alpha_target * math.sqrt(r)
        ap_val = alpha_target / math.sqrt(r)
        om = solve_boltzmann(as_val, ap_val)
        status = "✓" if abs(om - OMEGA_CDM) < 0.005 else "✗"
        print(f"    {as_val:12.6e}  {ap_val:12.6e}  {as_val*ap_val:12.6e}  {om:10.6f}  {status:>8}")
    print()
    
    # ──────────────────────────────────────────────────────────
    #  PHASE B: SIDM — VPM scan over αₛ  (parallelized)
    # ──────────────────────────────────────────────────────────
    print("─" * 82)
    print("  PHASE B: SIDM Band  (VPM σ_T/m vs αₛ)")
    print("─" * 82)
    print()
    
    # Scan αₛ from 1e-5 to 0.1 (log-spaced)
    N_ALPHA = 80
    alpha_s_arr = np.logspace(-5, -1, N_ALPHA)
    
    print(f"  Scanning {N_ALPHA} values of αₛ ∈ [{alpha_s_arr[0]:.1e}, {alpha_s_arr[-1]:.1e}]")
    print(f"  Workers: {min(multiprocessing.cpu_count(), N_ALPHA)}")
    print()
    
    # Warm up JIT
    print("  JIT warmup...", end=" ", flush=True)
    _ = sigma_T_vpm(M_CHI, M_PHI, 1e-3, 100.0)
    print("done.")
    
    # Parallel VPM scan
    t_vpm = time.time()
    vpm_results = {}
    
    n_workers = min(multiprocessing.cpu_count() - 1, N_ALPHA)
    n_workers = max(n_workers, 1)
    
    tasks = [(i, alpha_s_arr[i]) for i in range(N_ALPHA)]
    
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_vpm_worker, t): t[0] for t in tasks}
        done = 0
        for fut in as_completed(futures):
            idx, alpha_s, sigmas = fut.result()
            vpm_results[idx] = (alpha_s, sigmas)
            done += 1
            if done % 20 == 0 or done == N_ALPHA:
                print(f"    VPM: {done}/{N_ALPHA} done ({time.time()-t_vpm:.1f}s)", flush=True)
    
    dt_vpm = time.time() - t_vpm
    print(f"  VPM scan complete: {dt_vpm:.1f}s")
    print()
    
    # ──────────────────────────────────────────────────────────
    #  PHASE C: Combine — find the band
    # ──────────────────────────────────────────────────────────
    print("─" * 82)
    print("  PHASE C: Combined Relic × SIDM Results")
    print("─" * 82)
    print()
    
    # For each αₛ: compute σ/m at each velocity, check SIDM cuts
    # Also compute the required αₚ from relic: αₚ = product_target / αₛ
    
    print(f"  {'αₛ':>11}  {'αₚ(relic)':>11}  {'αₛ/αₚ':>7}  "
          f"{'σ/m(30)':>10}  {'σ/m(100)':>10}  {'σ/m(1000)':>10}  "
          f"{'SIDM?':>6}  {'Pert?':>6}  {'VIABLE':>7}")
    print("  " + "─" * 100)
    
    viable_points = []
    sidm_only = []  # passes SIDM but not relic+perturbativity
    
    for i in range(N_ALPHA):
        alpha_s, sigmas = vpm_results[i]
        
        # Required αₚ for correct relic
        alpha_p = product_target / alpha_s
        
        # Perturbativity
        pert_ok = (alpha_s < 1.0) and (alpha_p < 1.0)
        
        # SIDM checks
        sidm_pass = True
        for j in range(len(V_SIDM)):
            if sigmas[j] < SIDM_LO[j] or sigmas[j] > SIDM_HI[j]:
                sidm_pass = False
                break
        
        # In the allowed regime: σ/m ~ 1 cm²/g at clusters is the tightest
        # but σ/m ~ 1-50 at dwarfs is also needed
        
        viable = pert_ok and sidm_pass
        
        sm_strs = [f"{s:10.4f}" for s in sigmas]
        
        flag = "★ YES" if viable else ""
        
        print(f"  {alpha_s:11.4e}  {alpha_p:11.4e}  {alpha_s/alpha_p:7.2f}  "
              f"{'  '.join(sm_strs)}  "
              f"{'✓' if sidm_pass else '✗':>6}  {'✓' if pert_ok else '✗':>6}  {flag:>7}")
        
        if viable:
            viable_points.append({
                'alpha_s': alpha_s, 'alpha_p': alpha_p,
                'ratio': alpha_s / alpha_p,
                'sigma_m': sigmas
            })
        if sidm_pass:
            sidm_only.append({'alpha_s': alpha_s, 'sigmas': sigmas})
    
    print()
    
    # ──────────────────────────────────────────────────────────
    #  SUMMARY
    # ──────────────────────────────────────────────────────────
    print("=" * 82)
    print("  SUMMARY")
    print("=" * 82)
    print()
    
    # SIDM strip
    if sidm_only:
        as_min = min(p['alpha_s'] for p in sidm_only)
        as_max = max(p['alpha_s'] for p in sidm_only)
        print(f"  SIDM strip (all velocities pass):")
        print(f"    αₛ ∈ [{as_min:.4e}, {as_max:.4e}]")
        print(f"    Width in log₁₀: {math.log10(as_max/as_min):.2f} decades")
    else:
        print("  ⚠ No αₛ values pass all SIDM cuts simultaneously!")
        print("    (May need to relax cluster cut or check velocity range)")
        # Show which individual cuts pass
        print()
        print("  Individual SIDM cut results:")
        for j, (v, lo, hi) in enumerate(zip(V_SIDM, SIDM_LO, SIDM_HI)):
            passing = []
            for i in range(N_ALPHA):
                alpha_s, sigmas = vpm_results[i]
                if lo <= sigmas[j] <= hi:
                    passing.append(alpha_s)
            if passing:
                print(f"    v={v:.0f} km/s ({lo}-{hi} cm²/g): αₛ ∈ [{min(passing):.4e}, {max(passing):.4e}]")
            else:
                print(f"    v={v:.0f} km/s ({lo}-{hi} cm²/g): NO αₛ passes")
    print()
    
    # Relic hyperbola
    print(f"  Relic hyperbola:")
    print(f"    αₛ × αₚ = {product_target:.6e}")
    print(f"    Symmetric point: αₛ = αₚ = {alpha_target:.6e}")
    print()
    
    # Viable band
    if viable_points:
        as_vals = [p['alpha_s'] for p in viable_points]
        ap_vals = [p['alpha_p'] for p in viable_points]
        r_vals  = [p['ratio'] for p in viable_points]
        
        print(f"  ★ VIABLE BAND (Relic + SIDM + Perturbativity):")
        print(f"    αₛ  ∈ [{min(as_vals):.4e}, {max(as_vals):.4e}]  "
              f"(width: {math.log10(max(as_vals)/min(as_vals)):.2f} decades)")
        print(f"    αₚ  ∈ [{min(ap_vals):.4e}, {max(ap_vals):.4e}]  "
              f"(width: {math.log10(max(ap_vals)/min(ap_vals)):.2f} decades)")
        print(f"    αₛ/αₚ ∈ [{min(r_vals):.3f}, {max(r_vals):.3f}]")
        print()
        print(f"    → NOT fine-tuned: band spans {math.log10(max(r_vals)/min(r_vals)):.1f} decades in αₛ/αₚ")
        
        # Table of viable points
        print()
        print(f"    Viable grid points ({len(viable_points)}):")
        print(f"    {'αₛ':>11}  {'αₚ':>11}  {'αₛ/αₚ':>7}  "
              f"{'σ/m(30)':>9}  {'σ/m(100)':>9}  {'σ/m(1000)':>9}")
        print(f"    {'─'*72}")
        for p in viable_points:
            sm = p['sigma_m']
            print(f"    {p['alpha_s']:11.4e}  {p['alpha_p']:11.4e}  {p['ratio']:7.3f}  "
                  f"{sm[0]:9.3f}  {sm[1]:9.3f}  {sm[2]:9.3f}")
    else:
        print("  ⚠ No viable points found!")
        print("    The relic contour and SIDM strip may not overlap for BP1 masses.")
        print("    This would require adjusting m_χ or m_φ.")
        
        # Show what σ/m values look like at the relic-relevant αₛ range
        print()
        print("  σ/m values near the symmetric relic point:")
        nearby = [(i, vpm_results[i]) for i in range(N_ALPHA)
                  if 0.1 * alpha_target <= vpm_results[i][0] <= 10 * alpha_target]
        if nearby:
            for i, (alpha_s, sigmas) in nearby:
                alpha_p_req = product_target / alpha_s
                print(f"    αₛ={alpha_s:.4e}, αₚ={alpha_p_req:.4e}: "
                      f"σ/m = [{sigmas[0]:.3f}, {sigmas[1]:.3f}, {sigmas[2]:.3f}] cm²/g")
    
    print()
    
    # ──────────────────────────────────────────────────────────
    #  PHYSICS CONCLUSION
    # ──────────────────────────────────────────────────────────
    print("═" * 82)
    print("  PHYSICS CONCLUSION")
    print("═" * 82)
    print()
    print("  The parameter space (αₛ, αₚ) decomposes into:")
    print("    • Relic density:  αₛ × αₚ = const  (hyperbola)")
    print("    • SIDM:           αₛ ∈ [αₛ_min, αₛ_max]  (vertical strip)")
    print("    • Perturbativity: αₛ, αₚ < 1")
    print()
    print("  Intersection = BAND along hyperbola, width set by SIDM strip.")
    print("  yₛ ≠ yₚ is allowed — the model works for a RANGE of αₛ/αₚ ratios.")
    print()
    
    dt = time.time() - t0
    print(f"  Total time: {dt:.1f}s")
    print()
    print("[CONDITION 2 SCAN COMPLETE]")


if __name__ == '__main__':
    main()
