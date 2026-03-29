#!/usr/bin/env python3
# ⚠️ SUPERSEDED — This script applies (8/9) reduction to α_CSV,
# which is already α_s. See test_alpha_convention.py (Test 12).
# The "without decomposition" column is the correct result.
# Buggy output archived in archived_buggy/fornax_gc_check_BUGGY.txt
"""
Fornax GC timing constraint (Read+2019, 1808.06634):
  σ/m < 1.5 cm²/g  at  v ~ 10-20 km/s

Test all 17 relic-viable BPs with α_s = (8/9)α.
Compare to results WITHOUT θ-decomposition (full α).
"""

import math, os, sys, csv
import numpy as np
from numba import jit

# =============================================================
#  Constants
# =============================================================
C_KM_S      = 2.998e5          # c in km/s
GEV2_TO_CM2 = 3.8938e-28
GEV_IN_G    = 1.78266e-24
THETA        = math.asin(1.0 / 3.0)
COS2         = math.cos(THETA)**2   # 8/9
FORNAX_BOUND = 1.5                  # cm²/g

VELOCITIES = [10.0, 12.0, 15.0, 20.0]  # km/s

# =============================================================
#  VPM solver (copied from consistency_17bp.py, unmodified)
# =============================================================
@jit(nopython=True, cache=True)
def sph_jn_numba(l, z):
    if z < 1e-30: return 1.0 if l == 0 else 0.0
    j0 = math.sin(z) / z
    if l == 0: return j0
    j1 = math.sin(z) / (z * z) - math.cos(z) / z
    if l == 1: return j1
    j_prev, j_curr = j0, j1
    for n in range(1, l):
        j_next = (2 * n + 1) / z * j_curr - j_prev
        j_prev, j_curr = j_curr, j_next
        if abs(j_curr) < 1e-300: return 0.0
    return j_curr

@jit(nopython=True, cache=True)
def sph_yn_numba(l, z):
    if z < 1e-30: return -1e300
    y0 = -math.cos(z) / z
    if l == 0: return y0
    y1 = -math.cos(z) / (z * z) - math.sin(z) / z
    if l == 1: return y1
    y_prev, y_curr = y0, y1
    for n in range(1, l):
        y_next = (2 * n + 1) / z * y_curr - y_prev
        y_prev, y_curr = y_curr, y_next
        if abs(y_curr) > 1e200: return y_curr
    return y_curr

@jit(nopython=True, cache=True)
def _vpm_rhs(l, kappa, lam, x, delta):
    if x < 1e-20: return 0.0
    z = kappa * x
    if z < 1e-20: return 0.0
    jl = sph_jn_numba(l, z)
    nl = sph_yn_numba(l, z)
    j_hat = z * jl
    n_hat = -z * nl
    cd, sd = math.cos(delta), math.sin(delta)
    bracket = j_hat * cd - n_hat * sd
    if not math.isfinite(bracket): return 0.0
    pot = lam * math.exp(-x) / (kappa * x)
    val = pot * bracket * bracket
    return val if math.isfinite(val) else 0.0

@jit(nopython=True, cache=True)
def vpm_phase_shift(l, kappa, lam, x_max=50.0, N_steps=4000):
    if lam < 1e-30 or kappa < 1e-30: return 0.0
    x_min = max(1e-5, 0.05 / (kappa + 0.01))
    if l > 0:
        x_barrier = l / kappa
        if x_barrier > x_min: x_min = x_barrier
    h = (x_max - x_min) / N_steps
    delta = 0.0
    for i in range(N_steps):
        x = x_min + i * h
        k1 = _vpm_rhs(l, kappa, lam, x, delta)
        k2 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k1)
        k3 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k2)
        k4 = _vpm_rhs(l, kappa, lam, x + h, delta + h*k3)
        delta += h * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return delta

@jit(nopython=True, cache=True)
def sigma_T_vpm(m_chi, m_phi, alpha, v_km_s):
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15: return 0.0
    if kappa < 5:     x_max, N_steps = 50.0, 4000
    elif kappa < 50:  x_max, N_steps = 80.0, 8000
    else:             x_max, N_steps = 100.0, 12000
    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)
    sigma_sum = 0.0
    peak = 0.0
    n_small = 0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        w = 1.0 if l % 2 == 0 else 3.0
        c = w * (2*l + 1) * math.sin(delta)**2
        sigma_sum += c
        if c > peak: peak = c
        if peak > 0 and c / peak < 1e-4:
            n_small += 1
            if n_small >= 5: break
        else: n_small = 0
    sigma_GeV2 = 2.0 * math.pi * sigma_sum / (k * k)
    return sigma_GeV2 * GEV2_TO_CM2 / (m_chi * GEV_IN_G)

# =============================================================
#  Load 17 BPs
# =============================================================
CSV = os.path.join(os.path.dirname(__file__), "..",
    "Secluded-Majorana-SIDM", "predictions", "output", "sweep_17bp_results.csv")

bps = []
with open(CSV, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        bps.append({
            "label":    row["BP"],
            "m_chi":    float(row["m_chi_GeV"]),
            "m_phi":    float(row["m_phi_MeV"]) * 1e-3,
            "alpha":    float(row["alpha"]),
            "lam_orig": float(row["lambda"]),
        })

# =============================================================
#  Main
# =============================================================
def main():
    print("=" * 100)
    print("  FORNAX GC TIMING CONSTRAINT — Read+2019 (1808.06634)")
    print(f"  Bound: σ/m < {FORNAX_BOUND} cm²/g at v = 10-20 km/s")
    print(f"  θ = arcsin(1/3) = {math.degrees(THETA):.2f}°  →  α_s = (8/9)α")
    print("=" * 100)

    # JIT warmup
    print("\n  [JIT warmup] ", end="", flush=True)
    _ = sigma_T_vpm(20.0, 0.01, 1e-3, 30.0)
    print("done.\n")

    # Header
    hdr = f"{'BP':>6} | {'m_χ':>7} | {'m_φ':>7} | {'α':>9} | {'α_s':>9} | {'λ_s':>6}"
    for v in VELOCITIES:
        hdr += f" | σ/m({int(v):>2})"
    hdr += " |  worst  | margin | Pass?"
    print(hdr)
    print("-" * len(hdr))

    n_pass_decomp = 0
    n_pass_full   = 0
    results = []

    for bp in bps:
        m_chi, m_phi, alpha = bp["m_chi"], bp["m_phi"], bp["alpha"]
        alpha_s = alpha * COS2
        lam_s   = alpha_s * m_chi / m_phi

        # σ/m at each velocity — WITH θ-decomposition
        sm_decomp = {}
        for v in VELOCITIES:
            sm_decomp[v] = sigma_T_vpm(m_chi, m_phi, alpha_s, v)

        # σ/m at each velocity — WITHOUT θ-decomposition (full α)
        sm_full = {}
        for v in VELOCITIES:
            sm_full[v] = sigma_T_vpm(m_chi, m_phi, alpha, v)

        worst_decomp = max(sm_decomp.values())
        worst_full   = max(sm_full.values())
        margin_decomp = FORNAX_BOUND - worst_decomp
        pass_decomp = worst_decomp < FORNAX_BOUND
        pass_full   = worst_full < FORNAX_BOUND

        if pass_decomp: n_pass_decomp += 1
        if pass_full:   n_pass_full += 1

        results.append({
            "label": bp["label"], "m_chi": m_chi, "m_phi": m_phi,
            "alpha": alpha, "alpha_s": alpha_s, "lam_s": lam_s,
            "sm_decomp": sm_decomp, "sm_full": sm_full,
            "worst_decomp": worst_decomp, "worst_full": worst_full,
            "margin": margin_decomp, "pass_decomp": pass_decomp,
            "pass_full": pass_full,
        })

        line = f"{bp['label']:>6} | {m_chi:7.1f} | {m_phi*1e3:6.2f}M | {alpha:.2e} | {alpha_s:.2e} | {lam_s:6.2f}"
        for v in VELOCITIES:
            val = sm_decomp[v]
            flag = "!" if val > FORNAX_BOUND else " "
            line += f" | {val:6.3f}{flag}"
        line += f" | {worst_decomp:7.3f} | {margin_decomp:+6.3f} | {'✓' if pass_decomp else '✗'}"
        print(line)

    # Summary
    print("\n" + "=" * 100)
    print("  SUMMARY")
    print("=" * 100)
    print(f"  With θ-decomposition (α_s = 8/9 α): {n_pass_decomp}/17 pass Fornax GC")
    print(f"  Without θ-decomposition (full α):     {n_pass_full}/17 pass Fornax GC")
    print(f"  Improvement from θ-decomposition:      {n_pass_decomp - n_pass_full:+d} BPs saved")

    # Comparison table
    print(f"\n{'BP':>6} | {'worst(α_s)':>10} | {'worst(α)':>10} | {'reduction':>9} | {'α_s pass':>8} | {'α pass':>8}")
    print("-" * 72)
    for r in results:
        red = 1.0 - r["worst_decomp"] / r["worst_full"] if r["worst_full"] > 0 else 0
        print(f"{r['label']:>6} | {r['worst_decomp']:10.3f} | {r['worst_full']:10.3f} | "
              f"{red*100:8.1f}% | {'✓' if r['pass_decomp'] else '✗':>8} | {'✓' if r['pass_full'] else '✗':>8}")

    # Cross-reference with σ/m(30) viability
    print(f"\n{'BP':>6} | {'σ/m(30)_αs':>11} | {'σ/m(1k)_αs':>11} | {'Fornax':>6} | {'SIDM+Fornax':>11}")
    print("-" * 65)
    for r in results:
        sm30  = sigma_T_vpm(r["m_chi"], r["m_phi"], r["alpha_s"], 30.0)
        sm1k  = sigma_T_vpm(r["m_chi"], r["m_phi"], r["alpha_s"], 1000.0)
        sidm_ok = 0.5 <= sm30 <= 10.0 and sm1k < 0.47
        both = sidm_ok and r["pass_decomp"]
        print(f"{r['label']:>6} | {sm30:11.3f} | {sm1k:11.4f} | {'✓' if r['pass_decomp'] else '✗':>6} | {'✓' if both else '✗':>11}")

    n_both = sum(1 for r in results
                 if r["pass_decomp"]
                 and 0.5 <= sigma_T_vpm(r["m_chi"], r["m_phi"], r["alpha_s"], 30.0) <= 10.0
                 and sigma_T_vpm(r["m_chi"], r["m_phi"], r["alpha_s"], 1000.0) < 0.47)
    print(f"\n  BPs passing ALL constraints (SIDM + cluster + Fornax): {n_both}/17")

if __name__ == "__main__":
    main()
