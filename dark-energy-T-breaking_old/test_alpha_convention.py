#!/usr/bin/env python3
"""
Compare two scenarios for 17 BPs:
  A) α_s = (8/9)α  (our original interpretation: α_CSV = α_total)
  B) α_s = α        (Copilot's claim: α_CSV = α_Yukawa already)

This settles the question empirically.
"""
import math, os, sys, csv
import numpy as np
from numba import jit

C_KM_S      = 2.998e5
GEV2_TO_CM2 = 3.8938e-28
GEV_IN_G    = 1.78266e-24
THETA        = math.asin(1.0 / 3.0)
COS2         = math.cos(THETA)**2   # 8/9

# VPM solver (same as consistency_17bp.py)
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

# Load BPs
CSV = os.path.join(os.path.dirname(__file__), "..",
    "Secluded-Majorana-SIDM", "predictions", "output", "sweep_17bp_results.csv")

bps = []
with open(CSV, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        bps.append({
            "label":     row["BP"],
            "m_chi":     float(row["m_chi_GeV"]),
            "m_phi":     float(row["m_phi_MeV"]) * 1e-3,
            "alpha":     float(row["alpha"]),
            "sm30_csv":  float(row["sigma_m_30"]),
            "sm1k_csv":  float(row["sigma_m_1000"]),
        })

def main():
    print("=" * 120)
    print("  CONVENTION TEST: α_s = (8/9)α  vs  α_s = α")
    print("=" * 120)

    print("\n  [JIT warmup] ", end="", flush=True)
    _ = sigma_T_vpm(20.0, 0.01, 1e-3, 30.0)
    print("done.\n")

    # Header
    print(f"  {'BP':>4} | {'m_χ':>7} | {'α_CSV':>10} | "
          f"{'σ/m(30)_CSV':>11} | {'σ/m(30)_A':>10} | {'σ/m(30)_B':>10} | "
          f"{'σ/m(1k)_A':>10} | {'σ/m(1k)_B':>10} | "
          f"{'SIDM_A':>6} | {'SIDM_B':>6} | {'CSV_match':>9}")
    print("  " + "-" * 116)

    n_pass_A = 0  # (8/9)α
    n_pass_B = 0  # α
    n_csv_match = 0

    for bp in bps:
        m_chi, m_phi, alpha = bp["m_chi"], bp["m_phi"], bp["alpha"]

        # Scenario A: α_s = (8/9)α
        alpha_A = alpha * COS2
        sm30_A  = sigma_T_vpm(m_chi, m_phi, alpha_A, 30.0)
        sm1k_A  = sigma_T_vpm(m_chi, m_phi, alpha_A, 1000.0)
        pass_A  = 0.5 <= sm30_A <= 10.0 and sm1k_A < 0.47

        # Scenario B: α_s = α (no reduction)
        alpha_B = alpha
        sm30_B  = sigma_T_vpm(m_chi, m_phi, alpha_B, 30.0)
        sm1k_B  = sigma_T_vpm(m_chi, m_phi, alpha_B, 1000.0)
        pass_B  = 0.5 <= sm30_B <= 10.0 and sm1k_B < 0.47

        if pass_A: n_pass_A += 1
        if pass_B: n_pass_B += 1

        # Does scenario B reproduce the CSV values?
        csv_match_30 = abs(sm30_B / bp["sm30_csv"] - 1) < 0.02  # within 2%
        csv_match_1k = abs(sm1k_B / bp["sm1k_csv"] - 1) < 0.02
        csv_match = csv_match_30 and csv_match_1k
        if csv_match: n_csv_match += 1

        print(f"  {bp['label']:>4} | {m_chi:7.1f} | {alpha:.3e} | "
              f"{bp['sm30_csv']:11.4f} | {sm30_A:10.4f} | {sm30_B:10.4f} | "
              f"{sm1k_A:10.5f} | {sm1k_B:10.5f} | "
              f"{'✓' if pass_A else '✗':>6} | {'✓' if pass_B else '✗':>6} | "
              f"{'✓' if csv_match else '✗':>9}")

    print("\n" + "=" * 120)
    print(f"  Scenario A  (α_s = 8/9 α):  {n_pass_A}/17 pass SIDM")
    print(f"  Scenario B  (α_s = α):       {n_pass_B}/17 pass SIDM")
    print(f"  CSV match   (B ≈ CSV ±2%):   {n_csv_match}/17")
    print()
    if n_csv_match >= 15:
        print("  → Scenario B matches CSV. α_CSV = α_Yukawa confirmed.")
        print("  → The (8/9) reduction was a convention error.")
    elif n_csv_match <= 2:
        print("  → Scenario B does NOT match CSV. α_CSV ≠ α_Yukawa.")
        print("  → Need to investigate what α_CSV actually represents.")
    else:
        print("  → Mixed results. Need deeper investigation.")

if __name__ == "__main__":
    main()
