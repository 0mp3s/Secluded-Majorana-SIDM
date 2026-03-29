#!/usr/bin/env python3
# ⚠️ SUPERSEDED — This script applies (8/9) reduction to α_CSV,
# which is already α_s. See test_alpha_convention.py (Test 12).
# Correct results: all 17 BPs pass with α_s = α_CSV (no reduction).
# Buggy output archived in archived_buggy/consistency_17bp_BUGGY.txt
"""
dark-energy-T-breaking/consistency_17bp.py
==========================================
Run all 17 relic benchmark points through θ-decomposition.

Input: sweep_17bp_results.csv from Secluded-Majorana-SIDM project.
Test:  α → α_s = (8/9)α, recompute σ_T/m via VPM, check SIDM cuts.

Goal: Find PATTERN in which BPs survive and which don't.
"""
import sys, math, os
import numpy as np

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

try:
    from numba import jit
except ImportError:
    def jit(**kwargs):
        def decorator(func): return func
        return decorator

# ==============================================================
#  Constants & angle (same as consistency_check_sidm.py)
# ==============================================================
THETA = math.asin(1.0 / 3.0)
COS2 = math.cos(THETA)**2      # 8/9
SIN2 = math.sin(THETA)**2      # 1/9

GEV2_TO_CM2 = 3.8938e-28
GEV_IN_G    = 1.78266e-24
C_KM_S      = 299792.458
M_PL        = 1.220890e19
S_0         = 2891.2
RHO_CRIT_H2 = 1.0539e-5

SIGMA_M_30_LO  = 0.5
SIGMA_M_30_HI  = 10.0
SIGMA_M_1000_HI = 0.47

# ==============================================================
#  VPM solver (copied from Secluded-Majorana-SIDM, unmodified)
# ==============================================================
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

# ==============================================================
#  Relic (Kolb-Turner, from v27_boltzmann_relic.py)
# ==============================================================
_G_STAR_TABLE = np.array([
    [1e4,106.75,106.75],[200,106.75,106.75],[80,86.25,86.25],
    [10,86.25,86.25],[1,75.75,75.75],[0.3,61.75,61.75],
    [0.2,17.25,17.25],[0.15,14.25,14.25],[0.1,10.75,10.75],
    [0.01,10.75,10.75],[0.001,10.75,10.75],[0.0005,10.75,10.75],
    [0.0001,3.36,3.91],[1e-5,3.36,3.91],[1e-8,3.36,3.91],
])
_LOG_T = np.log(_G_STAR_TABLE[:, 0])
_G_RHO = _G_STAR_TABLE[:, 1]
_G_S   = _G_STAR_TABLE[:, 2]

def g_star_rho(T):
    return float(np.interp(math.log(T) if T > 0 else -50, _LOG_T[::-1], _G_RHO[::-1]))
def g_star_S(T):
    return float(np.interp(math.log(T) if T > 0 else -50, _LOG_T[::-1], _G_S[::-1]))

def kolb_turner(m_chi, alpha_s):
    sv0 = math.pi * alpha_s**2 / (4.0 * m_chi**2)
    x_fo = 20.0
    for _ in range(50):
        gr = g_star_rho(m_chi / x_fo)
        arg = 0.0764 * 2 * M_PL * m_chi * sv0 / math.sqrt(gr * x_fo)
        if arg <= 0: break
        x_new = math.log(arg) - 0.5 * math.log(x_fo)
        if abs(x_new - x_fo) < 0.001: x_fo = x_new; break
        x_fo = 0.5 * x_fo + 0.5 * x_new
    gs = g_star_S(m_chi / x_fo)
    gr = g_star_rho(m_chi / x_fo)
    lam = math.sqrt(math.pi / 45.0) * gs / math.sqrt(gr) * M_PL * m_chi
    Y_inf = x_fo / (lam * sv0)
    omega = m_chi * Y_inf * S_0 / RHO_CRIT_H2
    return x_fo, omega

# ==============================================================
#  Load 17 BPs
# ==============================================================
CSV = os.path.join(os.path.dirname(__file__), "..",
    "Secluded-Majorana-SIDM", "predictions", "output", "sweep_17bp_results.csv")

import csv
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
            "sm30_orig": float(row["sigma_m_30"]),
            "sm1k_orig": float(row["sigma_m_1000"]),
        })

# Yukawa resonance positions
LAM_CRITS = [1.68, 6.45, 14.7, 26.0]

def nearest_resonance(lam):
    dists = [(abs(lam - lc), i+1, lc) for i, lc in enumerate(LAM_CRITS)]
    dists.sort()
    return dists[0][1], dists[0][2], dists[0][0]

# ==============================================================
#  Main
# ==============================================================
def main():
    print("=" * 110)
    print("  17 RELIC BPs — θ-DECOMPOSITION CONSISTENCY CHECK")
    print("=" * 110)
    print(f"  θ = arcsin(1/3) = {math.degrees(THETA):.2f}°  →  α_s = (8/9)α  →  λ_new = (8/9)λ")
    print()

    # JIT warmup
    print("  [JIT warmup] ", end="", flush=True)
    _ = sigma_T_vpm(20.0, 0.01, 1e-3, 30.0)
    print("done.")
    print()

    # Compute for all 17
    results = []
    for i, bp in enumerate(bps):
        alpha_s = bp["alpha"] * COS2
        m_chi, m_phi = bp["m_chi"], bp["m_phi"]
        lam_new = alpha_s * m_chi / m_phi

        sm30  = sigma_T_vpm(m_chi, m_phi, alpha_s, 30.0)
        sm1k  = sigma_T_vpm(m_chi, m_phi, alpha_s, 1000.0)

        _, omega = kolb_turner(m_chi, alpha_s)

        pass_30 = SIGMA_M_30_LO <= sm30 <= SIGMA_M_30_HI
        pass_1k = sm1k < SIGMA_M_1000_HI

        res_n, res_lam, res_dist = nearest_resonance(lam_new)
        crossed = nearest_resonance(bp["lam_orig"])[0] != res_n

        pct_30 = (sm30 / bp["sm30_orig"] - 1) * 100
        pct_1k = (sm1k / bp["sm1k_orig"] - 1) * 100

        results.append({
            "label": bp["label"], "m_chi": m_chi, "m_phi": m_phi*1e3,
            "alpha_orig": bp["alpha"], "alpha_s": alpha_s,
            "lam_orig": bp["lam_orig"], "lam_new": lam_new,
            "sm30_orig": bp["sm30_orig"], "sm30_new": sm30,
            "sm1k_orig": bp["sm1k_orig"], "sm1k_new": sm1k,
            "pct_30": pct_30, "pct_1k": pct_1k,
            "omega": omega,
            "pass_30": pass_30, "pass_1k": pass_1k,
            "pass_all": pass_30 and pass_1k,
            "res_n": res_n, "res_lam": res_lam, "res_dist": res_dist,
            "crossed_res": crossed,
        })
        print(f"  [{i+1:2d}/17] {bp['label']:>4}  λ {bp['lam_orig']:.2f}→{lam_new:.2f}  "
              f"σ/m(30) {bp['sm30_orig']:.3f}→{sm30:.3f}  "
              f"{'✓' if pass_30 and pass_1k else '✗'}")

    # ==============================================================
    #  Results table
    # ==============================================================
    print()
    print("=" * 110)
    print("  FULL RESULTS")
    print("=" * 110)
    hdr = (f"  {'BP':>4}  {'m_χ':>7}  {'λ_old':>6}  {'λ_new':>6}  {'res#':>4}  "
           f"{'σ/m30_old':>9}  {'σ/m30_new':>9}  {'Δ%':>6}  "
           f"{'σ/m1k_new':>9}  {'Ωh²_KT':>8}  {'SIDM':>4}  {'cross?':>6}")
    print(hdr)
    print("  " + "-" * 106)

    n_pass = 0
    n_crossed = 0
    for r in results:
        tag = "✓" if r["pass_all"] else "✗"
        cross_tag = "YES!" if r["crossed_res"] else ""
        if r["pass_all"]: n_pass += 1
        if r["crossed_res"]: n_crossed += 1
        print(f"  {r['label']:>4}  {r['m_chi']:7.2f}  {r['lam_orig']:6.2f}  {r['lam_new']:6.2f}  "
              f"  R{r['res_n']}  "
              f"{r['sm30_orig']:9.4f}  {r['sm30_new']:9.4f}  {r['pct_30']:+5.1f}%  "
              f"{r['sm1k_new']:9.5f}  {r['omega']:8.4f}  {tag:>4}  {cross_tag:>6}")

    # ==============================================================
    #  Pattern analysis
    # ==============================================================
    print()
    print("=" * 110)
    print("  PATTERN ANALYSIS")
    print("=" * 110)
    print()
    print(f"  Survived: {n_pass}/17")
    print(f"  Resonance crossings: {n_crossed}/17")
    print()

    # Failed BPs
    failed = [r for r in results if not r["pass_all"]]
    if failed:
        print("  FAILED BPs:")
        for r in failed:
            reason = []
            if not r["pass_30"]:
                reason.append(f"σ/m(30)={r['sm30_new']:.3f} {'<0.5' if r['sm30_new']<0.5 else '>10'}")
            if not r["pass_1k"]:
                reason.append(f"σ/m(1000)={r['sm1k_new']:.4f} >0.47")
            print(f"    {r['label']}: {', '.join(reason)}  (λ: {r['lam_orig']:.2f}→{r['lam_new']:.2f})")
    else:
        print("  ALL 17 BPs PASSED!")
    print()

    # Pattern by λ
    print("  σ/m(30) change vs λ:")
    for r in sorted(results, key=lambda x: x["lam_new"]):
        bar = "#" * max(1, int(abs(r["pct_30"]) / 2))
        sign = "▼" if r["pct_30"] < 0 else "▲"
        tag = " ✗" if not r["pass_all"] else ""
        print(f"    λ={r['lam_new']:5.2f}  {sign}{bar}  {r['pct_30']:+.1f}%{tag}")
    print()

    # Pattern by resonance proximity
    print("  Distance to nearest Yukawa resonance:")
    for r in sorted(results, key=lambda x: x["res_dist"]):
        tag = " ✗" if not r["pass_all"] else ""
        print(f"    {r['label']:>4}: λ_new={r['lam_new']:.2f}  "
              f"nearest R{r['res_n']}(λ_c={r['res_lam']:.2f})  "
              f"dist={r['res_dist']:.2f}{tag}")
    print()

    # Relic density pattern
    relic_good = [r for r in results if 0.115 <= r["omega"] <= 0.125]
    print(f"  KT relic density 0.115–0.125:  {len(relic_good)}/17 pass")
    for r in results:
        tag = "✓" if 0.115 <= r["omega"] <= 0.125 else " "
        print(f"    {r['label']:>4}: Ωh²_KT = {r['omega']:.4f}  {tag}")
    print()

    # Summary
    print("=" * 110)
    print("  SUMMARY")
    print("=" * 110)
    print(f"  θ = 19.47° decomposition: {n_pass}/17 BPs survive SIDM cuts")
    print(f"  σ/m(30) shifts: {min(r['pct_30'] for r in results):+.1f}% to {max(r['pct_30'] for r in results):+.1f}%")
    print(f"  Resonance crossings: {n_crossed}/17")
    if failed:
        lams = [r["lam_new"] for r in failed]
        print(f"  Failed λ_new range: {min(lams):.2f}–{max(lams):.2f}")
        passing = [r for r in results if r["pass_all"]]
        if passing:
            print(f"  Passing λ_new range: {min(r['lam_new'] for r in passing):.2f}–{max(r['lam_new'] for r in passing):.2f}")
    print()

if __name__ == "__main__":
    main()
