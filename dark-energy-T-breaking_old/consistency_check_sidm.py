#!/usr/bin/env python3
# ⚠️ SUPERSEDED — This script applies (8/9) reduction to α_CSV,
# which is already α_s. See test_alpha_convention.py (Test 12).
# Correct results: all 17 BPs pass with α_s = α_CSV (no reduction).
# Buggy output archived in archived_buggy/consistency_check_sidm_BUGGY.txt
"""
dark-energy-T-breaking/consistency_check_sidm.py
=================================================
Consistency check: dark axion θ-decomposition vs SIDM benchmark points.

QUESTION
--------
If the Yukawa coupling y decomposes as:
    y → y_s = y cosθ  (scalar)  +  y_p = y sinθ  (pseudoscalar)
    with  θ = arcsin(1/3) ≈ 19.47°
do the existing SIDM benchmark points survive?

PHYSICS
-------
SIDM scattering (NR, v ~ 30 km/s):
  - Scalar exchange: V(r) = -α_s e^{-m_φ r}/r,  α_s = α cos²θ = (8/9)α
  - Pseudoscalar exchange: V_p ~ O(v²/c²) → negligible at SIDM velocities
  ⟹ effective SIDM coupling reduced by factor 8/9

Relic density (freeze-out, v ~ 0.3c):
  - s-wave dominated by scalar: ⟨σv⟩ = π α_s² / (4 m_χ²)
  - pseudoscalar contribution is p-wave suppressed
  ⟹ relic ⟨σv⟩ reduced by factor (8/9)² = 64/81

SCENARIOS TESTED
----------------
A. "Direct decomposition": α_s = (8/9) α_original
   → SIDM weakens, relic Ωh² increases (less annihilation)

B. "Rescaled total": α_total = (9/8) α_original so that α_s = α_original
   → SIDM unchanged, relic uses α_s = α_original → Ωh² changes

C. "Relic-matched": find α_total so that Ωh² stays at 0.120
   → Check if the new α_s still passes SIDM cuts

CODE REUSE (from Secluded-Majorana-SIDM)
-----------------------------------------
  VPM solver:       core/v22_raw_scan.py  (sigma_T_vpm + helpers)
  Relic density:    core/v27_boltzmann_relic.py  (kolb_turner_swave, Y_to_omega_h2)
  Constants & BPs:  data/global_config.json
  
  All copied functions are UNMODIFIED except where marked "# DARK AXION EDIT".
"""
import sys, math
import numpy as np

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

try:
    from numba import jit
except ImportError:
    print("WARNING: numba not installed. VPM solver will be SLOW.")
    def jit(**kwargs):
        def decorator(func):
            return func
        return decorator


# ==============================================================
#  DARK AXION ANGLE  (the ONE new parameter)
# ==============================================================
THETA = math.asin(1.0 / 3.0)           # 19.47°
COS2_THETA = math.cos(THETA)**2        # 8/9
SIN2_THETA = math.sin(THETA)**2        # 1/9


# ==============================================================
#  Constants (from Secluded-Majorana-SIDM/data/global_config.json)
# ==============================================================
GEV2_TO_CM2 = 3.8938e-28
GEV_IN_G    = 1.78266e-24
C_KM_S      = 299792.458
M_PL        = 1.220890e19    # GeV (Planck mass)
S_0         = 2891.2         # cm⁻³ (entropy density today)
RHO_CRIT_H2 = 1.0539e-5     # GeV/cm³


# ==============================================================
#  Benchmark points (from global_config.json)
# ==============================================================
BENCHMARKS = {
    "BP1":       {"m_chi": 20.69,  "m_phi_MeV": 11.34, "alpha": 1.048e-3},
    "BP9":       {"m_chi": 37.93,  "m_phi_MeV": 12.98, "alpha": 1.858e-3},
    "BP16":      {"m_chi": 20.70,  "m_phi_MeV":  9.91, "alpha": 1.048e-3},
    "MAP":       {"m_chi": 94.07,  "m_phi_MeV": 11.10, "alpha": 5.734e-3},
    "MAP_relic": {"m_chi": 94.07,  "m_phi_MeV": 11.10, "alpha": 4.523e-3},
}

# SIDM cuts (from global_config.json, sourced from literature)
SIGMA_M_30_LO  = 0.5     # cm²/g  (Elbert+2015)
SIGMA_M_30_HI  = 10.0    # cm²/g  (KTY16)
SIGMA_M_1000_HI = 0.47   # cm²/g  (Harvey+2015)


# ==============================================================
#  VPM SOLVER — copied verbatim from core/v22_raw_scan.py
# ==============================================================

@jit(nopython=True, cache=True)
def sph_jn_numba(l, z):
    """Spherical Bessel j_l(z) via upward recurrence."""
    if z < 1e-30:
        return 1.0 if l == 0 else 0.0
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
    """Spherical Bessel y_l(z) via upward recurrence."""
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
    """RHS of VPM ODE: d(delta_l)/dx."""
    if x < 1e-20: return 0.0
    z = kappa * x
    if z < 1e-20: return 0.0
    jl = sph_jn_numba(l, z)
    nl = sph_yn_numba(l, z)
    j_hat = z * jl
    n_hat = -z * nl
    cd = math.cos(delta)
    sd = math.sin(delta)
    bracket = j_hat * cd - n_hat * sd
    if not math.isfinite(bracket): return 0.0
    pot = lam * math.exp(-x) / (kappa * x)
    val = pot * bracket * bracket
    return val if math.isfinite(val) else 0.0


@jit(nopython=True, cache=True)
def vpm_phase_shift(l, kappa, lam, x_max=50.0, N_steps=4000):
    """Phase shift delta_l via RK4."""
    if lam < 1e-30 or kappa < 1e-30: return 0.0
    x_min = max(1e-5, 0.05 / (kappa + 0.01))
    if l > 0:
        x_barrier = l / kappa
        if x_barrier > x_min:
            x_min = x_barrier
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
    """σ_T/m [cm²/g] for identical Majorana fermions (full VPM).
    
    Copied from Secluded-Majorana-SIDM/core/v22_raw_scan.py.
    Unmodified — the dark axion coupling change enters through
    the 'alpha' argument, not inside this function.
    """
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15: return 0.0

    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000

    l_max_hard = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)

    sigma_sum = 0.0
    peak_contrib = 0.0
    n_small = 0
    for l in range(l_max_hard + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        weight = 1.0 if l % 2 == 0 else 3.0
        contrib = weight * (2*l + 1) * math.sin(delta)**2
        sigma_sum += contrib
        if contrib > peak_contrib:
            peak_contrib = contrib
        if peak_contrib > 0.0 and contrib / peak_contrib < 1e-4:
            n_small += 1
            if n_small >= 5: break
        else:
            n_small = 0

    sigma_GeV2 = 2.0 * math.pi * sigma_sum / (k * k)
    sigma_cm2 = sigma_GeV2 * GEV2_TO_CM2
    return sigma_cm2 / (m_chi * GEV_IN_G)


# ==============================================================
#  RELIC DENSITY — copied from core/v27_boltzmann_relic.py
# ==============================================================

_G_STAR_TABLE = np.array([
    [1e4,    106.75, 106.75],
    [200,    106.75, 106.75],
    [80,      86.25,  86.25],
    [10,      86.25,  86.25],
    [1,       75.75,  75.75],
    [0.3,     61.75,  61.75],
    [0.2,     17.25,  17.25],
    [0.15,    14.25,  14.25],
    [0.1,     10.75,  10.75],
    [0.01,    10.75,  10.75],
    [0.001,   10.75,  10.75],
    [0.0005,  10.75,  10.75],
    [0.0001,   3.36,   3.91],
    [1e-5,     3.36,   3.91],
    [1e-8,     3.36,   3.91],
])
_LOG_T = np.log(_G_STAR_TABLE[:, 0])
_G_RHO = _G_STAR_TABLE[:, 1]
_G_S   = _G_STAR_TABLE[:, 2]


def g_star_rho(T):
    logT = math.log(T) if T > 0 else -50
    return float(np.interp(logT, _LOG_T[::-1], _G_RHO[::-1]))


def g_star_S(T):
    logT = math.log(T) if T > 0 else -50
    return float(np.interp(logT, _LOG_T[::-1], _G_S[::-1]))


def sigma_v_swave(alpha_d, m_chi):
    """⟨σv⟩ = πα²/(4m_χ²) for Majorana χχ→φφ (s-wave, scalar only)."""
    return math.pi * alpha_d**2 / (4.0 * m_chi**2)


def kolb_turner_swave(m_chi, sv0, g_chi=2):
    """Kolb-Turner analytic freeze-out for s-wave.
    Returns (x_fo, Y_inf).
    Copied from v27_boltzmann_relic.py — unmodified.
    """
    x_fo = 20.0
    for _ in range(50):
        T_fo = m_chi / x_fo
        gr = g_star_rho(T_fo)
        arg = 0.0764 * g_chi * M_PL * m_chi * sv0 / math.sqrt(gr * x_fo)
        if arg <= 0: break
        x_fo_new = math.log(arg) - 0.5 * math.log(x_fo)
        if abs(x_fo_new - x_fo) < 0.001:
            x_fo = x_fo_new
            break
        x_fo = 0.5 * x_fo + 0.5 * x_fo_new

    T_fo = m_chi / x_fo
    gs = g_star_S(T_fo)
    gr = g_star_rho(T_fo)
    g_eff = gs / math.sqrt(gr)
    lambda_val = math.sqrt(math.pi / 45.0) * g_eff * M_PL * m_chi
    Y_inf = x_fo / (lambda_val * sv0)
    return x_fo, Y_inf


def Y_to_omega_h2(Y_inf, m_chi):
    """Convert Y_∞ to Ωh². Copied from v27_boltzmann_relic.py."""
    return m_chi * Y_inf * S_0 / RHO_CRIT_H2


# ==============================================================
#  CONSISTENCY CHECK
# ==============================================================

def compute_bp(label, m_chi, m_phi_GeV, alpha, tag=""):
    """Compute σ_T/m and Ωh² for a single benchmark point."""
    sv0 = sigma_v_swave(alpha, m_chi)
    x_fo, Y_inf = kolb_turner_swave(m_chi, sv0)
    omega = Y_to_omega_h2(Y_inf, m_chi)
    lam = alpha * m_chi / m_phi_GeV
    
    sm30   = sigma_T_vpm(m_chi, m_phi_GeV, alpha, 30.0)
    sm1000 = sigma_T_vpm(m_chi, m_phi_GeV, alpha, 1000.0)
    
    pass_30  = SIGMA_M_30_LO <= sm30 <= SIGMA_M_30_HI
    pass_1k  = sm1000 < SIGMA_M_1000_HI
    
    return {
        "label": label, "tag": tag,
        "alpha": alpha, "lambda": lam,
        "sv0": sv0, "x_fo": x_fo, "omega_h2": omega,
        "sm30": sm30, "sm1000": sm1000,
        "pass_30": pass_30, "pass_1000": pass_1k,
        "pass_all": pass_30 and pass_1k,
    }


def print_results_table(results, title):
    """Print a comparison table."""
    print()
    print("=" * 100)
    print(f"  {title}")
    print("=" * 100)
    hdr = (f"  {'BP':>10}  {'α':>10}  {'λ':>6}  {'Ωh²':>8}  {'x_fo':>5}  "
           f"{'σ/m(30)':>10}  {'σ/m(1k)':>10}  {'SIDM?':>5}  {'relic?':>6}")
    print(hdr)
    print("  " + "-" * 96)
    
    for r in results:
        relic_ok = "✓" if 0.115 <= r["omega_h2"] <= 0.125 else "✗"
        sidm_ok = "✓" if r["pass_all"] else "✗"
        print(f"  {r['label']:>10}  {r['alpha']:.3e}  {r['lambda']:5.2f}  "
              f"{r['omega_h2']:.4f}  {r['x_fo']:5.1f}  "
              f"{r['sm30']:10.4f}  {r['sm1000']:10.4f}  "
              f"{sidm_ok:>5}  {relic_ok:>6}")


def main():
    print("=" * 100)
    print("  CONSISTENCY CHECK: Dark Axion θ-Decomposition vs SIDM Benchmarks")
    print("=" * 100)
    print()
    print(f"  Dark axion angle θ = arcsin(1/3) = {math.degrees(THETA):.2f}°")
    print(f"  cos²θ = {COS2_THETA:.6f}  (= 8/9 = {8/9:.6f})")
    print(f"  sin²θ = {SIN2_THETA:.6f}  (= 1/9 = {1/9:.6f})")
    print()
    print("  SIDM scattering: only SCALAR exchange matters (pseudoscalar is v²-suppressed)")
    print(f"  → effective α_SIDM = α_s = α × cos²θ = α × {COS2_THETA:.4f}")
    print()
    print("  Relic density: s-wave from scalar coupling only")
    print(f"  → ⟨σv⟩ = π α_s² / (4m²) = π (αcos²θ)² / (4m²)")
    print(f"  → reduction factor: cos⁴θ = {COS2_THETA**2:.6f} (= 64/81 = {64/81:.6f})")
    print()

    # Warmup numba JIT
    print("  [JIT warmup] ", end="", flush=True)
    _ = sigma_T_vpm(20.0, 0.01, 1e-3, 30.0)
    print("done.")
    print()

    # ------------------------------------------------------------------
    #  SCENARIO A: Original benchmarks (baseline)
    # ------------------------------------------------------------------
    print("  Computing Scenario A: ORIGINAL benchmarks (pure scalar)...")
    original = []
    for label, bp in BENCHMARKS.items():
        m_phi_GeV = bp["m_phi_MeV"] * 1e-3
        r = compute_bp(label, bp["m_chi"], m_phi_GeV, bp["alpha"], "original")
        original.append(r)
    print_results_table(original, "SCENARIO A — ORIGINAL (pure scalar, α as-is)")

    # ------------------------------------------------------------------
    #  SCENARIO B: Direct decomposition  α → α_s = (8/9)α
    # ------------------------------------------------------------------
    print()
    print("  Computing Scenario B: α_SIDM = (8/9)α, ⟨σv⟩_relic uses α_s = (8/9)α ...")
    decomposed = []
    for label, bp in BENCHMARKS.items():
        m_phi_GeV = bp["m_phi_MeV"] * 1e-3
        alpha_s = bp["alpha"] * COS2_THETA  # DARK AXION EDIT: α → (8/9)α
        r = compute_bp(label, bp["m_chi"], m_phi_GeV, alpha_s, "decomposed")
        decomposed.append(r)
    print_results_table(decomposed, 
        "SCENARIO B — DIRECT DECOMPOSITION (α_s = 8α/9 for both SIDM and relic)")

    # ------------------------------------------------------------------
    #  SCENARIO C: Rescaled total  α_total = (9/8)α  so α_s = α_original
    # ------------------------------------------------------------------
    print()
    print("  Computing Scenario C: α_total = (9/8)α so α_s = α_original...")
    rescaled = []
    for label, bp in BENCHMARKS.items():
        m_phi_GeV = bp["m_phi_MeV"] * 1e-3
        # SIDM uses α_s = α_total × cos²θ = (9/8)α × (8/9) = α  → same as original
        # Relic uses α_s = (9/8)α × cos²θ = α  → same ⟨σv⟩ → same Ωh²
        # NOTE: This scenario preserves EVERYTHING — it's just a rescaling
        # The total coupling is 12.5% larger, but SIDM and relic both see α_s = α
        alpha_total = bp["alpha"] * 9.0 / 8.0
        alpha_s_for_sidm = alpha_total * COS2_THETA  # = α_original
        r = compute_bp(label, bp["m_chi"], m_phi_GeV, alpha_s_for_sidm, "rescaled")
        r["alpha_total"] = alpha_total
        r["alpha_p"] = alpha_total * SIN2_THETA
        rescaled.append(r)
    print_results_table(rescaled, 
        "SCENARIO C — RESCALED (α_total = 9α/8, α_s = α_original, α_p = α/8)")

    # ------------------------------------------------------------------
    #  SCENARIO D: Find α_total so Ωh² = 0.120 with decomposed coupling
    #              then check if SIDM survives
    # ------------------------------------------------------------------
    print()
    print("  Computing Scenario D: Find α_total giving Ωh² = 0.120, check SIDM...")
    print()
    relic_matched = []
    for label, bp in BENCHMARKS.items():
        m_chi = bp["m_chi"]
        m_phi_GeV = bp["m_phi_MeV"] * 1e-3
        
        # Bisect to find α_s that gives Ωh² = 0.120  (s-wave relic uses α_s)
        # In the original model: α_original gives Ωh² ≈ 0.120
        # With decomposition: α_s < α_original → Ωh² > 0.120
        # Need larger α_s: bisect in [α_original, 2×α_original]
        def omega_for_alpha_s(a_s):
            sv = sigma_v_swave(a_s, m_chi)
            _, Y = kolb_turner_swave(m_chi, sv)
            return Y_to_omega_h2(Y, m_chi)
        
        a_lo, a_hi = bp["alpha"] * 0.5, bp["alpha"] * 3.0
        for _ in range(60):
            a_mid = (a_lo + a_hi) / 2.0
            om = omega_for_alpha_s(a_mid)
            if om > 0.120:
                a_lo = a_mid  # need more annihilation → larger α
            else:
                a_hi = a_mid
        
        alpha_s_matched = (a_lo + a_hi) / 2.0
        alpha_total_matched = alpha_s_matched / COS2_THETA  # α_total = α_s / cos²θ
        
        # Now check SIDM with this α_s
        r = compute_bp(label, m_chi, m_phi_GeV, alpha_s_matched, "relic-matched")
        r["alpha_total"] = alpha_total_matched
        r["alpha_s"] = alpha_s_matched
        r["alpha_p"] = alpha_total_matched * SIN2_THETA
        r["ratio_to_original"] = alpha_s_matched / bp["alpha"]
        relic_matched.append(r)
    
    print_results_table(relic_matched, 
        "SCENARIO D — RELIC-MATCHED (α_s tuned so Ωh² = 0.120, then check SIDM)")
    
    # Additional detail for Scenario D
    print()
    print("  Scenario D detail — coupling adjustments:")
    print(f"  {'BP':>10}  {'α_orig':>10}  {'α_s(new)':>10}  {'ratio':>7}  "
          f"{'α_total':>10}  {'α_p':>10}  {'λ_new':>6}")
    print("  " + "-" * 80)
    for r in relic_matched:
        bp = BENCHMARKS[r["label"]]
        print(f"  {r['label']:>10}  {bp['alpha']:.3e}  {r['alpha_s']:.3e}  "
              f"{r['ratio_to_original']:.4f}  {r['alpha_total']:.3e}  "
              f"{r['alpha_p']:.3e}  {r['lambda']:5.2f}")

    # ------------------------------------------------------------------
    #  SUMMARY
    # ------------------------------------------------------------------
    print()
    print()
    print("=" * 100)
    print("  SUMMARY")
    print("=" * 100)
    print()
    print("  θ = arcsin(1/3) = 19.47°  →  cos²θ = 8/9,  sin²θ = 1/9")
    print()
    
    # Scenario B summary
    b_pass = sum(1 for r in decomposed if r["pass_all"])
    print(f"  Scenario B (direct α → 8α/9):  {b_pass}/{len(decomposed)} BPs pass SIDM")
    for r in decomposed:
        bp = BENCHMARKS[r["label"]]
        pct_change_30  = (r["sm30"] / [o for o in original if o["label"]==r["label"]][0]["sm30"] - 1) * 100
        pct_change_1k = (r["sm1000"] / [o for o in original if o["label"]==r["label"]][0]["sm1000"] - 1) * 100
        print(f"    {r['label']:>10}: σ/m(30) {pct_change_30:+.1f}%  σ/m(1000) {pct_change_1k:+.1f}%  "
              f"Ωh² = {r['omega_h2']:.4f}")

    print()
    
    # Scenario D summary
    d_pass = sum(1 for r in relic_matched if r["pass_all"])
    print(f"  Scenario D (relic-matched α_s):  {d_pass}/{len(relic_matched)} BPs pass SIDM")
    for r in relic_matched:
        print(f"    {r['label']:>10}: α_s/α_orig = {r['ratio_to_original']:.4f}  "
              f"σ/m(30) = {r['sm30']:.3f}  σ/m(1000) = {r['sm1000']:.4f}  "
              f"{'✓' if r['pass_all'] else '✗'}")

    print()
    print("  KEY QUESTION: Does the 8/9 reduction in effective coupling")
    print("  break or preserve the SIDM phenomenology?")
    print()


if __name__ == "__main__":
    main()
