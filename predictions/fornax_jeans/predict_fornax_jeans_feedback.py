#!/usr/bin/env python3
"""
predictions/fornax_jeans/predict_fornax_jeans_feedback.py
=========================================================
Fornax dSph — Jeans Analysis with SIDM + Baryonic Feedback (SN heating)

Extends predict_fornax_jeans.py by incorporating supernova-driven
baryonic feedback that creates a DM core *independently* of SIDM.

Physics:
  1. Start with NFW halo
  2. Apply coreNFW transform (Read, Agertz & Collins 2016, MNRAS 459, 2573):
       M_cNFW(<r) = M_NFW(<r) × f^n(r)
       f(r) = tanh(r / r_c),  r_c = η × r_half
       n = tanh(κ × M_star / M_halo)   [0 → NFW, 1 → max core]
     SN outflows produce rapid potential fluctuations that heat DM
  3. On top of cNFW, apply SIDM isothermal coring (Kaplinghat+2016):
       find r_1 where σ/m × ρ_cNFW(r_1) × v_rel × t_age = 1
  4. Solve spherical Jeans equation (β = 0) with combined potential
  5. Compare 4 scenarios: NFW / feedback-only / SIDM-only / SIDM+feedback

Key question: Does BP1 (σ/m = 0.5 cm²/g) + feedback match Fornax data
              better than MAP (σ/m = 1.35) alone?

References:
  - Read, Agertz & Collins 2016, MNRAS 459, 2573 (coreNFW profile)
  - Di Cintio et al. 2014, MNRAS 441, 2986 (M_*/M_halo → core size)
  - Peñarrubia et al. 2012, ApJ 759, L42 (feedback + DM core interplay)
  - Pontzen & Governato 2012, MNRAS 421, 3464 (SN cusp-core transform)
  - Sameie et al. 2018, MNRAS 479, 359 (SIDM + baryons N-body)
  - Burger et al. 2022, MNRAS 513, 3458 (SIDM + feedback synergy)
  - Walker et al. 2009, ApJ 704, 1274 (Fornax kinematics)
"""
import sys, os, math, json
import numpy as np

# ---------- path bootstrap ----------
_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
# ------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Physical constants ──
G_KPC_MSUN = 4.302e-6          # G [kpc (km/s)² / M_sun]
MSUN_G     = 1.989e33          # g
KPC_CM     = 3.086e21          # cm
GYR_S      = 3.156e16          # s
KM_S_CM_S  = 1e5               # cm/s


# ════════════════════════════════════════════════════════════════════
#  NFW profile
# ════════════════════════════════════════════════════════════════════

def nfw_params(M200, c200):
    """Return (rho_s [M_sun/kpc³], r_s [kpc], R200 [kpc])."""
    h = GC.cosmological_constants()["h_hubble"]
    rho_crit = GC.cosmological_constants()["rho_crit_Msun_kpc3"]  # M_sun / kpc³
    R200 = (3 * M200 / (4 * math.pi * 200 * rho_crit))**(1.0 / 3.0)
    r_s = R200 / c200
    gc = math.log(1 + c200) - c200 / (1 + c200)
    rho_s = M200 / (4 * math.pi * r_s**3 * gc)
    return rho_s, r_s, R200


def rho_nfw(r, rho_s, r_s):
    """NFW density [M_sun/kpc³]."""
    x = r / r_s
    return rho_s / (x * (1 + x)**2)


def M_nfw(r, rho_s, r_s):
    """NFW enclosed mass [M_sun]."""
    x = r / r_s
    return 4 * math.pi * rho_s * r_s**3 * (math.log(1 + x) - x / (1 + x))


# ════════════════════════════════════════════════════════════════════
#  coreNFW profile — Read, Agertz & Collins 2016
# ════════════════════════════════════════════════════════════════════

def core_nfw_n(M_star, M_halo, kappa=80.0):
    """
    Core formation efficiency parameter n ∈ [0, 1].
    Read+2016 eq. 11: n = tanh(κ × M_star / M_halo)
    κ ~ 60–100 calibrated on FIRE/NIHAO sims.
    n → 1 when M_*/M_halo ~ 10⁻² (Fornax-like), n → 0 for UFDs.
    """
    return math.tanh(kappa * M_star / M_halo)


def M_core_nfw(r, rho_s, r_s, r_c, n):
    """
    coreNFW enclosed mass (Read+2016 eq. 9):
      M_cNFW(<r) = M_NFW(<r) × f^n(r)
      f(r) = tanh(r / r_c)
    """
    f = math.tanh(r / r_c)
    return M_nfw(r, rho_s, r_s) * f**n


def rho_core_nfw(r, rho_s, r_s, r_c, n, dr_frac=1e-4):
    """
    coreNFW density from ρ = (1/4πr²) dM/dr, computed numerically.
    """
    dr = max(r * dr_frac, 1e-6)
    M_plus = M_core_nfw(r + dr, rho_s, r_s, r_c, n)
    M_minus = M_core_nfw(r - dr if r > dr else 0.0, rho_s, r_s, r_c, n)
    dMdr = (M_plus - M_minus) / (2 * dr if r > dr else dr + r)
    return max(dMdr / (4 * math.pi * r**2), 0.0)


# ════════════════════════════════════════════════════════════════════
#  SIDM coring on top of an arbitrary base profile
# ════════════════════════════════════════════════════════════════════

def sidm_matching_on_profile(r_arr, rho_base_arr, sigma_over_m,
                              sigma_v_km_s, t_age_Gyr):
    """
    Find r_1 [kpc] where N_scatter = σ/m × ρ_base(r_1) × v_rel × t_age = 1
    on an arbitrary base density profile (NFW or coreNFW).
    Returns (r_1, rho_0 = rho_base(r_1)).
    """
    rho_conv = MSUN_G / KPC_CM**3
    v_rel = sigma_v_km_s * math.sqrt(2) * KM_S_CM_S
    t_age = t_age_Gyr * GYR_S

    # Walk from outside inward to find where N_scat crosses 1
    for i in range(len(r_arr) - 1, 0, -1):
        rho_cgs = rho_base_arr[i] * rho_conv
        n_scat = sigma_over_m * rho_cgs * v_rel * t_age
        if n_scat >= 1.0:
            # Interpolate between this point and the previous
            rho_cgs_prev = rho_base_arr[i + 1] * rho_conv if i + 1 < len(r_arr) else 0
            n_prev = sigma_over_m * rho_cgs_prev * v_rel * t_age
            if n_prev < 1.0 and n_scat > 1.0:
                # Linear interpolation in log space
                frac = (1.0 - n_prev) / (n_scat - n_prev) if (n_scat - n_prev) > 0 else 0.5
                r1 = r_arr[i + 1] + frac * (r_arr[i] - r_arr[i + 1]) if i + 1 < len(r_arr) else r_arr[i]
            else:
                r1 = r_arr[i]
            rho0 = np.interp(r1, r_arr, rho_base_arr)
            return r1, rho0
    # If never reaches 1 → no SIDM core
    return 0.0, rho_base_arr[0]


def build_combined_profile(r_arr, rho_base_arr, M_base_arr,
                            sigma_over_m, sigma_v_km_s, t_age_Gyr):
    """
    Build combined (feedback + SIDM) profile:
      - Outside r_1: base profile (NFW or coreNFW)
      - Inside r_1: King-like isothermal core ρ₀/(1 + (r/r_1)²)
    Returns (rho_combined, M_combined, r_core, rho_0).
    """
    if sigma_over_m < 1e-10:
        return rho_base_arr.copy(), M_base_arr.copy(), 0.0, rho_base_arr[0]

    r1, rho0 = sidm_matching_on_profile(
        r_arr, rho_base_arr, sigma_over_m, sigma_v_km_s, t_age_Gyr)

    if r1 < r_arr[1]:
        return rho_base_arr.copy(), M_base_arr.copy(), 0.0, rho_base_arr[0]

    rho_comb = np.zeros_like(r_arr)
    for i, r in enumerate(r_arr):
        if r <= r1:
            rho_comb[i] = rho0 / (1.0 + (r / r1)**2)
        else:
            rho_comb[i] = rho_base_arr[i]

    # Numerical enclosed mass
    M_comb = np.zeros_like(r_arr)
    for i in range(1, len(r_arr)):
        dr = r_arr[i] - r_arr[i - 1]
        rmid = 0.5 * (r_arr[i] + r_arr[i - 1])
        rhomid = 0.5 * (rho_comb[i] + rho_comb[i - 1])
        M_comb[i] = M_comb[i - 1] + 4 * math.pi * rmid**2 * rhomid * dr

    return rho_comb, M_comb, r1, rho0


# ════════════════════════════════════════════════════════════════════
#  Stellar (Plummer) profile
# ════════════════════════════════════════════════════════════════════

def rho_plummer(r, M_star, a):
    return (3 * M_star) / (4 * math.pi * a**3) * (1 + (r / a)**2)**(-2.5)


def M_plummer(r, M_star, a):
    return M_star * r**3 / (r**2 + a**2)**1.5


# ════════════════════════════════════════════════════════════════════
#  Jeans equation solver (isotropic β = 0)
# ════════════════════════════════════════════════════════════════════

def solve_jeans_isotropic(r_arr, rho_star_arr, M_tot_arr):
    """σ_r²(r) = (1/ρ_★) ∫_r^∞ ρ_★ G M_tot / r'² dr'. Returns (km/s)²."""
    integrand = rho_star_arr * G_KPC_MSUN * M_tot_arr / r_arr**2
    cum_integral = np.zeros(len(r_arr))
    for i in range(len(r_arr) - 2, -1, -1):
        dr = r_arr[i + 1] - r_arr[i]
        cum_integral[i] = cum_integral[i + 1] + \
                           0.5 * (integrand[i] + integrand[i + 1]) * dr
    sigma_r2 = np.zeros_like(r_arr)
    mask = rho_star_arr > 0
    sigma_r2[mask] = cum_integral[mask] / rho_star_arr[mask]
    return sigma_r2


def project_sigma_los(R_proj, r_arr, rho_star_arr, sigma_r2_arr):
    """Abel projection → σ_los(R) [km/s]."""
    sigma_los = np.zeros(len(R_proj))
    for j, R in enumerate(R_proj):
        mask = r_arr > R * 1.001
        r_sel = r_arr[mask]
        rho_sel = rho_star_arr[mask]
        sr2_sel = sigma_r2_arr[mask]
        if len(r_sel) < 3:
            continue
        integrand = rho_sel * sr2_sel * r_sel / np.sqrt(r_sel**2 - R**2)
        I_num = np.trapezoid(integrand, r_sel)
        integrand_sigma = rho_sel * r_sel / np.sqrt(r_sel**2 - R**2)
        Sigma_star = 2.0 * np.trapezoid(integrand_sigma, r_sel)
        if Sigma_star > 0:
            sigma_los[j] = math.sqrt(2.0 * I_num / Sigma_star)
    return sigma_los


def chi2_calc(slos_model, R_proj, R_obs, s_obs, s_err):
    """Interpolate model σ_los at observed R, return χ²."""
    slos_at_obs = np.interp(R_obs, R_proj, slos_model)
    return np.sum(((s_obs - slos_at_obs) / s_err)**2), slos_at_obs


# ════════════════════════════════════════════════════════════════════
#  Main
# ════════════════════════════════════════════════════════════════════

def main():
    cfg = load_config(__file__)
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = GC.benchmarks_from_labels(cfg['benchmark_labels'])
    fnx = GC.fornax_halo()
    obs = GC.walker2009_fornax()

    M200   = fnx['M200_Msun']
    c200   = fnx['c200']
    sv     = fnx['sigma_v_km_s']
    r_half = fnx['r_half_pc'] / 1000.0     # kpc
    M_star = fnx['M_star_Msun']
    t_age  = fnx['t_age_Gyr']

    rho_s, r_s, R200 = nfw_params(M200, c200)

    # Observational data
    R_obs = np.array(obs['R_pc']) / 1000.0
    s_obs = np.array(obs['sigma_los_km_s'])
    s_err = np.array(obs['sigma_err_km_s'])
    ndof = len(s_obs) - 1

    # Radial grid
    Nr = 3000
    r_arr = np.logspace(np.log10(0.001), np.log10(20.0), Nr)
    R_proj = np.linspace(0.020, 2.5, 120)

    # Stellar arrays
    rho_star_arr = np.array([rho_plummer(r, M_star, r_half) for r in r_arr])
    M_star_arr   = np.array([M_plummer(r, M_star, r_half) for r in r_arr])

    # ── coreNFW parameters (Read+2016) ──
    ETA = 1.75        # r_c = η × r_half (Read+2016 eq. 10)
    r_c = ETA * r_half
    n_fb = core_nfw_n(M_star, M200)
    X_ratio = M_star / M200

    print("=" * 96)
    print("  Fornax dSph — SIDM + Baryonic Feedback (SN heating) Jeans Analysis")
    print("=" * 96)
    print()
    print(f"  Halo:     M200 = {M200:.2e} M_sun, c200 = {c200}, R200 = {R200:.2f} kpc")
    print(f"  NFW:      rho_s = {rho_s:.3e} M_sun/kpc³, r_s = {r_s:.3f} kpc")
    print(f"  Stars:    M_★ = {M_star:.1e} M_sun, r_half = {r_half*1000:.0f} pc")
    print(f"  Feedback: M_★/M_halo = {X_ratio:.2e}, n = {n_fb:.3f}, "
          f"r_c = {r_c*1000:.0f} pc")
    print(f"  (n → 0: pure NFW; n → 1: max SN core; "
          f"Fornax has enough SFH for n = {n_fb:.2f})")
    print()

    # ── Build base profiles ──
    # 1. Pure NFW
    rho_nfw_arr = np.array([rho_nfw(r, rho_s, r_s) for r in r_arr])
    M_nfw_arr   = np.array([M_nfw(r, rho_s, r_s) for r in r_arr])

    # 2. coreNFW (feedback only, no SIDM)
    rho_cnfw_arr = np.array([rho_core_nfw(r, rho_s, r_s, r_c, n_fb) for r in r_arr])
    M_cnfw_arr   = np.array([M_core_nfw(r, rho_s, r_s, r_c, n_fb) for r in r_arr])

    # Feedback-only core size (where ρ_cNFW drops to half of central value)
    rho0_cnfw = rho_cnfw_arr[rho_cnfw_arr > 0][0] if rho_cnfw_arr[0] > 0 else rho_cnfw_arr[1]
    r_core_fb_idx = np.searchsorted(-rho_cnfw_arr, -rho0_cnfw * 0.5)
    r_core_fb = r_arr[min(r_core_fb_idx, len(r_arr) - 1)] * 1000  # pc

    # ════════════════════════════════════════════════════════════════
    #  Run 4 scenarios for each benchmark
    # ════════════════════════════════════════════════════════════════

    # --- Scenario A: Pure NFW ---
    sr2 = solve_jeans_isotropic(r_arr, rho_star_arr, M_nfw_arr + M_star_arr)
    slos_nfw = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2)
    chi2_nfw, slos_nfw_obs = chi2_calc(slos_nfw, R_proj, R_obs, s_obs, s_err)

    # --- Scenario B: Feedback only (coreNFW, no SIDM) ---
    sr2 = solve_jeans_isotropic(r_arr, rho_star_arr, M_cnfw_arr + M_star_arr)
    slos_fb = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2)
    chi2_fb, slos_fb_obs = chi2_calc(slos_fb, R_proj, R_obs, s_obs, s_err)

    print(f"  {'A) NFW (baseline)':32s}: χ²/dof = {chi2_nfw:6.1f}/{ndof}")
    print(f"    σ_los: {np.array2string(slos_nfw_obs, precision=1, separator=', ')}")
    print(f"  {'B) Feedback only (coreNFW)':32s}: χ²/dof = {chi2_fb:6.1f}/{ndof}  "
          f"[r_core_fb ≈ {r_core_fb:.0f} pc, n = {n_fb:.2f}]")
    print(f"    σ_los: {np.array2string(slos_fb_obs, precision=1, separator=', ')}")
    print()

    # --- Scenarios C (SIDM only) and D (SIDM + feedback) per benchmark ---
    all_results = []

    colors_bp = {'BP1': '#2196F3', 'BP9': '#4CAF50', 'MAP': '#E91E63', 'MAP_relic': '#9C27B0'}

    for bp in bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        lam   = alpha * m_chi / m_phi
        v_rel = sv * math.sqrt(2)
        sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_rel)

        print(f"  ─── {label}: m_χ={m_chi} GeV, m_φ={bp['m_phi_MeV']} MeV, "
              f"α={alpha:.3e}, λ={lam:.2f}, σ/m={sigma_m:.3f} cm²/g ───")

        # C) SIDM only (on NFW base)
        rho_sidm, M_sidm, r1_sidm, rho0_sidm = build_combined_profile(
            r_arr, rho_nfw_arr, M_nfw_arr, sigma_m, sv, t_age)
        sr2 = solve_jeans_isotropic(r_arr, rho_star_arr, M_sidm + M_star_arr)
        slos_sidm = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2)
        chi2_sidm, slos_sidm_obs = chi2_calc(slos_sidm, R_proj, R_obs, s_obs, s_err)

        # D) SIDM + feedback (on coreNFW base)
        rho_both, M_both, r1_both, rho0_both = build_combined_profile(
            r_arr, rho_cnfw_arr, M_cnfw_arr, sigma_m, sv, t_age)
        sr2 = solve_jeans_isotropic(r_arr, rho_star_arr, M_both + M_star_arr)
        slos_both = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2)
        chi2_both, slos_both_obs = chi2_calc(slos_both, R_proj, R_obs, s_obs, s_err)

        # Effective combined core: where ρ drops to half-central
        rho0_eff = rho_both[rho_both > 0][0] if rho_both[0] > 0 else rho_both[1]
        idx_eff = np.searchsorted(-rho_both, -rho0_eff * 0.5)
        r_core_eff = r_arr[min(idx_eff, len(r_arr) - 1)] * 1000

        print(f"    C) SIDM only:     r_core = {r1_sidm*1000:6.0f} pc  "
              f"χ²/dof = {chi2_sidm:6.1f}/{ndof}")
        print(f"       σ_los: {np.array2string(slos_sidm_obs, precision=1, separator=', ')}")
        print(f"    D) SIDM+feedback: r_core = {r_core_eff:6.0f} pc  "
              f"χ²/dof = {chi2_both:6.1f}/{ndof}")
        print(f"       σ_los: {np.array2string(slos_both_obs, precision=1, separator=', ')}")

        improvement = (chi2_sidm - chi2_both) / chi2_sidm * 100
        print(f"       Δχ² = {chi2_sidm - chi2_both:+.1f}  "
              f"({improvement:+.0f}% from adding feedback)")
        print()

        all_results.append({
            'label': label, 'sigma_m': sigma_m, 'lam': lam,
            'r_core_sidm_pc': r1_sidm * 1000,
            'r_core_eff_pc': r_core_eff,
            'chi2_sidm': chi2_sidm, 'chi2_both': chi2_both,
            'slos_sidm': slos_sidm, 'slos_both': slos_both,
            'slos_sidm_obs': slos_sidm_obs, 'slos_both_obs': slos_both_obs,
            'rho_sidm': rho_sidm, 'rho_both': rho_both,
        })

    # ════════════════════════════════════════════════════════════════
    #  Summary table
    # ════════════════════════════════════════════════════════════════
    print()
    print("  ┌──────────────────────────┬──────────┬──────────┬──────────┐")
    print("  │ Scenario                 │ r_core   │  χ²/dof  │ χ²/dof   │")
    print("  │                          │   [pc]   │          │  (per)   │")
    print("  ├──────────────────────────┼──────────┼──────────┼──────────┤")
    print(f"  │ A) NFW (baseline)        │ {'0':>8s} │"
          f" {chi2_nfw:6.1f}/{ndof} │ {chi2_nfw/ndof:6.2f}    │")
    print(f"  │ B) Feedback only         │ {r_core_fb:>6.0f}   │"
          f" {chi2_fb:6.1f}/{ndof} │ {chi2_fb/ndof:6.2f}    │")
    for res in all_results:
        print(f"  │ C) {res['label']:3s} SIDM only       │"
              f" {res['r_core_sidm_pc']:>6.0f}   │"
              f" {res['chi2_sidm']:6.1f}/{ndof} │"
              f" {res['chi2_sidm']/ndof:6.2f}    │")
        tag = "★" if res['chi2_both'] == min(r['chi2_both'] for r in all_results) else " "
        print(f"  │ D) {res['label']:3s} SIDM+feedback {tag}  │"
              f" {res['r_core_eff_pc']:>6.0f}   │"
              f" {res['chi2_both']:6.1f}/{ndof} │"
              f" {res['chi2_both']/ndof:6.2f}    │")
    print("  └──────────────────────────┴──────────┴──────────┴──────────┘")
    print()

    # Best result
    best_both = min(all_results, key=lambda x: x['chi2_both'])
    best_sidm = min(all_results, key=lambda x: x['chi2_sidm'])
    print(f"  Best SIDM-only:     {best_sidm['label']} "
          f"(χ²/dof = {best_sidm['chi2_sidm']:.1f}/{ndof} = "
          f"{best_sidm['chi2_sidm']/ndof:.2f})")
    print(f"  Best SIDM+feedback: {best_both['label']} "
          f"(χ²/dof = {best_both['chi2_both']:.1f}/{ndof} = "
          f"{best_both['chi2_both']/ndof:.2f})")
    print()

    # Key physics verdict
    bp1_res = next(r for r in all_results if r['label'] == 'BP1')
    map_res = next(r for r in all_results if r['label'] == 'MAP')
    print("  ── Physical Verdict ──")
    print(f"  BP1 SIDM-only χ²/dof = {bp1_res['chi2_sidm']/ndof:.2f}  →  "
          f"BP1 + feedback = {bp1_res['chi2_both']/ndof:.2f}")
    print(f"  MAP SIDM-only χ²/dof = {map_res['chi2_sidm']/ndof:.2f}  →  "
          f"MAP + feedback = {map_res['chi2_both']/ndof:.2f}")
    if bp1_res['chi2_both'] < map_res['chi2_sidm']:
        print("  ✓ BP1 + feedback BEATS MAP alone — "
              "baryonic heating resolves the tension!")
    else:
        print("  ○ BP1 + feedback does NOT beat MAP alone.")
    print()

    # ════════════════════════════════════════════════════════════════
    #  Plot: 3-panel figure
    # ════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))

    # Panel 1: Density profiles (BP1 focus)
    ax = axes[0]
    ax.loglog(r_arr * 1000, rho_nfw_arr, 'k:', lw=1.5, label='NFW', alpha=0.5)
    ax.loglog(r_arr * 1000, rho_cnfw_arr, 'k--', lw=1.5,
              label=f'coreNFW (feedback, n={n_fb:.2f})', alpha=0.7)
    for res in all_results:
        c = colors_bp[res['label']]
        ax.loglog(r_arr * 1000, res['rho_sidm'], color=c, ls='-', lw=1.5,
                  label=f"{res['label']} SIDM only", alpha=0.6)
        ax.loglog(r_arr * 1000, res['rho_both'], color=c, ls='--', lw=2.5,
                  label=f"{res['label']} SIDM+fb")
    ax.axvline(r_half * 1000, color='gray', ls=':', lw=0.8, alpha=0.5)
    ax.set_xlabel('$r$ [pc]', fontsize=12)
    ax.set_ylabel(r'$\rho_{\rm DM}$ [M$_\odot$ kpc$^{-3}$]', fontsize=12)
    ax.set_title('DM Density Profiles', fontsize=12)
    ax.legend(fontsize=7, loc='lower left', ncol=1)
    ax.set_xlim(1, 1e4)
    ax.set_ylim(1e4, 1e9)

    # Panel 2: σ_los comparison (all scenarios, all BPs)
    ax = axes[1]
    ax.errorbar(R_obs * 1000, s_obs, yerr=s_err,
                fmt='ko', ms=6, capsize=3, label='Walker+2009', zorder=10)
    ax.plot(R_proj * 1000, slos_nfw, 'k:', lw=1.5, label='NFW', alpha=0.5)
    ax.plot(R_proj * 1000, slos_fb, 'k--', lw=1.5,
            label='Feedback only', alpha=0.7)
    for res in all_results:
        c = colors_bp[res['label']]
        ax.plot(R_proj * 1000, res['slos_sidm'], color=c, ls='-', lw=1.5,
                label=f"{res['label']} SIDM ($\\chi^2$={res['chi2_sidm']:.0f})",
                alpha=0.6)
        ax.plot(R_proj * 1000, res['slos_both'], color=c, ls='--', lw=2.5,
                label=f"{res['label']} +fb ($\\chi^2$={res['chi2_both']:.0f})")
    ax.axvline(r_half * 1000, color='gray', ls=':', lw=0.8, alpha=0.5)
    ax.set_xlabel('Projected radius $R$ [pc]', fontsize=12)
    ax.set_ylabel(r'$\sigma_{\rm los}$ [km/s]', fontsize=12)
    ax.set_title(r'Stellar Velocity Dispersion ($\beta = 0$)', fontsize=12)
    ax.legend(fontsize=7, loc='upper right', ncol=1)
    ax.set_xlim(0, 2500)
    ax.set_ylim(0, 18)

    # Panel 3: χ² bar chart
    ax = axes[2]
    labels_bar = ['NFW', 'Feedback\nonly']
    chi2_bars = [chi2_nfw, chi2_fb]
    bar_colors = ['gray', '#9E9E9E']
    for res in all_results:
        labels_bar.extend([f"{res['label']}\nSIDM", f"{res['label']}\n+fb"])
        chi2_bars.extend([res['chi2_sidm'], res['chi2_both']])
        bar_colors.extend([colors_bp[res['label']], colors_bp[res['label']]])

    alphas = [0.4, 0.6] + [0.5, 1.0] * len(all_results)
    bars = ax.bar(range(len(chi2_bars)), chi2_bars, color=bar_colors, alpha=1)
    for i, (bar, a) in enumerate(zip(bars, alphas)):
        bar.set_alpha(a)
    ax.set_xticks(range(len(labels_bar)))
    ax.set_xticklabels(labels_bar, fontsize=8)
    ax.set_ylabel(r'$\chi^2$', fontsize=12)
    ax.set_title(r'Goodness of Fit ($\chi^2$, 8 dof)', fontsize=12)
    ax.axhline(ndof, color='green', ls='--', lw=1, alpha=0.7, label=f'dof = {ndof}')
    ax.axhline(2 * ndof, color='orange', ls='--', lw=1, alpha=0.5, label=f'2×dof = {2*ndof}')
    for i, v in enumerate(chi2_bars):
        ax.text(i, v + 2, f'{v:.0f}', ha='center', fontsize=7)
    ax.legend(fontsize=8)
    ax.set_ylim(0, max(chi2_bars) * 1.2)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'fornax_jeans_feedback.png')
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # ════════════════════════════════════════════════════════════════
    #  Sensitivity: scan over κ (feedback strength parameter)
    # ════════════════════════════════════════════════════════════════
    print()
    print("  ── Sensitivity to κ (feedback strength) ──")
    print(f"  {'κ':>5s} {'n_fb':>6s} {'BP1 SIDM':>10s} {'BP1 +fb':>10s} "
          f"{'MAP SIDM':>10s} {'MAP +fb':>10s}")

    bp1_bp = next(b for b in bps if b['label'] == 'BP1')
    map_bp = next(b for b in bps if b['label'] == 'MAP')
    sm_bp1 = sigma_T_vpm(bp1_bp['m_chi_GeV'], bp1_bp['m_phi_MeV']/1000,
                          bp1_bp['alpha'], sv * math.sqrt(2))
    sm_map = sigma_T_vpm(map_bp['m_chi_GeV'], map_bp['m_phi_MeV']/1000,
                          map_bp['alpha'], sv * math.sqrt(2))

    for kappa in [20, 40, 60, 80, 100, 150, 200]:
        n_test = core_nfw_n(M_star, M200, kappa)
        rho_cn = np.array([rho_core_nfw(r, rho_s, r_s, r_c, n_test)
                           for r in r_arr])
        M_cn   = np.array([M_core_nfw(r, rho_s, r_s, r_c, n_test)
                           for r in r_arr])

        # BP1 + feedback
        _, M_b1, _, _ = build_combined_profile(
            r_arr, rho_cn, M_cn, sm_bp1, sv, t_age)
        sr2 = solve_jeans_isotropic(r_arr, rho_star_arr, M_b1 + M_star_arr)
        slos = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2)
        chi2_b1, _ = chi2_calc(slos, R_proj, R_obs, s_obs, s_err)

        # BP1 SIDM only (constant)
        chi2_b1_sidm = bp1_res['chi2_sidm']

        # MAP + feedback
        _, M_bm, _, _ = build_combined_profile(
            r_arr, rho_cn, M_cn, sm_map, sv, t_age)
        sr2 = solve_jeans_isotropic(r_arr, rho_star_arr, M_bm + M_star_arr)
        slos = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2)
        chi2_bm, _ = chi2_calc(slos, R_proj, R_obs, s_obs, s_err)

        # MAP SIDM only (constant)
        chi2_bm_sidm = map_res['chi2_sidm']

        print(f"  {kappa:>5d} {n_test:>6.3f} {chi2_b1_sidm:>10.1f} "
              f"{chi2_b1:>10.1f} {chi2_bm_sidm:>10.1f} {chi2_bm:>10.1f}")

    print()


if __name__ == '__main__':
    main()
