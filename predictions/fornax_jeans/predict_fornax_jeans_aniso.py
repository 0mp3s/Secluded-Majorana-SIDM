#!/usr/bin/env python3
"""
predictions/fornax_jeans/predict_fornax_jeans_aniso.py
======================================================
Fornax dSph — Full Jeans Analysis: SIDM + Feedback + Anisotropy

Extends the feedback script by adding Osipkov-Merritt velocity anisotropy:

    β(r) = r² / (r² + r_a²)

which gives β → 0 (isotropic) at r ≪ r_a and β → 1 (radial) at r ≫ r_a.

The anisotropic Jeans equation (Binney & Tremaine 2008, eq. 4.215):

    d/dr [ρ_★ σ_r²] + (2β/r) ρ_★ σ_r² = -ρ_★ G M_tot / r²

For Osipkov-Merritt, the solution uses the substitution Q = r² + r_a²:

    σ_r²(r) = [1/(ρ_★ Q)] ∫_r^∞ ρ_★ Q G M_tot / r'² dr'

    where Q(r) = 1 + (r/r_a)²

The projected σ_los includes the anisotropy kernel (Mamon & Łokas 2005):

    σ_los²(R) = (2/Σ_★) ∫_R^∞ K(r/R, r_a/R) ρ_★ σ_r² r dr / √(r²-R²)

    K(u, s) = (1 + s⁻²) × [1 - √(1-u⁻²) if u>1 else π/2]
              minus the tangential correction term.

We implement the exact Mamon & Łokas (2005) kernel for Osipkov-Merritt.

Scenarios tested (all with best-fit r_a):
  A) NFW + β(r)
  B) Feedback only + β(r)
  C) SIDM only + β(r)          [BP1, BP9, MAP]
  D) SIDM + Feedback + β(r)    [BP1, BP9, MAP]

Also: scan r_a to find best-fit for each scenario.

References:
  - Osipkov 1979, PAZh 5, 77; Merritt 1985, AJ 90, 1027
  - Mamon & Łokas 2005, MNRAS 363, 705 (anisotropic Jeans projection)
  - Read et al. 2019, MNRAS 484, 1401 (GravSphere; β profiles for dSphs)
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

sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)  # JIT warmup

_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Physical constants ──
G_KPC_MSUN = 4.302e-6          # G [kpc (km/s)² / M_sun]
MSUN_G     = 1.989e33          # g
KPC_CM     = 3.086e21          # cm
GYR_S      = 3.156e16          # s
KM_S_CM_S  = 1e5               # cm/s


# ════════════════════════════════════════════════════════════════
#  NFW
# ════════════════════════════════════════════════════════════════

def nfw_params(M200, c200):
    """Compute NFW scale density rho_s [M_sun/kpc^3], scale radius r_s [kpc], and R200 [kpc]."""
    h = 0.674
    rho_crit = 277.5 * h**2
    R200 = (3 * M200 / (4 * math.pi * 200 * rho_crit))**(1.0/3.0)
    r_s = R200 / c200
    gc = math.log(1 + c200) - c200 / (1 + c200)
    rho_s = M200 / (4 * math.pi * r_s**3 * gc)
    return rho_s, r_s, R200

def rho_nfw(r, rho_s, r_s):
    """NFW density profile rho(r) [M_sun/kpc^3]."""
    x = r / r_s
    return rho_s / (x * (1 + x)**2)

def M_nfw(r, rho_s, r_s):
    """NFW enclosed mass M(<r) [M_sun]."""
    x = r / r_s
    return 4 * math.pi * rho_s * r_s**3 * (math.log(1 + x) - x / (1 + x))


# ════════════════════════════════════════════════════════════════
#  coreNFW (Read+2016) — feedback
# ════════════════════════════════════════════════════════════════

def core_nfw_n(M_star, M_halo, kappa=80.0):
    """coreNFW feedback index n = tanh(kappa * M_star/M_halo) (Read+2016)."""
    return math.tanh(kappa * M_star / M_halo)

def M_core_nfw(r, rho_s, r_s, r_c, n):
    """coreNFW enclosed mass: M_NFW(r) * tanh(r/r_c)^n."""
    f = math.tanh(r / r_c)
    return M_nfw(r, rho_s, r_s) * f**n

def rho_core_nfw(r, rho_s, r_s, r_c, n, dr_frac=1e-4):
    """coreNFW density via numerical differentiation of M_core_nfw."""
    dr = max(r * dr_frac, 1e-6)
    M_plus  = M_core_nfw(r + dr, rho_s, r_s, r_c, n)
    M_minus = M_core_nfw(r - dr if r > dr else 0.0, rho_s, r_s, r_c, n)
    dMdr = (M_plus - M_minus) / (2 * dr if r > dr else dr + r)
    return max(dMdr / (4 * math.pi * r**2), 0.0)


# ════════════════════════════════════════════════════════════════
#  SIDM coring on arbitrary base profile
# ════════════════════════════════════════════════════════════════

def sidm_matching_on_profile(r_arr, rho_base_arr, sigma_over_m,
                              sigma_v_km_s, t_age_Gyr):
    """Find SIDM coring radius r1 where rho * (sigma/m) * v * t = 1."""
    rho_conv = MSUN_G / KPC_CM**3
    v_rel = sigma_v_km_s * math.sqrt(2) * KM_S_CM_S
    t_age = t_age_Gyr * GYR_S
    for i in range(len(r_arr) - 1, 0, -1):
        rho_cgs = rho_base_arr[i] * rho_conv
        n_scat = sigma_over_m * rho_cgs * v_rel * t_age
        if n_scat >= 1.0:
            rho_cgs_prev = rho_base_arr[i+1] * rho_conv if i+1 < len(r_arr) else 0
            n_prev = sigma_over_m * rho_cgs_prev * v_rel * t_age
            if n_prev < 1.0 and n_scat > 1.0:
                frac = (1.0 - n_prev) / (n_scat - n_prev) if (n_scat - n_prev) > 0 else 0.5
                r1 = r_arr[i+1] + frac * (r_arr[i] - r_arr[i+1]) if i+1 < len(r_arr) else r_arr[i]
            else:
                r1 = r_arr[i]
            rho0 = np.interp(r1, r_arr, rho_base_arr)
            return r1, rho0
    return 0.0, rho_base_arr[0]

def build_combined_profile(r_arr, rho_base_arr, M_base_arr,
                            sigma_over_m, sigma_v_km_s, t_age_Gyr):
    """Build SIDM-cored density/mass profile: isothermal core inside r1, NFW outside."""
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
    M_comb = np.zeros_like(r_arr)
    for i in range(1, len(r_arr)):
        dr = r_arr[i] - r_arr[i-1]
        rmid = 0.5 * (r_arr[i] + r_arr[i-1])
        rhomid = 0.5 * (rho_comb[i] + rho_comb[i-1])
        M_comb[i] = M_comb[i-1] + 4 * math.pi * rmid**2 * rhomid * dr
    return rho_comb, M_comb, r1, rho0


# ════════════════════════════════════════════════════════════════
#  Stellar (Plummer)
# ════════════════════════════════════════════════════════════════

def rho_plummer(r, M_star, a):
    """Plummer stellar density profile [M_sun/kpc^3]."""
    return (3 * M_star) / (4 * math.pi * a**3) * (1 + (r / a)**2)**(-2.5)

def M_plummer(r, M_star, a):
    """Plummer enclosed stellar mass [M_sun]."""
    return M_star * r**3 / (r**2 + a**2)**1.5


# ════════════════════════════════════════════════════════════════
#  Anisotropic Jeans equation: Osipkov-Merritt
# ════════════════════════════════════════════════════════════════

def solve_jeans_OM(r_arr, rho_star_arr, M_tot_arr, r_a):
    """
    Solve the anisotropic Jeans equation with Osipkov-Merritt β(r).
    β(r) = r²/(r² + r_a²).

    Solution (BT08 §4.3.2):
      σ_r²(r) = 1/(ρ_★ Q(r)) × ∫_r^∞ ρ_★(r') Q(r') G M_tot(r') / r'² dr'
    where Q(r) = 1 + (r/r_a)².

    If r_a → ∞: Q → 1, recovers isotropic case.
    Returns σ_r² in (km/s)².
    """
    Q = 1.0 + (r_arr / r_a)**2
    integrand = rho_star_arr * Q * G_KPC_MSUN * M_tot_arr / r_arr**2

    # Backward cumulative integral (trapezoid)
    cum = np.zeros(len(r_arr))
    for i in range(len(r_arr) - 2, -1, -1):
        dr = r_arr[i+1] - r_arr[i]
        cum[i] = cum[i+1] + 0.5 * (integrand[i] + integrand[i+1]) * dr

    sigma_r2 = np.zeros_like(r_arr)
    mask = (rho_star_arr > 0) & (Q > 0)
    sigma_r2[mask] = cum[mask] / (rho_star_arr[mask] * Q[mask])
    return sigma_r2


def project_sigma_los_OM(R_proj, r_arr, rho_star_arr, sigma_r2_arr, r_a):
    """
    Abel projection for Osipkov-Merritt anisotropy.

    Mamon & Łokas (2005, MNRAS 363, 705), eq. 23:

      σ_los²(R) = (2/Σ_★) ∫_R^∞ K_OM(r, R, r_a) ρ_★ σ_r² dr

    where the OM kernel:
      K_OM(r, R, r_a) = [r / √(r²-R²)] × [1 - R²/r² × (r² + r_a²)/(R² + r_a²)]

    The factor (r² + r_a²)/(R² + r_a²) = Q(r)/Q(R) handles the
    anisotropy contribution to the tangential dispersion.
    """
    sigma_los = np.zeros(len(R_proj))
    for j, R in enumerate(R_proj):
        mask = r_arr > R * 1.001
        r_sel = r_arr[mask]
        rho_sel = rho_star_arr[mask]
        sr2_sel = sigma_r2_arr[mask]
        if len(r_sel) < 3:
            continue

        # Anisotropy factor: 1 - β(r) R²/r²
        # β(r) = r²/(r²+r_a²)  →  1 - R²/(r²+r_a²)
        factor = 1.0 - R**2 / (r_sel**2 + r_a**2)
        K = factor * r_sel / np.sqrt(r_sel**2 - R**2)

        integrand = K * rho_sel * sr2_sel
        I_num = np.trapezoid(integrand, r_sel)

        # Surface brightness denominator
        integrand_Sigma = rho_sel * r_sel / np.sqrt(r_sel**2 - R**2)
        Sigma_star = 2.0 * np.trapezoid(integrand_Sigma, r_sel)

        if Sigma_star > 0:
            val = 2.0 * I_num / Sigma_star
            sigma_los[j] = math.sqrt(max(val, 0.0))
    return sigma_los


def chi2_calc(slos_model, R_proj, R_obs, s_obs, s_err):
    slos_at_obs = np.interp(R_obs, R_proj, slos_model)
    return np.sum(((s_obs - slos_at_obs) / s_err)**2), slos_at_obs


# ════════════════════════════════════════════════════════════════
#  Main
# ════════════════════════════════════════════════════════════════

def main():
    cfg = load_config(__file__)
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps  = GC.benchmarks_from_labels(cfg['benchmark_labels'])
    fnx  = GC.fornax_halo()
    obs  = GC.walker2009_fornax()

    M200   = fnx['M200_Msun']
    c200   = fnx['c200']
    sv     = fnx['sigma_v_km_s']
    r_half = fnx['r_half_pc'] / 1000.0
    M_star = fnx['M_star_Msun']
    t_age  = fnx['t_age_Gyr']

    rho_s, r_s, R200 = nfw_params(M200, c200)

    R_obs = np.array(obs['R_pc']) / 1000.0
    s_obs = np.array(obs['sigma_los_km_s'])
    s_err = np.array(obs['sigma_err_km_s'])
    ndof_base = len(s_obs)

    # Radial grid
    Nr = 3000
    r_arr = np.logspace(np.log10(0.001), np.log10(20.0), Nr)
    R_proj = np.linspace(0.020, 2.5, 120)

    # Stars
    rho_star_arr = np.array([rho_plummer(r, M_star, r_half) for r in r_arr])
    M_star_arr   = np.array([M_plummer(r, M_star, r_half) for r in r_arr])

    # coreNFW parameters
    ETA = 1.75
    r_c = ETA * r_half
    n_fb = core_nfw_n(M_star, M200)

    # NFW arrays
    rho_nfw_arr = np.array([rho_nfw(r, rho_s, r_s) for r in r_arr])
    M_nfw_arr   = np.array([M_nfw(r, rho_s, r_s) for r in r_arr])

    # coreNFW arrays (feedback)
    rho_cnfw_arr = np.array([rho_core_nfw(r, rho_s, r_s, r_c, n_fb) for r in r_arr])
    M_cnfw_arr   = np.array([M_core_nfw(r, rho_s, r_s, r_c, n_fb) for r in r_arr])

    print("=" * 96)
    print("  Fornax dSph — Anisotropic Jeans: SIDM + Feedback + Osipkov-Merritt β(r)")
    print("=" * 96)
    print()
    print(f"  Halo: M200={M200:.2e}, c200={c200}, ρ_s={rho_s:.3e} M_sun/kpc³, r_s={r_s:.3f} kpc")
    print(f"  Stars: M_★={M_star:.1e}, r_half={r_half*1000:.0f} pc, n_fb={n_fb:.3f}")
    print(f"  β(r) = r²/(r² + r_a²) — scan r_a ∈ [0.3, 3.0] kpc")
    print()

    # ── r_a scan values ──
    r_a_scan = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                          1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0])
    # r_a=50 is effectively isotropic (β<0.001 everywhere relevant)

    # ────────────────────────────────────────────────────────────
    #  Build DM profiles for all scenarios
    # ────────────────────────────────────────────────────────────
    scenarios = []

    # A) NFW
    scenarios.append({
        'tag': 'NFW', 'label': 'A) NFW',
        'M_dm': M_nfw_arr, 'color': 'gray', 'ls': ':'
    })

    # B) Feedback only (coreNFW)
    scenarios.append({
        'tag': 'FB', 'label': 'B) Feedback only',
        'M_dm': M_cnfw_arr, 'color': '#607D8B', 'ls': '--'
    })

    # Per-BP scenarios
    colors_bp = {'BP1': '#2196F3', 'BP9': '#4CAF50', 'MAP': '#E91E63', 'MAP_relic': '#9C27B0'}
    for bp in bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        v_rel = sv * math.sqrt(2)
        sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_rel)

        # C) SIDM only (on NFW)
        _, M_sidm, r1_s, _ = build_combined_profile(
            r_arr, rho_nfw_arr, M_nfw_arr, sigma_m, sv, t_age)
        scenarios.append({
            'tag': f'{label}_SIDM', 'label': f'C) {label} SIDM',
            'M_dm': M_sidm, 'color': colors_bp[label], 'ls': '-',
            'sigma_m': sigma_m, 'r1': r1_s
        })

        # D) SIDM + feedback (on coreNFW)
        _, M_both, r1_b, _ = build_combined_profile(
            r_arr, rho_cnfw_arr, M_cnfw_arr, sigma_m, sv, t_age)
        scenarios.append({
            'tag': f'{label}_SIDM_FB', 'label': f'D) {label} SIDM+fb',
            'M_dm': M_both, 'color': colors_bp[label], 'ls': '--',
            'sigma_m': sigma_m, 'r1': r1_b
        })

    # ────────────────────────────────────────────────────────────
    #  For each scenario: scan r_a, find best-fit
    # ────────────────────────────────────────────────────────────
    print(f"  {'Scenario':28s} │ β=0 χ²  │ best r_a │ best χ² │ χ²/dof │ Δχ²")
    print("  " + "─" * 28 + "─┼─────────┼──────────┼─────────┼────────┼──────")

    best_results = []

    for sc in scenarios:
        M_tot = sc['M_dm'] + M_star_arr

        # Isotropic baseline (r_a → ∞)
        sr2_iso = solve_jeans_OM(r_arr, rho_star_arr, M_tot, 50.0)
        slos_iso = project_sigma_los_OM(R_proj, r_arr, rho_star_arr, sr2_iso, 50.0)
        chi2_iso, _ = chi2_calc(slos_iso, R_proj, R_obs, s_obs, s_err)

        # Scan r_a
        best_chi2 = chi2_iso
        best_ra = 50.0
        best_slos = slos_iso
        chi2_vs_ra = []

        for r_a in r_a_scan:
            sr2 = solve_jeans_OM(r_arr, rho_star_arr, M_tot, r_a)
            slos = project_sigma_los_OM(R_proj, r_arr, rho_star_arr, sr2, r_a)
            c2, _ = chi2_calc(slos, R_proj, R_obs, s_obs, s_err)
            chi2_vs_ra.append(c2)
            if c2 < best_chi2:
                best_chi2 = c2
                best_ra = r_a
                best_slos = slos.copy()

        ndof = ndof_base - 2  # 2 fit params: DM profile + r_a
        _, best_slos_obs = chi2_calc(best_slos, R_proj, R_obs, s_obs, s_err)
        delta = chi2_iso - best_chi2

        print(f"  {sc['label']:28s} │ {chi2_iso:7.1f} │ {best_ra:6.2f}   │"
              f" {best_chi2:7.1f} │ {best_chi2/ndof:6.2f} │ {delta:+5.1f}")

        best_results.append({
            **sc,
            'chi2_iso': chi2_iso,
            'best_ra': best_ra,
            'best_chi2': best_chi2,
            'best_slos': best_slos,
            'best_slos_obs': best_slos_obs,
            'chi2_vs_ra': chi2_vs_ra,
            'ndof': ndof,
        })

    print()

    # ────────────────────────────────────────────────────────────
    #  Highlight key comparison
    # ────────────────────────────────────────────────────────────
    bp1_fb = next(r for r in best_results if r['tag'] == 'BP1_SIDM_FB')
    bp1_sidm = next(r for r in best_results if r['tag'] == 'BP1_SIDM')
    map_sidm = next(r for r in best_results if r['tag'] == 'MAP_SIDM')
    map_fb = next(r for r in best_results if r['tag'] == 'MAP_SIDM_FB')
    nfw_res = next(r for r in best_results if r['tag'] == 'NFW')

    print("  ── Key Comparisons ──")
    print()
    print(f"  NFW isotropic:               χ² = {nfw_res['chi2_iso']:.1f}")
    print(f"  NFW + β(r_a={nfw_res['best_ra']:.1f} kpc):    χ² = {nfw_res['best_chi2']:.1f}")
    print()
    print(f"  BP1 SIDM only (β=0):         χ² = {bp1_sidm['chi2_iso']:.1f}")
    print(f"  BP1 SIDM+fb (β=0):           χ² = {bp1_fb['chi2_iso']:.1f}")
    print(f"  BP1 SIDM+fb + β(r_a={bp1_fb['best_ra']:.1f}): χ² = {bp1_fb['best_chi2']:.1f}  ★")
    print()
    print(f"  MAP SIDM only (β=0):         χ² = {map_sidm['chi2_iso']:.1f}")
    print(f"  MAP SIDM+fb + β(r_a={map_fb['best_ra']:.1f}): χ² = {map_fb['best_chi2']:.1f}")
    print()

    # σ_los detail for MAP SIDM+fb (best overall)
    print(f"  σ_los at observed radii [km/s] — MAP SIDM+fb + β(r_a={map_fb['best_ra']:.1f}):")
    print(f"    R [pc]:  {np.array2string(R_obs*1000, precision=0, separator=', ')}")
    print(f"    Data:    {np.array2string(s_obs, precision=1, separator=', ')}")
    _, map_fb_obs = chi2_calc(map_fb['best_slos'], R_proj, R_obs, s_obs, s_err)
    print(f"    Model:   {np.array2string(map_fb_obs, precision=1, separator=', ')}")
    print(f"    Error:   {np.array2string(s_err, precision=1, separator=', ')}")
    residuals_map = (s_obs - map_fb_obs) / s_err
    print(f"    Residual (σ): {np.array2string(residuals_map, precision=2, separator=', ')}")
    print()

    # σ_los detail for BP1 SIDM+fb (for comparison)
    print(f"  σ_los at observed radii [km/s] — BP1 SIDM+fb (β=0, r_a=∞):")
    print(f"    Model:   {np.array2string(bp1_fb['best_slos_obs'], precision=1, separator=', ')}")
    residuals_bp1 = (s_obs - bp1_fb['best_slos_obs']) / s_err
    print(f"    Residual (σ): {np.array2string(residuals_bp1, precision=2, separator=', ')}")
    print()

    # χ² scan printout for key scenarios
    print("  ── χ² vs r_a scan ──")
    hdr = f"  {'r_a [kpc]':>10s}"
    for tag in ('BP1_SIDM_FB', 'MAP_SIDM_FB', 'NFW'):
        hdr += f" │ {tag:>14s}"
    print(hdr)
    for k, ra in enumerate(r_a_scan):
        line = f"  {ra:10.2f}"
        for tag in ('BP1_SIDM_FB', 'MAP_SIDM_FB', 'NFW'):
            res = next(r for r in best_results if r['tag'] == tag)
            line += f" │ {res['chi2_vs_ra'][k]:14.1f}"
        print(line)
    print()

    # Best overall
    overall_best = min(best_results, key=lambda x: x['best_chi2'])
    print(f"  ★ Overall best: {overall_best['label']}, "
          f"r_a = {overall_best['best_ra']:.2f} kpc, "
          f"χ²/dof = {overall_best['best_chi2']:.1f}/{overall_best['ndof']} = "
          f"{overall_best['best_chi2']/overall_best['ndof']:.2f}")
    verdict = "EXCELLENT" if overall_best['best_chi2']/overall_best['ndof'] < 2 else \
              "GOOD" if overall_best['best_chi2']/overall_best['ndof'] < 3 else "MARGINAL"
    print(f"    Verdict: {verdict}")
    print()

    # ════════════════════════════════════════════════════════════
    #  Plot 1: σ_los panel — best-fit β for each scenario
    # ════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))

    # Panel 1: σ_los with best-fit β
    ax = axes[0]
    ax.errorbar(R_obs*1000, s_obs, yerr=s_err, fmt='ko', ms=6, capsize=3,
                label='Walker+2009', zorder=10)
    for res in best_results:
        lbl = f"{res['label']} ($r_a$={res['best_ra']:.1f}, $\\chi^2$={res['best_chi2']:.0f})"
        ax.plot(R_proj*1000, res['best_slos'], color=res['color'],
                ls=res['ls'], lw=2 if 'FB' in res['tag'] else 1.5,
                label=lbl, alpha=0.9 if 'FB' in res['tag'] else 0.6)
    ax.axvline(r_half*1000, color='gray', ls=':', lw=0.8, alpha=0.5)
    ax.set_xlabel('Projected radius $R$ [pc]', fontsize=12)
    ax.set_ylabel(r'$\sigma_{\rm los}$ [km/s]', fontsize=12)
    ax.set_title(r'Best-fit $\sigma_{\rm los}$ with OM $\beta(r)$', fontsize=12)
    ax.legend(fontsize=6.5, loc='upper right', ncol=1)
    ax.set_xlim(0, 2500)
    ax.set_ylim(0, 18)

    # Panel 2: χ² vs r_a for key scenarios
    ax = axes[1]
    for res in best_results:
        if res['tag'] in ('NFW', 'FB', 'BP1_SIDM', 'BP1_SIDM_FB',
                          'MAP_SIDM', 'MAP_SIDM_FB'):
            ax.plot(r_a_scan, res['chi2_vs_ra'], color=res['color'],
                    ls=res['ls'], lw=2, marker='o', ms=3,
                    label=f"{res['label']}")
    ax.axhline(ndof_base - 2, color='green', ls='--', lw=1, alpha=0.5,
               label=f'dof = {ndof_base - 2}')
    ax.axhline(2*(ndof_base-2), color='orange', ls='--', lw=1, alpha=0.5,
               label=f'2×dof = {2*(ndof_base-2)}')
    ax.set_xlabel(r'Anisotropy radius $r_a$ [kpc]', fontsize=12)
    ax.set_ylabel(r'$\chi^2$', fontsize=12)
    ax.set_title(r'$\chi^2$ vs Anisotropy Radius', fontsize=12)
    ax.set_xscale('log')
    ax.set_xlim(0.25, 60)
    ax.legend(fontsize=7, loc='upper left', ncol=1)

    # Panel 3: β(r) profile + bar chart comparison
    ax = axes[2]
    labels_bar = []
    chi2_iso_bars = []
    chi2_aniso_bars = []
    bar_colors = []
    for res in best_results:
        labels_bar.append(res['tag'].replace('_', '\n'))
        chi2_iso_bars.append(res['chi2_iso'])
        chi2_aniso_bars.append(res['best_chi2'])
        bar_colors.append(res['color'])
    x_pos = np.arange(len(labels_bar))
    w = 0.35
    bars1 = ax.bar(x_pos - w/2, chi2_iso_bars, w, color=bar_colors, alpha=0.3,
                   label=r'$\beta = 0$')
    bars2 = ax.bar(x_pos + w/2, chi2_aniso_bars, w, color=bar_colors, alpha=0.9,
                   label=r'Best $\beta(r)$')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels_bar, fontsize=6.5)
    ax.set_ylabel(r'$\chi^2$', fontsize=12)
    ax.set_title(r'Isotropic vs Anisotropic $\chi^2$', fontsize=12)
    ax.axhline(ndof_base-2, color='green', ls='--', lw=1, alpha=0.5)
    for i, (v1, v2) in enumerate(zip(chi2_iso_bars, chi2_aniso_bars)):
        ax.text(i - w/2, v1 + 1.5, f'{v1:.0f}', ha='center', fontsize=6, alpha=0.5)
        ax.text(i + w/2, v2 + 1.5, f'{v2:.0f}', ha='center', fontsize=6.5, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(0, max(chi2_iso_bars) * 1.15)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'fornax_jeans_aniso.png')
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # ════════════════════════════════════════════════════════════
    #  Plot 2: β(r) profile illustration
    # ════════════════════════════════════════════════════════════
    fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))
    r_plot = np.linspace(0.01, 3.0, 200)  # kpc
    for r_a_val in [0.5, 0.8, 1.0, 1.5, 3.0]:
        beta_OM = r_plot**2 / (r_plot**2 + r_a_val**2)
        ax2.plot(r_plot*1000, beta_OM, lw=2, label=f'$r_a$ = {r_a_val} kpc')
    ax2.axvline(r_half*1000, color='gray', ls=':', lw=1, label=f'$r_{{half}}$ = {r_half*1000:.0f} pc')
    ax2.set_xlabel('$r$ [pc]', fontsize=12)
    ax2.set_ylabel(r'$\beta(r)$', fontsize=12)
    ax2.set_title('Osipkov-Merritt Anisotropy Profile', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.set_xlim(0, 3000)
    ax2.set_ylim(-0.05, 1.05)
    fig2.tight_layout()
    fig2_path = os.path.join(out_dir, 'beta_OM_profiles.png')
    fig2.savefig(fig2_path, dpi=150)
    plt.close(fig2)
    print(f"  Saved: {fig2_path}")


if __name__ == '__main__':
    main()
