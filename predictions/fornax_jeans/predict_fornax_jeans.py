#!/usr/bin/env python3
"""
predictions/fornax_jeans/predict_fornax_jeans.py
=================================================
§7.2 extended: Fornax dSph Stellar Velocity Dispersion Profile — Jeans Analysis

Physics:
  1. Compute SIDM density profile via Kaplinghat+2016 isothermal matching:
     - Find matching radius r_1 where σ/m × ρ_NFW(r_1) × v_rel × t_age = 1
     - Inside r_1: cored isothermal profile ρ ∝ 1/(1 + (r/r_c)²)
     - Outside r_1: NFW profile
  2. Stellar distribution: Plummer profile with r_half = 710 pc (Walker+2009)
  3. Solve spherical Jeans equation (isotropic, β=0):
       σ_r²(r) = (1/ρ_★) ∫_r^∞ ρ_★(r') GM(r')/r'² dr'
  4. Project to line-of-sight:
       σ_los²(R) = (2/Σ_★) ∫_R^∞ ρ_★ σ_r² r/√(r²−R²) dr
  5. Compare with Walker+2009 binned kinematic data

Model: Majorana χ + real scalar φ with (y_s + iy_p γ_5) coupling
VPM: σ_T/m(v) from core/v22_raw_scan.py with α = α_s (scattering coupling)

References:
  - Walker et al. 2009, ApJ 704, 1274 (Fornax kinematics, 2633 stars)
  - Walker & Peñarrubia 2011, ApJ 742, 20 (mass estimators)
  - Read et al. 2019, MNRAS 484, 1401 (halo parameters)
  - Kaplinghat, Tulin & Yu 2016, PRL 116, 041302 (SIDM isothermal matching)
  - Amorisco & Evans 2012, MNRAS 419, 184 (Fornax DM profile)
"""
import sys, os, math, json
import numpy as np
from scipy import integrate

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
G_KPC_MSUN = 4.302e-6          # G [kpc (km/s)² / M_sun]  (4.302e-3 pc → ÷1000)
MSUN_G     = 1.989e33          # g
KPC_CM     = 3.086e21          # cm
PC_CM      = 3.086e18          # cm
GYR_S      = 3.156e16          # s
KM_S_CM_S  = 1e5               # cm/s


# ════════════════════════════════════════════════════════════════════
#  NFW profile
# ════════════════════════════════════════════════════════════════════

def nfw_params(M200, c200):
    """Return (rho_s [M_sun/kpc³], r_s [kpc], R200 [kpc]) from M200, c200."""
    h = GC.cosmological_constants()["h_hubble"]
    rho_crit = GC.cosmological_constants()["rho_crit_Msun_kpc3"]  # M_sun / kpc³
    R200 = (3 * M200 / (4 * math.pi * 200 * rho_crit))**(1.0/3.0)
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
#  SIDM cored profile (Kaplinghat+2016 isothermal matching)
# ════════════════════════════════════════════════════════════════════

def sidm_matching_radius(rho_s, r_s, sigma_over_m, sigma_v_km_s, t_age_Gyr):
    """
    Find r_1 [kpc] where N_scatter = σ/m × ρ_NFW × v_rel × t_age = 1.
    """
    rho_conv = MSUN_G / KPC_CM**3        # M_sun/kpc³ → g/cm³
    v_rel = sigma_v_km_s * math.sqrt(2) * KM_S_CM_S   # cm/s
    t_age = t_age_Gyr * GYR_S                          # s

    r_lo, r_hi = 1e-4, 20.0       # kpc
    for _ in range(120):
        r_mid = math.sqrt(r_lo * r_hi)
        rho_cgs = rho_nfw(r_mid, rho_s, r_s) * rho_conv
        n_scat = sigma_over_m * rho_cgs * v_rel * t_age
        if n_scat > 1.0:
            r_lo = r_mid
        else:
            r_hi = r_mid
    return math.sqrt(r_lo * r_hi)


def rho_sidm(r, rho_s, r_s, r1, rho0):
    """
    SIDM density [M_sun/kpc³].
    Inside r_1: King-like core  ρ₀/(1 + (r/r_1)²)
    Outside r_1: NFW
    """
    if r <= r1:
        return rho0 / (1.0 + (r / r1)**2)
    return rho_nfw(r, rho_s, r_s)


def M_sidm_enclosed(r_arr, rho_s, r_s, r1, rho0):
    """Enclosed DM mass [M_sun] by numerical integration of the SIDM profile."""
    M = np.zeros(len(r_arr))
    for i in range(1, len(r_arr)):
        dr = r_arr[i] - r_arr[i-1]
        rmid = 0.5 * (r_arr[i] + r_arr[i-1])
        rhomid = rho_sidm(rmid, rho_s, r_s, r1, rho0)
        M[i] = M[i-1] + 4 * math.pi * rmid**2 * rhomid * dr
    return M


# ════════════════════════════════════════════════════════════════════
#  Stellar (Plummer) profile
# ════════════════════════════════════════════════════════════════════

def rho_plummer(r, M_star, a):
    """3D Plummer density [M_sun/kpc³].  a = projected half-light radius."""
    return (3 * M_star) / (4 * math.pi * a**3) * (1 + (r / a)**2)**(-2.5)


def M_plummer(r, M_star, a):
    """Plummer enclosed stellar mass [M_sun]."""
    return M_star * r**3 / (r**2 + a**2)**1.5


def sigma_plummer(R, M_star, a):
    """Plummer projected surface density [M_sun/kpc²]."""
    return M_star / (math.pi * a**2) * (1 + (R / a)**2)**(-2)


# ════════════════════════════════════════════════════════════════════
#  Jeans equation solver (isotropic β = 0)
# ════════════════════════════════════════════════════════════════════

def solve_jeans_isotropic(r_arr, rho_star_arr, M_tot_arr):
    """
    Solve σ_r²(r) = (1/ρ_★(r)) ∫_r^∞ ρ_★(r') G M_tot(r') / r'² dr'.

    Returns σ_r² in (km/s)² at each r.
    Uses G in [km²/s² kpc / M_sun], so r in kpc and M in M_sun give km²/s².
    """
    n = len(r_arr)
    sigma_r2 = np.zeros(n)

    # Integrate from outside in (backward cumulative trapezoid)
    integrand = rho_star_arr * G_KPC_MSUN * M_tot_arr / r_arr**2
    # Cumulative integral from r to r_max (reversed)
    cum_integral = np.zeros(n)
    for i in range(n - 2, -1, -1):
        dr = r_arr[i + 1] - r_arr[i]
        cum_integral[i] = cum_integral[i + 1] + 0.5 * (integrand[i] + integrand[i + 1]) * dr

    mask = rho_star_arr > 0
    sigma_r2[mask] = cum_integral[mask] / rho_star_arr[mask]
    return sigma_r2  # (km/s)²


def project_sigma_los(R_proj, r_arr, rho_star_arr, sigma_r2_arr):
    """
    Project σ_r²(r) to line-of-sight σ_los²(R) assuming isotropy (β=0):
      σ_los²(R) = (2/Σ_★(R)) ∫_R^∞ ρ_★(r) σ_r²(r)  r / √(r²−R²) dr

    R_proj: array of projected radii [kpc]
    Returns σ_los in km/s.
    """
    sigma_los = np.zeros(len(R_proj))

    for j, R in enumerate(R_proj):
        # Select r > R (avoid singularity at r = R)
        mask = r_arr > R * 1.001
        r_sel = r_arr[mask]
        rho_sel = rho_star_arr[mask]
        sr2_sel = sigma_r2_arr[mask]

        if len(r_sel) < 3:
            continue

        integrand = rho_sel * sr2_sel * r_sel / np.sqrt(r_sel**2 - R**2)
        I_num = np.trapezoid(integrand, r_sel)

        # Σ_★(R) by Abel projection
        integrand_sigma = rho_sel * r_sel / np.sqrt(r_sel**2 - R**2)
        Sigma_star = 2.0 * np.trapezoid(integrand_sigma, r_sel)

        if Sigma_star > 0:
            sigma_los[j] = math.sqrt(2.0 * I_num / Sigma_star)

    return sigma_los


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
    r_half = fnx['r_half_pc'] / 1000.0    # → kpc
    M_star = fnx['M_star_Msun']
    t_age  = fnx['t_age_Gyr']

    rho_s, r_s, R200 = nfw_params(M200, c200)

    print("=" * 90)
    print("  Fornax dSph — Stellar Velocity Dispersion from Jeans Equation (SIDM)")
    print("=" * 90)
    print()
    print(f"  Halo:  M200 = {M200:.2e} M_sun, c200 = {c200}, R200 = {R200:.2f} kpc")
    print(f"  NFW:   rho_s = {rho_s:.3e} M_sun/kpc³, r_s = {r_s:.3f} kpc")
    print(f"  Stars: M_★ = {M_star:.1e} M_sun, r_half = {r_half*1000:.0f} pc")
    print(f"  σ_v = {sv} km/s, t_age = {t_age} Gyr")
    print()

    # ── Radial grid ──
    # Fine grid from 1 pc to 20 kpc (extended to capture full integral)
    Nr = 3000
    r_arr = np.logspace(np.log10(0.001), np.log10(20.0), Nr)   # kpc

    # Stellar density on grid
    rho_star_arr = np.array([rho_plummer(r, M_star, r_half) for r in r_arr])

    # Observational data
    R_obs = np.array(obs['R_pc']) / 1000.0        # → kpc
    s_obs = np.array(obs['sigma_los_km_s'])
    s_err = np.array(obs['sigma_err_km_s'])

    # Projection radii for model curves
    R_proj = np.linspace(0.020, 2.5, 120)          # kpc

    # ── NFW baseline (no SIDM) ──
    M_nfw_arr = np.array([M_nfw(r, rho_s, r_s) for r in r_arr])
    M_star_arr = np.array([M_plummer(r, M_star, r_half) for r in r_arr])
    M_tot_nfw = M_nfw_arr + M_star_arr

    sr2_nfw = solve_jeans_isotropic(r_arr, rho_star_arr, M_tot_nfw)
    slos_nfw = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2_nfw)

    # χ² for NFW baseline
    slos_nfw_at_obs = np.interp(R_obs, R_proj, slos_nfw)
    chi2_nfw = np.sum(((s_obs - slos_nfw_at_obs) / s_err)**2)
    ndof = len(s_obs) - 1
    print(f"  {'NFW (no SIDM)':18s}: χ²/dof = {chi2_nfw:.1f}/{ndof}")
    print(f"    σ_los at data R: {np.array2string(slos_nfw_at_obs, precision=1, separator=', ')} km/s")

    # ── SIDM benchmarks ──
    results = []
    colors = {'BP1': '#2196F3', 'BP9': '#4CAF50', 'MAP': '#E91E63', 'MAP_relic': '#9C27B0'}
    styles = {'BP1': '-', 'BP9': '--', 'MAP': '-.', 'MAP_relic': ':'}

    for bp in bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0       # → GeV
        alpha = bp['alpha']
        lam   = alpha * m_chi / m_phi

        v_rel = sv * math.sqrt(2)
        sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_rel)

        # Core formation
        if sigma_m < 1e-10:
            r1 = 0.0
            rho0 = rho_nfw(1e-4, rho_s, r_s)
        else:
            r1 = sidm_matching_radius(rho_s, r_s, sigma_m, sv, t_age)
            rho0 = rho_nfw(r1, rho_s, r_s)

        # SIDM enclosed mass on grid
        M_dm_arr = M_sidm_enclosed(r_arr, rho_s, r_s, r1, rho0)
        M_tot = M_dm_arr + M_star_arr

        # Solve Jeans
        sr2 = solve_jeans_isotropic(r_arr, rho_star_arr, M_tot)
        slos = project_sigma_los(R_proj, r_arr, rho_star_arr, sr2)

        # χ²
        slos_at_obs = np.interp(R_obs, R_proj, slos)
        chi2 = np.sum(((s_obs - slos_at_obs) / s_err)**2)

        results.append({
            'label': label, 'lam': lam, 'sigma_m': sigma_m,
            'r_core_pc': r1 * 1000, 'rho0': rho0,
            'slos': slos, 'chi2': chi2,
            'slos_at_obs': slos_at_obs,
        })

        print(f"  {label:18s}: σ/m = {sigma_m:.3f} cm²/g, "
              f"r_core = {r1*1000:.0f} pc, χ²/dof = {chi2:.1f}/{ndof}")
        print(f"    σ_los at data R: {np.array2string(slos_at_obs, precision=1, separator=', ')} km/s")

    print()

    # ── Density profile plot ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    rho_nfw_plot = np.array([rho_nfw(r, rho_s, r_s) for r in r_arr])
    ax.loglog(r_arr * 1000, rho_nfw_plot, 'k:', lw=1.5, label='NFW', alpha=0.6)
    for bp, res in zip(bps, results):
        r1 = res['r_core_pc'] / 1000.0
        rho0 = res['rho0']
        rho_plot = np.array([rho_sidm(r, rho_s, r_s, r1, rho0) for r in r_arr])
        ax.loglog(r_arr * 1000, rho_plot, styles[bp['label']],
                  color=colors[bp['label']], lw=2,
                  label=f"{bp['label']} ($r_{{\\rm core}}$ = {res['r_core_pc']:.0f} pc)")
        ax.axvline(res['r_core_pc'], color=colors[bp['label']],
                   ls=':', lw=0.8, alpha=0.5)
    ax.axvline(r_half * 1000, color='gray', ls='--', lw=1, alpha=0.6)
    ax.text(r_half * 1050, ax.get_ylim()[0] * 3, '$r_{1/2}$', fontsize=9, color='gray')
    ax.set_xlabel('$r$ [pc]', fontsize=12)
    ax.set_ylabel(r'$\rho_{\rm DM}$ [M$_\odot$ kpc$^{-3}$]', fontsize=12)
    ax.set_title('Dark Matter Density Profile', fontsize=12)
    ax.legend(fontsize=9, loc='lower left')
    ax.set_xlim(1, 1e4)

    # ── σ_los profile plot ──
    ax = axes[1]
    ax.errorbar(R_obs * 1000, s_obs, yerr=s_err,
                fmt='ko', ms=6, capsize=3, label='Walker+2009', zorder=10)
    ax.plot(R_proj * 1000, slos_nfw, 'k:', lw=1.5, label='NFW (no SIDM)', alpha=0.6)
    for bp, res in zip(bps, results):
        ax.plot(R_proj * 1000, res['slos'], styles[bp['label']],
                color=colors[bp['label']], lw=2,
                label=f"{bp['label']} ($\\chi^2$ = {res['chi2']:.1f})")
    ax.axvline(r_half * 1000, color='gray', ls='--', lw=1, alpha=0.6)
    ax.text(r_half * 1050, 5, '$r_{1/2}$', fontsize=9, color='gray')
    ax.set_xlabel('Projected radius $R$ [pc]', fontsize=12)
    ax.set_ylabel(r'$\sigma_{\rm los}$ [km/s]', fontsize=12)
    ax.set_title(r'Fornax Stellar Velocity Dispersion (Jeans, $\beta = 0$)', fontsize=12)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xlim(0, 2500)
    ax.set_ylim(0, 18)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'fornax_jeans_prediction.png')
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # ── Summary table ──
    print()
    print("  ┌────────┬────────────┬───────────┬────────────┬────────────┐")
    print("  │ Model  │ σ/m cm²/g  │ r_core pc │ χ²/dof     │ Verdict    │")
    print("  ├────────┼────────────┼───────────┼────────────┼────────────┤")
    print(f"  │ {'NFW':6s} │ {'—':10s} │ {'0':9s} │ {chi2_nfw:6.1f}/{ndof:<3d}  │"
          f" {'CUSP':10s} │")
    for bp, res in zip(bps, results):
        verdict = "✓" if res['chi2'] < chi2_nfw else "✗ worse"
        print(f"  │ {bp['label']:6s} │ {res['sigma_m']:10.3f} │ "
              f"{res['r_core_pc']:7.0f}   │ {res['chi2']:6.1f}/{ndof:<3d}  │"
              f" {verdict:10s} │")
    print("  └────────┴────────────┴───────────┴────────────┴────────────┘")
    print()

    # ── Key physical insight ──
    best = min(results, key=lambda x: x['chi2'])
    print(f"  Best fit: {best['label']} (r_core = {best['r_core_pc']:.0f} pc, "
          f"χ² = {best['chi2']:.1f})")
    print()
    print("  Physical insight:")
    print("    SIDM core → flatter σ_los(R) profile inside r_core")
    print("    NFW cusp → σ_los rises towards center (over-concentrated)")
    print("    Data shows flat ~11.5 km/s profile → consistent with cored DM")
    print()


if __name__ == '__main__':
    main()
