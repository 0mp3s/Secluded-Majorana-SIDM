#!/usr/bin/env python3
"""
model_validations/rar_mcgaugh/rar_analysis.py
==============================================
§7.3 — Radial Acceleration Relation (RAR) from SIDM rotation curves.

The RAR (McGaugh+2016) relates observed gravitational acceleration g_obs
to baryonic acceleration g_bar via:
    g_obs = g_bar / (1 - exp(-sqrt(g_bar / g†)))
with g† = 1.2e-10 m/s².

This script:
  1. Loads SPARC rotation curve data (V_obs, V_bar at each radius)
  2. Computes g_obs and g_bar from the observed and baryonic velocities
  3. Fits the SIDM rotation curves (using fit_sparc_baryons.py logic)
  4. Computes SIDM-predicted g_DM at each radius
  5. Overlays: (a) data, (b) McGaugh empirical curve, (c) SIDM prediction
  6. Computes scatter for dwarfs vs spirals (key prediction)

Output: model_validations/rar_mcgaugh/output/
"""
import sys, os, csv, math
import numpy as np

# ── path bootstrap ──
_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
RC_DIR = os.path.join(_ROOT, 'predictions', 'rotation_curves')
sys.path.insert(0, RC_DIR)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v22_raw_scan import sigma_T_vpm
from fit_sparc_baryons import (
    load_rotation_data, nfw_params_with_cM, fit_galaxy_ac_sidm,
    G_N, KPC_CM, MSUN_G, SEC_PER_GYR, KM_S_TO_CM_S
)
from config_loader import load_config

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ── constants ──
G_SI = 6.674e-11           # m³/(kg·s²)
MSUN_KG = 1.989e30
KPC_M = 3.086e19
G_DAGGER = 1.2e-10         # m/s² — McGaugh+2016 critical acceleration

# BP1 parameters from config
_CFG = load_config(__file__)
def _get_bp(label):
    for bp in _CFG.get("benchmark_points", []):
        if bp["label"] == label:
            return bp
    return {}
_BP1 = _get_bp("BP1")
M_CHI = _BP1.get("m_chi_GeV", 20.69)
M_PHI = _BP1.get("m_phi_MeV", 11.34) * 1e-3   # → GeV
ALPHA = _BP1.get("alpha", 1.048e-3)
HALO_AGE_GYR = _CFG.get("halo_age_Gyr", 10.0)

# ── output directory ──
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── galaxy metadata ──
GALAXY_META = {
    'DDO_154':  {'V_max': 47,  'category': 'dwarf'},
    'IC_2574':  {'V_max': 66,  'category': 'dwarf'},
    'NGC_2366': {'V_max': 55,  'category': 'dwarf'},
    'NGC_2403': {'V_max': 136, 'category': 'spiral'},
    'NGC_2976': {'V_max': 90,  'category': 'spiral'},
    'NGC_3198': {'V_max': 150, 'category': 'spiral'},
    'UGC_128':  {'V_max': 58,  'category': 'dwarf'},
}


def mcgaugh_rar(g_bar):
    """McGaugh+2016 empirical RAR: g_obs = g_bar / (1 - exp(-sqrt(g_bar/g†)))."""
    x = np.sqrt(np.maximum(g_bar / G_DAGGER, 1e-30))
    return g_bar / (1.0 - np.exp(-x))


def v_to_g(v_km_s, r_kpc):
    """Convert circular velocity to centripetal acceleration in m/s²."""
    v_m_s = v_km_s * 1e3
    r_m = r_kpc * KPC_M
    return v_m_s**2 / r_m


def main():
    print("=" * 80)
    print("§7.3 — RADIAL ACCELERATION RELATION (RAR)")
    print("=" * 80)
    print(f"  BP1: m_χ = {M_CHI} GeV, m_φ = {M_PHI*1e3:.2f} MeV, α = {ALPHA:.3e}")
    print(f"  g† = {G_DAGGER:.1e} m/s² (McGaugh+2016)")
    print()

    # ── load SPARC data ──
    data_csv = os.path.join(RC_DIR, 'sparc_rotation_data.csv')
    galaxies = load_rotation_data(data_csv)
    print(f"  Loaded {len(galaxies)} galaxies from SPARC data")

    # ── fit each galaxy with SIDM and extract RAR data ──
    all_g_bar = []
    all_g_obs = []
    all_g_sidm = []
    all_cat = []

    results_by_cat = {'dwarf': {'g_bar': [], 'g_obs': [], 'g_sidm': []},
                      'spiral': {'g_bar': [], 'g_obs': [], 'g_sidm': []}}

    for name in sorted(galaxies.keys()):
        if name not in GALAXY_META:
            print(f"  WARNING: {name} not in metadata, skipping")
            continue

        meta = GALAXY_META[name]
        cat = meta['category']
        V_max = meta['V_max']
        g = galaxies[name]

        print(f"\n  {name} ({cat}, V_max={V_max} km/s):")

        # NFW parameters from V_max with c-M relation
        rho_s, r_s, c, M_200 = nfw_params_with_cM(V_max)

        # σ/m at representative velocity (V_max/√2 ≈ characteristic velocity)
        v_char = V_max / math.sqrt(2)
        sigma_m = sigma_T_vpm(M_CHI, M_PHI, ALPHA, v_char)
        sigma_v = V_max / math.sqrt(3)  # 1D velocity dispersion

        t_age_s = HALO_AGE_GYR * SEC_PER_GYR

        # Fit SIDM rotation curve
        upsilon, chi2, ndof, r_1, V_DM = fit_galaxy_ac_sidm(
            g['r'], g['V_obs'], g['V_err'], g['V_bar'],
            rho_s, r_s, sigma_m, sigma_v, t_age_s
        )

        chi2_dof = chi2 / max(ndof, 1)
        print(f"    Υ_* = {upsilon:.3f}, χ²/dof = {chi2_dof:.2f}, r_1 = {r_1:.3f} kpc")

        # Compute g_bar and g_obs at each data point
        for i in range(len(g['r'])):
            r = g['r'][i]
            if r < 0.1:  # skip very inner points where g_bar → ∞
                continue

            g_obs_i = v_to_g(g['V_obs'][i], r)
            g_bar_i = v_to_g(g['V_bar'][i] * math.sqrt(upsilon), r)

            # SIDM prediction: g_obs = g_bar + g_DM
            g_dm_i = v_to_g(V_DM[i], r)
            g_sidm_i = g_bar_i + g_dm_i

            if g_bar_i > 0 and g_obs_i > 0:
                all_g_bar.append(g_bar_i)
                all_g_obs.append(g_obs_i)
                all_g_sidm.append(g_sidm_i)
                all_cat.append(cat)

                results_by_cat[cat]['g_bar'].append(g_bar_i)
                results_by_cat[cat]['g_obs'].append(g_obs_i)
                results_by_cat[cat]['g_sidm'].append(g_sidm_i)

    all_g_bar = np.array(all_g_bar)
    all_g_obs = np.array(all_g_obs)
    all_g_sidm = np.array(all_g_sidm)
    all_cat = np.array(all_cat)

    # ── RAR statistics ──
    print("\n" + "=" * 80)
    print("RAR STATISTICS:")
    print("=" * 80)

    # Residuals relative to McGaugh curve
    g_mcgaugh = mcgaugh_rar(all_g_bar)
    resid_obs = np.log10(all_g_obs / g_mcgaugh)
    resid_sidm = np.log10(all_g_sidm / g_mcgaugh)

    print(f"  All galaxies:")
    print(f"    N data points: {len(all_g_bar)}")
    print(f"    Observed scatter (dex): {np.std(resid_obs):.4f}")
    print(f"    SIDM scatter (dex):     {np.std(resid_sidm):.4f}")

    for cat in ['dwarf', 'spiral']:
        mask = all_cat == cat
        if np.sum(mask) == 0:
            continue
        r_obs = np.log10(all_g_obs[mask] / mcgaugh_rar(all_g_bar[mask]))
        r_sidm = np.log10(all_g_sidm[mask] / mcgaugh_rar(all_g_bar[mask]))
        print(f"\n  {cat}s ({np.sum(mask)} points):")
        print(f"    Observed scatter (dex): {np.std(r_obs):.4f}")
        print(f"    SIDM scatter (dex):     {np.std(r_sidm):.4f}")

    # ── KEY PREDICTION ──
    dwarf_mask = all_cat == 'dwarf'
    spiral_mask = all_cat == 'spiral'
    if np.sum(dwarf_mask) > 0 and np.sum(spiral_mask) > 0:
        scatter_dwarf = np.std(np.log10(all_g_sidm[dwarf_mask] / mcgaugh_rar(all_g_bar[dwarf_mask])))
        scatter_spiral = np.std(np.log10(all_g_sidm[spiral_mask] / mcgaugh_rar(all_g_bar[spiral_mask])))
        print(f"\n  ★ KEY PREDICTION:")
        print(f"    SIDM scatter dwarfs:  {scatter_dwarf:.4f} dex")
        print(f"    SIDM scatter spirals: {scatter_spiral:.4f} dex")
        print(f"    Ratio: {scatter_spiral/scatter_dwarf:.2f}" if scatter_dwarf > 0
              else "    Cannot compute ratio")
        print(f"    Prediction: dwarfs should have LESS scatter (more thermalized)")

    # ── Plot ──
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))

    # McGaugh curve
    g_bar_curve = np.logspace(-13, -8, 500)
    g_obs_curve = mcgaugh_rar(g_bar_curve)
    ax.plot(g_bar_curve, g_obs_curve, 'k-', lw=2, label='McGaugh+2016 (empirical)')
    ax.plot(g_bar_curve, g_bar_curve, 'k--', lw=0.5, alpha=0.3, label='g_obs = g_bar (no DM)')

    # Data points by category
    colors = {'dwarf': '#2196F3', 'spiral': '#FF5722'}
    markers = {'dwarf': 'o', 'spiral': 's'}

    for cat in ['dwarf', 'spiral']:
        mask = all_cat == cat
        if np.sum(mask) == 0:
            continue
        ax.scatter(all_g_bar[mask], all_g_obs[mask],
                   c=colors[cat], marker=markers[cat], s=25, alpha=0.6,
                   edgecolors='none', label=f'{cat}s (observed)')
        ax.scatter(all_g_bar[mask], all_g_sidm[mask],
                   c=colors[cat], marker=markers[cat], s=25, alpha=0.6,
                   edgecolors='k', linewidths=0.5,
                   label=f'{cat}s (SIDM prediction)')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$g_{\rm bar}$ [m/s$^2$]', fontsize=13)
    ax.set_ylabel(r'$g_{\rm obs}$ [m/s$^2$]', fontsize=13)
    ax.set_title('Radial Acceleration Relation — SIDM vs McGaugh+2016', fontsize=14)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_xlim(1e-13, 1e-8)
    ax.set_ylim(1e-12, 1e-8)
    ax.grid(True, alpha=0.2)

    plot_path = os.path.join(OUT_DIR, 'rar_comparison.png')
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\n  Plot saved: {plot_path}")

    # ── Save CSV ──
    csv_path = os.path.join(OUT_DIR, 'rar_data.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['galaxy', 'category', 'g_bar_m_s2', 'g_obs_m_s2',
                     'g_sidm_m_s2', 'g_mcgaugh_m_s2'])
        for i in range(len(all_g_bar)):
            gal_idx = i  # simplified — we track by index
            w.writerow([
                '',  # galaxy name not tracked per-point in current structure
                all_cat[i],
                f"{all_g_bar[i]:.6e}",
                f"{all_g_obs[i]:.6e}",
                f"{all_g_sidm[i]:.6e}",
                f"{mcgaugh_rar(all_g_bar[i]):.6e}",
            ])
    print(f"  CSV saved: {csv_path}")

    print("\n✓ §7.3 RAR analysis DONE.")


if __name__ == '__main__':
    main()
