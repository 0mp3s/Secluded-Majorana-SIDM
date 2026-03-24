#!/usr/bin/env python3
"""
predictions/mc_diversity/mc_diversity.py
========================================
Monte Carlo halo diversity test: can SIDM + cosmological scatter
in halo concentration reproduce the observed scatter in V(2 kpc)
and central density?

Method
------
1. Draw N halos with V_max uniform-in-log over [30, 250] km/s.
2. For each halo, compute the median concentration c(M) from
   Dutton & Macciò (2014), then scatter c by σ_{log c} = 0.16 dex.
3. Build an NFW profile from (V_max, c) and compute the SIDM
   thermalized core via the Kaplinghat+2016 criterion:
       ρ(r_1) × (σ/m) × v × t_age = 1
4. Record V_SIDM(2 kpc) and ρ_core = ρ_NFW(r_1).
5. Plot the scatter cloud vs SPARC observed points.
6. Compare BP1 (λ = 1.91, near Ramsauer-Townsend) with MAP
   (λ = 48.6, deep resonant): BP1 should produce wider scatter
   because |d ln σ/d ln v| is steep near the resonance.

Output
------
- Figure 1: V(2 kpc) vs V_max (Oman+2015 diversity plot)
- Figure 2: ρ_core vs V_max
- CSV with all MC results
"""
import sys, os, csv, math, time
import numpy as np

# ---------- path bootstrap ----------
_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
# ------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config_loader import load_config
from v22_raw_scan import sigma_T_vpm

# JIT warmup
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ---- constants ----
G_N          = 4.302e-6       # kpc (km/s)^2 / M_sun
KPC_CM       = 3.086e21       # cm
MSUN_G       = 1.989e33       # g
SEC_PER_GYR  = 3.156e16       # s
KM_S_TO_CM_S = 1e5

# =========================================================================
#  NFW helpers
# =========================================================================

def nfw_rho(r_kpc, rho_s, r_s):
    x = max(r_kpc / r_s, 1e-8)
    return rho_s / (x * (1.0 + x)**2)


def nfw_mass(r_kpc, rho_s, r_s):
    x = r_kpc / r_s
    return 4.0 * math.pi * rho_s * r_s**3 * (math.log(1 + x) - x / (1 + x))


def nfw_v_circ(r_kpc, rho_s, r_s):
    if r_kpc < 1e-8:
        return 0.0
    return math.sqrt(G_N * nfw_mass(r_kpc, rho_s, r_s) / r_kpc)


def nfw_params_from_vmax(V_max, c=12.0):
    rho_crit = 126.0
    f_c = math.log(1.0 + c) - c / (1.0 + c)
    delta_c = (200.0 / 3.0) * c**3 / f_c
    rho_s = rho_crit * delta_c

    x_max = 2.163
    f_xmax = math.log(1 + x_max) - x_max / (1 + x_max)
    r_s = math.sqrt(V_max**2 * x_max /
                    (G_N * 4.0 * math.pi * rho_s * f_xmax))
    return rho_s, r_s


def concentration_mass(M_200_Msun, h=0.674):
    """Dutton & Macciò 2014 eq.7 (Planck cosmology, z=0)."""
    log_M = math.log10(max(M_200_Msun * h / 1e12, 1e-30))
    return 10**(0.905 - 0.101 * log_M)


def nfw_params_self_consistent(V_max):
    """NFW params with iterated c(M) relation.
    Returns (rho_s, r_s, c, M_200).
    """
    rho_crit = 126.0
    c = 12.0
    for _ in range(30):
        rho_s, r_s = nfw_params_from_vmax(V_max, c=c)
        R_200 = c * r_s
        M_200 = (4.0 / 3.0) * math.pi * 200.0 * rho_crit * R_200**3
        c_new = concentration_mass(M_200)
        if abs(c_new - c) / max(c, 0.1) < 1e-3:
            break
        c = 0.5 * (c + c_new)
    rho_s, r_s = nfw_params_from_vmax(V_max, c=c)
    return rho_s, r_s, c, M_200


# =========================================================================
#  SIDM core
# =========================================================================

def find_r1(rho_s, r_s, sigma_over_m, sigma_v_km_s, t_age_s):
    v_cm = sigma_v_km_s * KM_S_TO_CM_S
    target_rho_cgs = 1.0 / (sigma_over_m * v_cm * t_age_s)
    target_rho = target_rho_cgs * KPC_CM**3 / MSUN_G

    r_lo, r_hi = 1e-4, 10.0 * r_s
    if nfw_rho(r_lo, rho_s, r_s) < target_rho:
        return r_lo
    if nfw_rho(r_hi, rho_s, r_s) > target_rho:
        return r_hi
    for _ in range(100):
        mid = 0.5 * (r_lo + r_hi)
        if nfw_rho(mid, rho_s, r_s) > target_rho:
            r_lo = mid
        else:
            r_hi = mid
    return 0.5 * (r_lo + r_hi)


def sidm_v_circ(r_kpc, rho_s, r_s, r_1, rho_core):
    if r_kpc < 1e-8:
        return 0.0
    if r_kpc <= r_1:
        M_enc = (4.0 / 3.0) * math.pi * rho_core * r_kpc**3
    else:
        M_core = (4.0 / 3.0) * math.pi * rho_core * r_1**3
        M_enc = M_core + nfw_mass(r_kpc, rho_s, r_s) - nfw_mass(r_1, rho_s, r_s)
    return math.sqrt(G_N * M_enc / r_kpc)


# =========================================================================
#  SPARC loader
# =========================================================================

def load_sparc():
    csv_path = os.path.join(_ROOT, 'predictions', 'rotation_curves',
                            'sparc_subset.csv')
    gals = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            sigma_v = float(row['sigma_v_km_s'])
            r_core = float(row['r_core_kpc'])
            # Estimate ρ_central from isothermal core:
            #   ρ_core ≈ 9σ²/(4πG r²_core)  [M_sun/kpc³]
            rho_est = 9.0 * sigma_v**2 / (4.0 * math.pi * G_N * r_core**2)
            gals.append({
                'name':       row['name'],
                'V_max':      float(row['V_max_km_s']),
                'V_2kpc_obs': float(row['V_2kpc_km_s']),
                'V_2kpc_err': float(row['V_2kpc_err']),
                'r_core_obs': r_core,
                'sigma_v':    sigma_v,
                'rho_central': rho_est,
                'category':   row['category'],
            })
    return gals


# =========================================================================
#  Monte Carlo
# =========================================================================

def run():
    t0 = time.time()
    cfg = load_config(__file__)

    bps          = cfg['benchmark_points']
    n_halos      = cfg.get('n_halos', 1000)
    vmax_lo, vmax_hi = cfg.get('V_max_range_km_s', [30, 250])
    sigma_logc   = cfg.get('sigma_log10_c', 0.16)
    t_age_gyr    = cfg.get('halo_age_Gyr', 10.0)
    seed         = cfg.get('seed', 42)
    out_dir      = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    t_age_s = t_age_gyr * SEC_PER_GYR
    rng = np.random.default_rng(seed)

    # --- Draw V_max uniform-in-log ---
    log_vmax = rng.uniform(np.log10(vmax_lo), np.log10(vmax_hi), n_halos)
    V_max_arr = 10**log_vmax

    # --- Compute median c(M) for each V_max, then scatter ---
    c_median = np.empty(n_halos)
    M_200_arr = np.empty(n_halos)
    for i, vm in enumerate(V_max_arr):
        _, _, c_med, M200 = nfw_params_self_consistent(vm)
        c_median[i] = c_med
        M_200_arr[i] = M200

    # log-normal scatter: log10(c) ~ N(log10(c_median), sigma_logc)
    log_c_scattered = np.log10(c_median) + rng.normal(0, sigma_logc, n_halos)
    c_scattered = 10**log_c_scattered

    print("=" * 76)
    print("  Monte Carlo Halo Diversity Simulation")
    print(f"  {n_halos} halos, V_max ∈ [{vmax_lo}, {vmax_hi}] km/s")
    print(f"  σ(log₁₀ c) = {sigma_logc:.2f} dex, t_age = {t_age_gyr} Gyr")
    print(f"  seed = {seed}")
    print("=" * 76)

    # --- Load SPARC for overlay ---
    sparc = load_sparc()

    # --- Run SIDM for each BP ---
    all_csv = []
    mc_results = {}  # mc_results[label] = dict of arrays

    for bp in bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi_GeV = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        lam = alpha * m_chi / m_phi_GeV

        V2_sidm = np.empty(n_halos)
        V2_nfw  = np.empty(n_halos)
        rho_core_arr = np.empty(n_halos)
        r1_arr = np.empty(n_halos)
        sigma_m_arr = np.empty(n_halos)

        print(f"\n  {label}: λ = {lam:.2f}, processing {n_halos} halos ... ",
              end='', flush=True)
        t_bp = time.time()

        for i in range(n_halos):
            c = c_scattered[i]
            vm = V_max_arr[i]

            rho_s, r_s = nfw_params_from_vmax(vm, c=c)

            # Inner velocity dispersion: use V_circ at characteristic
            # inner radius (1 kpc or r_s, whichever is smaller).
            # This makes σ_v c-dependent: higher c → more mass at 1 kpc
            # → deeper potential → higher σ_v → different σ/m.
            r_inner = min(1.0, r_s)
            v_circ_inner = nfw_v_circ(r_inner, rho_s, r_s)
            sigma_v = max(v_circ_inner / math.sqrt(2.0), 3.0)  # floor 3 km/s
            v_rel = sigma_v * math.sqrt(2.0)

            sigma_m = sigma_T_vpm(m_chi, m_phi_GeV, alpha, v_rel)
            r_1 = find_r1(rho_s, r_s, sigma_m, sigma_v, t_age_s)
            rho_core = nfw_rho(r_1, rho_s, r_s)

            V2_nfw[i]  = nfw_v_circ(2.0, rho_s, r_s)
            V2_sidm[i] = sidm_v_circ(2.0, rho_s, r_s, r_1, rho_core)
            rho_core_arr[i] = rho_core
            r1_arr[i] = r_1
            sigma_m_arr[i] = sigma_m

            all_csv.append({
                'bp': label, 'V_max': vm, 'c_median': c_median[i],
                'c_draw': c, 'M_200': M_200_arr[i],
                'sigma_m': sigma_m, 'v_rel': v_rel,
                'r_1_kpc': r_1, 'rho_core_Msun_kpc3': rho_core,
                'V2_nfw': V2_nfw[i], 'V2_sidm': V2_sidm[i],
            })

        dt = time.time() - t_bp
        print(f"{dt:.1f} s")

        mc_results[label] = {
            'V_max': V_max_arr, 'c': c_scattered,
            'V2_sidm': V2_sidm, 'V2_nfw': V2_nfw,
            'rho_core': rho_core_arr, 'r_1': r1_arr,
            'sigma_m': sigma_m_arr,
        }

        # --- Summary stats in V_max bins ---
        edges = [30, 60, 100, 160, 250]
        print(f"\n  {'V_max bin':>15s} {'<V2_nfw>':>9s} {'<V2_sidm>':>10s} "
              f"{'σ(V2_sidm)':>11s} {'<ρ_core>':>12s} {'σ(log ρ)':>10s}")
        print("  " + "-" * 68)
        for j in range(len(edges) - 1):
            mask = (V_max_arr >= edges[j]) & (V_max_arr < edges[j+1])
            if mask.sum() == 0:
                continue
            v2n_m = np.mean(V2_nfw[mask])
            v2s_m = np.mean(V2_sidm[mask])
            v2s_s = np.std(V2_sidm[mask])
            rho_m = np.mean(rho_core_arr[mask])
            rho_ls = np.std(np.log10(np.maximum(rho_core_arr[mask], 1e-3)))
            print(f"  [{edges[j]:>3d},{edges[j+1]:>3d}) km/s "
                  f"{v2n_m:>9.1f} {v2s_m:>10.1f} {v2s_s:>11.1f} "
                  f"{rho_m:>12.1f} {rho_ls:>10.3f}")

    # =====================================================================
    #  Write CSV
    # =====================================================================
    csv_path = os.path.join(out_dir, "mc_diversity.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "bp", "V_max", "c_median", "c_draw", "M_200",
            "sigma_m", "v_rel", "r_1_kpc", "rho_core_Msun_kpc3",
            "V2_nfw", "V2_sidm",
        ])
        w.writeheader()
        w.writerows(all_csv)
    print(f"\n  CSV: {csv_path}")

    # =====================================================================
    #  Figure 1: V(2kpc) vs V_max — the Oman+2015 diversity plot
    # =====================================================================
    n_bp = len(bps)
    fig, axes = plt.subplots(1, n_bp, figsize=(6.5 * n_bp, 6), sharey=True,
                             squeeze=False)

    # NFW reference curve (no scatter, no SIDM)
    vm_line = np.linspace(30, 250, 100)
    v2_nfw_ref = []
    for v in vm_line:
        _, _, c_ref, _ = nfw_params_self_consistent(v)
        rs, rss = nfw_params_from_vmax(v, c=c_ref)
        v2_nfw_ref.append(nfw_v_circ(2.0, rs, rss))

    for col, bp in enumerate(bps):
        ax = axes[0][col]
        label = bp['label']
        lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000)
        mc = mc_results[label]

        # NFW reference
        ax.plot(vm_line, v2_nfw_ref, 'k--', lw=1.5, alpha=0.4,
                label='NFW (no SIDM)')

        # MC cloud
        ax.scatter(mc['V_max'], mc['V2_sidm'], s=6, alpha=0.25,
                   c='steelblue', edgecolors='none',
                   label=f'MC halos (N={n_halos})')

        # SPARC data
        sp_vm = [g['V_max'] for g in sparc]
        sp_v2 = [g['V_2kpc_obs'] for g in sparc]
        sp_e  = [g['V_2kpc_err'] for g in sparc]
        ax.errorbar(sp_vm, sp_v2, yerr=sp_e, fmt='ro', ms=6,
                    capsize=3, zorder=5, label='SPARC observed')

        # 1:1 line
        ax.plot([30, 250], [30, 250], 'gray', ls=':', alpha=0.3)

        ax.set_xlabel(r'$V_{\rm max}$ [km s$^{-1}$]', fontsize=12)
        if col == 0:
            ax.set_ylabel(r'$V(2\,{\rm kpc})$ [km s$^{-1}$]', fontsize=12)
        ax.set_title(f'{label}  ($\\lambda = {lam:.1f}$)', fontsize=13)
        ax.set_xlim(25, 260)
        ax.set_ylim(0, 180)
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(True, alpha=0.3)

    fig.suptitle('Rotation-Curve Diversity: SIDM + Concentration Scatter vs SPARC',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    p1 = os.path.join(out_dir, 'mc_v2kpc_vs_vmax.png')
    fig.savefig(p1, dpi=150, bbox_inches='tight')
    fig.savefig(p1.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Fig1: {p1}")

    # =====================================================================
    #  Figure 2: ρ_core vs V_max
    # =====================================================================
    # Convert ρ_core from M_sun/kpc³ to M_sun/pc³ for readability
    fig, axes = plt.subplots(1, n_bp, figsize=(6.5 * n_bp, 6), sharey=True,
                             squeeze=False)

    for col, bp in enumerate(bps):
        ax = axes[0][col]
        label = bp['label']
        lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000)
        mc = mc_results[label]

        rho_pc3 = mc['rho_core'] / 1e9  # M_sun/kpc³ → M_sun/pc³

        ax.scatter(mc['V_max'], rho_pc3, s=6, alpha=0.25,
                   c='steelblue', edgecolors='none',
                   label=f'MC halos (N={n_halos})')

        # SPARC observed ρ_central (estimated from σ_v, r_core)
        sp_vm = [g['V_max'] for g in sparc]
        sp_rho = [g['rho_central'] / 1e9 for g in sparc]  # → M_sun/pc³
        ax.scatter(sp_vm, sp_rho, c='red', s=50, marker='D',
                   edgecolors='k', linewidths=0.5, zorder=5,
                   label='SPARC (estimated)')
        for g in sparc:
            ax.annotate(g['name'], (g['V_max'], g['rho_central'] / 1e9),
                        fontsize=5, alpha=0.7,
                        xytext=(3, 3), textcoords='offset points')

        ax.set_xlabel(r'$V_{\rm max}$ [km s$^{-1}$]', fontsize=12)
        if col == 0:
            ax.set_ylabel(r'$\rho_{\rm core}$ [M$_\odot$ pc$^{-3}$]',
                          fontsize=12)
        ax.set_title(f'{label}  ($\\lambda = {lam:.1f}$)', fontsize=13)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(25, 260)
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(True, alpha=0.3, which='both')

    fig.suptitle('Central Density Diversity: SIDM + Concentration Scatter',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    p2 = os.path.join(out_dir, 'mc_rho_core_vs_vmax.png')
    fig.savefig(p2, dpi=150, bbox_inches='tight')
    fig.savefig(p2.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Fig2: {p2}")

    # =====================================================================
    #  Figure 3: overlay BP1 vs MAP on single panel
    # =====================================================================
    if len(bps) >= 2:
        fig, ax = plt.subplots(figsize=(9, 7))

        # NFW ref
        ax.plot(vm_line, v2_nfw_ref, 'k--', lw=2, alpha=0.4,
                label='NFW (no SIDM)')

        bp_styles = [
            ('steelblue', 0.15, 'o'),
            ('darkorange', 0.15, 's'),
            ('forestgreen', 0.15, '^'),
        ]
        for idx, bp in enumerate(bps):
            label_bp = bp['label']
            lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000)
            mc = mc_results[label_bp]
            c, a, m = bp_styles[idx % len(bp_styles)]
            ax.scatter(mc['V_max'], mc['V2_sidm'], s=8, alpha=a,
                       c=c, edgecolors='none', marker=m,
                       label=f'{label_bp} ($\\lambda = {lam:.1f}$)')

        # SPARC
        sp_vm = [g['V_max'] for g in sparc]
        sp_v2 = [g['V_2kpc_obs'] for g in sparc]
        sp_e  = [g['V_2kpc_err'] for g in sparc]
        ax.errorbar(sp_vm, sp_v2, yerr=sp_e, fmt='ko', ms=7,
                    capsize=3, zorder=5, label='SPARC observed')
        for g in sparc:
            ax.annotate(g['name'], (g['V_max'], g['V_2kpc_obs']),
                        fontsize=6, alpha=0.7,
                        xytext=(4, 4), textcoords='offset points')

        ax.plot([30, 250], [30, 250], 'gray', ls=':', alpha=0.3)
        ax.set_xlabel(r'$V_{\rm max}$ [km s$^{-1}$]', fontsize=13)
        ax.set_ylabel(r'$V(2\,{\rm kpc})$ [km s$^{-1}$]', fontsize=13)
        ax.set_title('Halo Diversity: Full MC Comparison', fontsize=14)
        ax.set_xlim(25, 260)
        ax.set_ylim(0, 180)
        ax.legend(fontsize=9, loc='upper left')
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        p3 = os.path.join(out_dir, 'mc_diversity_overlay.png')
        fig.savefig(p3, dpi=150, bbox_inches='tight')
        fig.savefig(p3.replace('.png', '.pdf'), bbox_inches='tight')
        plt.close(fig)
        print(f"  Fig3: {p3}")

    # =====================================================================
    #  Quantitative diversity comparison
    # =====================================================================
    print(f"\n{'='*76}")
    print("  DIVERSITY QUANTIFICATION: σ(V2_sidm) / σ(V2_nfw) in V_max bins")
    print(f"{'='*76}")
    edges = [30, 60, 100, 160, 250]
    print(f"  {'V_max bin':>15s}", end='')
    for bp in bps:
        print(f"  {bp['label']:>12s}", end='')
    print(f"  {'NFW (no SC)':>12s}")
    print("  " + "-" * (15 + (n_bp + 1) * 14))

    for j in range(len(edges) - 1):
        mask = (V_max_arr >= edges[j]) & (V_max_arr < edges[j+1])
        if mask.sum() < 10:
            continue
        # NFW scatter from c-scatter alone (no SIDM)
        v2nfw = mc_results[bps[0]['label']]['V2_nfw'][mask]
        nfw_std = np.std(v2nfw)
        print(f"  [{edges[j]:>3d},{edges[j+1]:>3d}) km/s ", end='')
        for bp in bps:
            v2s = mc_results[bp['label']]['V2_sidm'][mask]
            print(f"  {np.std(v2s):>12.1f}", end='')
        print(f"  {nfw_std:>12.1f}")

    # SPARC coverage check
    print(f"\n  SPARC data points inside MC cloud:")
    for bp in bps:
        label = bp['label']
        mc = mc_results[label]
        inside = 0
        for g in sparc:
            # Check if any MC halo within ΔV_max=15 has V2 within 2σ
            near = np.abs(mc['V_max'] - g['V_max']) < 15
            if near.sum() == 0:
                continue
            v2_near = mc['V2_sidm'][near]
            lo, hi = np.percentile(v2_near, [5, 95])
            if lo <= g['V_2kpc_obs'] <= hi:
                inside += 1
        print(f"    {label}: {inside}/{len(sparc)} galaxies "
              f"within 5-95th percentile band")

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f} s")
    print("=" * 76)


if __name__ == "__main__":
    run()
