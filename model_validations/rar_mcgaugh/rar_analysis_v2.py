#!/usr/bin/env python3
"""
model_validations/rar_mcgaugh/rar_analysis_v2.py
=================================================
§7.3 — RAR with gas-fraction fix for dwarfs.

Fix: V_bar in SPARC CSV was pre-computed with Υ_*_default = 0.5:
    V_bar² = V_gas² + 0.5 × V_disk²

For gas-dominated dwarfs, fitting Υ_* on *all* of V_bar is wrong.
We split using literature gas fractions (f_gas = M_gas/M_bar):
    V_gas²    = f_gas × V_bar²
    V_disk²   = 2(1 - f_gas) × V_bar²    (undo the 0.5 factor)

Then fit:
    V_tot² = V_gas² + Υ_* × V_disk² + V_DM²
           = V_bar² × [f_gas + 2Υ_*(1-f_gas)] + V_DM²

Gas fractions from: Oh+2015, de Blok+2008, Lelli+2016.
"""
import sys, os, csv, math
import numpy as np

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
RC_DIR = os.path.join(_ROOT, 'predictions', 'rotation_curves')
sys.path.insert(0, RC_DIR)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v22_raw_scan import sigma_T_vpm
from fit_sparc_baryons import (
    load_rotation_data, nfw_params_with_cM,
    G_N, KPC_CM, MSUN_G, SEC_PER_GYR, KM_S_TO_CM_S,
    F_B, nfw_mass, solve_ac_ri,
)
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)  # JIT warmup

# ── constants ──
G_SI = 6.674e-11
KPC_M = 3.086e19
G_DAGGER = 1.2e-10     # m/s² — McGaugh+2016

# Benchmarks
BENCHMARKS = {
    'BP1': {'m_chi': 20.69, 'm_phi': 11.34e-3, 'alpha': 1.048e-3},
    'MAP': {'m_chi': 94.07, 'm_phi': 11.10e-3, 'alpha': 5.734e-3},
}
HALO_AGE_GYR = 10.0

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── galaxy metadata with gas fractions ──
# f_gas = M_gas / M_baryonic from literature.
# Sources: Oh+2015 (LITTLE THINGS), de Blok+2008, Adams+2014, Lelli+2016 (SPARC)
GALAXY_META = {
    'DDO_154':  {'V_max': 47,  'category': 'dwarf',  'f_gas': 0.95},  # Oh+2015
    'IC_2574':  {'V_max': 66,  'category': 'dwarf',  'f_gas': 0.82},  # Oh+2015
    'NGC_2366': {'V_max': 55,  'category': 'dwarf',  'f_gas': 0.88},  # Oh+2015
    'NGC_2403': {'V_max': 136, 'category': 'spiral', 'f_gas': 0.25},  # Lelli+2016
    'NGC_2976': {'V_max': 90,  'category': 'spiral', 'f_gas': 0.15},  # Adams+2014
    'NGC_3198': {'V_max': 150, 'category': 'spiral', 'f_gas': 0.30},  # Lelli+2016
    'UGC_128':  {'V_max': 58,  'category': 'dwarf',  'f_gas': 0.75},  # de Blok+2008
}

UPS_DEFAULT = 0.5  # Υ_* used to pre-compute V_bar in SPARC CSV


def mcgaugh_rar(g_bar):
    x = np.sqrt(np.maximum(g_bar / G_DAGGER, 1e-30))
    return g_bar / (1.0 - np.exp(-x))


def v_to_g(v_km_s, r_kpc):
    v_m_s = v_km_s * 1e3
    r_m = r_kpc * KPC_M
    return v_m_s**2 / r_m


def fit_galaxy_with_gas_fraction(r_data, V_obs, V_err, V_bar,
                                  rho_s, r_s, sigma_m, sigma_v, t_age_s,
                                  f_gas):
    """
    Like fit_galaxy_ac_sidm but properly handles gas fraction.

    V_bar_csv was computed as: V_bar² = V_gas² + 0.5 × V_disk²
    We decompose:
        V_gas   = sqrt(f_gas) × V_bar
        V_disk  = sqrt(2(1-f_gas)) × V_bar
    Then fit: V_tot² = V_gas² + Υ_* × V_disk² + V_DM²
    """
    V_bar_func = interp1d(r_data, V_bar, kind='linear',
                          bounds_error=False,
                          fill_value=(V_bar[0], V_bar[-1]))

    # Decompose V_bar into gas and stellar components
    V_gas2 = f_gas * V_bar**2
    V_disk2 = 2.0 * (1.0 - f_gas) * V_bar**2  # undo the 0.5 default

    upsilon = 0.5
    best_chi2 = 1e30
    V_DM_result = np.zeros(len(r_data))
    r_1_result = 0.0

    for iteration in range(8):
        r_min = max(r_data.min() * 0.3, 0.02)
        r_max = r_data.max() * 2.0
        r_fine = np.geomspace(r_min, r_max, 200)

        # AC'd DM mass with baryonic contribution using CURRENT Υ_*
        M_DM_AC = np.empty(len(r_fine))
        for i, rf in enumerate(r_fine):
            Vb = float(V_bar_func(rf))
            # Effective baryonic V² with gas fraction separation
            Vb_eff2 = f_gas * Vb**2 + upsilon * 2.0 * (1.0 - f_gas) * Vb**2
            m_bar = max(Vb_eff2 * rf / G_N, 0.0)
            ri = solve_ac_ri(rf, m_bar, rho_s, r_s)
            M_DM_AC[i] = (1 - F_B) * nfw_mass(ri, rho_s, r_s)

        # Density from mass profile
        rho_AC = np.gradient(M_DM_AC, r_fine) / (4.0 * np.pi * r_fine**2)
        rho_AC = np.maximum(rho_AC, 1e-10)

        # SIDM core radius
        v_cm = sigma_v * KM_S_TO_CM_S
        target_rho_cgs = 1.0 / (sigma_m * v_cm * t_age_s)
        target_rho = target_rho_cgs * KPC_CM**3 / MSUN_G

        r_1 = r_fine[-1]
        for i in range(len(r_fine)):
            if rho_AC[i] < target_rho:
                r_1 = r_fine[i]
                break
        rho_core = target_rho

        # V_DM at data radii with SIDM coring
        M_DM_AC_at_r1 = float(np.interp(r_1, r_fine, M_DM_AC))
        M_DM_AC_data = np.interp(r_data, r_fine, M_DM_AC)

        V_DM = np.empty(len(r_data))
        for i, r in enumerate(r_data):
            if r <= r_1:
                M_enc = (4.0 / 3.0) * math.pi * rho_core * r**3
            else:
                M_core = (4.0 / 3.0) * math.pi * rho_core * r_1**3
                M_enc = M_core + (M_DM_AC_data[i] - M_DM_AC_at_r1)
            V_DM[i] = math.sqrt(max(G_N * M_enc / r, 0))

        V_DM2 = V_DM**2

        # FIT: V_tot² = V_gas² + Υ_* × V_disk² + V_DM²
        def chi2_func(ups):
            V_tot2 = V_gas2 + ups * V_disk2 + V_DM2
            V_tot = np.sqrt(np.maximum(V_tot2, 1e-10))
            return np.sum(((V_obs - V_tot) / V_err)**2)

        result = minimize_scalar(chi2_func, bounds=(0.01, 3.0), method='bounded')
        new_upsilon = result.x
        new_chi2 = result.fun

        V_DM_result = V_DM.copy()
        r_1_result = r_1
        best_chi2 = new_chi2

        if abs(new_upsilon - upsilon) < 0.01:
            upsilon = new_upsilon
            break
        upsilon = new_upsilon

    ndof = len(V_obs) - 1
    return upsilon, best_chi2, ndof, r_1_result, V_DM_result


def main():
    print("=" * 80)
    print("§7.3 — RAR (v2: gas-fraction fix)")
    print("=" * 80)
    print(f"  g† = {G_DAGGER:.1e} m/s² (McGaugh+2016)")
    print(f"  Gas fractions: literature values per galaxy")
    print()

    data_csv = os.path.join(RC_DIR, 'sparc_rotation_data.csv')
    galaxies = load_rotation_data(data_csv)
    print(f"  Loaded {len(galaxies)} galaxies")

    for bp_name, bp in BENCHMARKS.items():
        print(f"\n{'='*80}")
        print(f"BENCHMARK: {bp_name}")
        print(f"  m_χ={bp['m_chi']} GeV, m_φ={bp['m_phi']*1e3:.2f} MeV, α={bp['alpha']:.3e}")
        print(f"{'='*80}")

        all_g_bar = []
        all_g_obs = []
        all_g_sidm = []
        all_cat = []
        fit_results = []

        for name in sorted(galaxies.keys()):
            if name not in GALAXY_META:
                continue

            meta = GALAXY_META[name]
            cat = meta['category']
            V_max = meta['V_max']
            f_gas = meta['f_gas']
            g = galaxies[name]

            rho_s, r_s, c, M_200 = nfw_params_with_cM(V_max)
            v_char = V_max / math.sqrt(2)
            sigma_m = sigma_T_vpm(bp['m_chi'], bp['m_phi'], bp['alpha'], v_char)
            sigma_v = V_max / math.sqrt(3)
            t_age_s = HALO_AGE_GYR * SEC_PER_GYR

            upsilon, chi2, ndof, r_1, V_DM = fit_galaxy_with_gas_fraction(
                g['r'], g['V_obs'], g['V_err'], g['V_bar'],
                rho_s, r_s, sigma_m, sigma_v, t_age_s,
                f_gas
            )

            chi2_dof = chi2 / max(ndof, 1)
            print(f"\n  {name} ({cat}, V_max={V_max}, f_gas={f_gas:.2f}):")
            print(f"    Υ_* = {upsilon:.3f}, χ²/dof = {chi2_dof:.2f}, r_1 = {r_1:.3f} kpc")

            fit_results.append({
                'name': name, 'category': cat, 'f_gas': f_gas,
                'upsilon': upsilon, 'chi2_dof': chi2_dof, 'r_1': r_1
            })

            # Compute RAR data points
            V_gas2 = f_gas * g['V_bar']**2
            V_disk2 = 2.0 * (1.0 - f_gas) * g['V_bar']**2

            for i in range(len(g['r'])):
                r = g['r'][i]
                if r < 0.1:
                    continue

                g_obs_i = v_to_g(g['V_obs'][i], r)
                # g_bar with proper Υ_* treatment
                V_bar_eff = math.sqrt(V_gas2[i] + upsilon * V_disk2[i])
                g_bar_i = v_to_g(V_bar_eff, r)

                g_dm_i = v_to_g(V_DM[i], r)
                g_sidm_i = g_bar_i + g_dm_i

                if g_bar_i > 0 and g_obs_i > 0:
                    all_g_bar.append(g_bar_i)
                    all_g_obs.append(g_obs_i)
                    all_g_sidm.append(g_sidm_i)
                    all_cat.append(cat)

        all_g_bar = np.array(all_g_bar)
        all_g_obs = np.array(all_g_obs)
        all_g_sidm = np.array(all_g_sidm)
        all_cat = np.array(all_cat)

        # RAR statistics
        print(f"\n{'='*80}")
        print(f"RAR STATISTICS ({bp_name}):")
        print(f"{'='*80}")

        g_mcgaugh = mcgaugh_rar(all_g_bar)
        resid_obs = np.log10(all_g_obs / g_mcgaugh)
        resid_sidm = np.log10(all_g_sidm / g_mcgaugh)

        print(f"  All galaxies ({len(all_g_bar)} points):")
        print(f"    Observed scatter: {np.std(resid_obs):.4f} dex")
        print(f"    SIDM scatter:     {np.std(resid_sidm):.4f} dex")

        for cat in ['dwarf', 'spiral']:
            mask = all_cat == cat
            if np.sum(mask) == 0:
                continue
            r_obs = np.log10(all_g_obs[mask] / mcgaugh_rar(all_g_bar[mask]))
            r_sidm = np.log10(all_g_sidm[mask] / mcgaugh_rar(all_g_bar[mask]))
            print(f"\n  {cat}s ({np.sum(mask)} points):")
            print(f"    Observed scatter: {np.std(r_obs):.4f} dex")
            print(f"    SIDM scatter:     {np.std(r_sidm):.4f} dex")

        # Summary table
        print(f"\n  {'Galaxy':>12s}  {'cat':>6s}  {'f_gas':>5s}  {'Υ_*':>6s}  {'χ²/dof':>8s}  {'r_1':>6s}")
        print(f"  {'-'*55}")
        for fr in fit_results:
            print(f"  {fr['name']:>12s}  {fr['category']:>6s}  {fr['f_gas']:>5.2f}  "
                  f"{fr['upsilon']:>6.3f}  {fr['chi2_dof']:>8.2f}  {fr['r_1']:>6.3f}")

        # ── Save CSV ──
        csv_path = os.path.join(OUT_DIR, f'rar_v2_{bp_name}.csv')
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['log10_g_bar', 'log10_g_obs', 'log10_g_sidm', 'category'])
            for j in range(len(all_g_bar)):
                w.writerow([f'{np.log10(all_g_bar[j]):.6f}',
                            f'{np.log10(all_g_obs[j]):.6f}',
                            f'{np.log10(all_g_sidm[j]):.6f}',
                            all_cat[j]])
        print(f"\n  CSV saved → {csv_path}")

        # ── RAR Plot ──
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Left: g_obs vs g_bar
        ax = axes[0]
        g_line = np.logspace(-13, -8, 200)
        ax.plot(g_line, mcgaugh_rar(g_line), 'k-', lw=2, label='McGaugh+2016')
        ax.plot(g_line, g_line, 'k--', lw=0.5, alpha=0.3)

        mask_d = all_cat == 'dwarf'
        mask_s = all_cat == 'spiral'
        if np.any(mask_d):
            ax.scatter(all_g_bar[mask_d], all_g_obs[mask_d],
                       c='royalblue', s=20, alpha=0.7, label='obs dwarf')
            ax.scatter(all_g_bar[mask_d], all_g_sidm[mask_d],
                       c='royalblue', s=20, alpha=0.7, marker='x', label='SIDM dwarf')
        if np.any(mask_s):
            ax.scatter(all_g_bar[mask_s], all_g_obs[mask_s],
                       c='firebrick', s=20, alpha=0.7, label='obs spiral')
            ax.scatter(all_g_bar[mask_s], all_g_sidm[mask_s],
                       c='firebrick', s=20, alpha=0.7, marker='x', label='SIDM spiral')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$g_{\rm bar}$ [m/s$^2$]')
        ax.set_ylabel(r'$g_{\rm obs}$ or $g_{\rm SIDM}$ [m/s$^2$]')
        ax.set_title(f'{bp_name} — RAR (v2 gas-frac fix)')
        ax.legend(fontsize=8)
        ax.set_xlim(1e-13, 1e-8)
        ax.set_ylim(1e-12, 1e-8)

        # Right: residuals
        ax = axes[1]
        g_mcg = mcgaugh_rar(all_g_bar)
        resid_obs_all = np.log10(all_g_obs / g_mcg)
        resid_sidm_all = np.log10(all_g_sidm / g_mcg)

        if np.any(mask_d):
            ax.scatter(np.log10(all_g_bar[mask_d]), resid_obs_all[mask_d],
                       c='royalblue', s=20, alpha=0.5, label=f'obs dwarf (σ={np.std(resid_obs_all[mask_d]):.3f})')
            ax.scatter(np.log10(all_g_bar[mask_d]), resid_sidm_all[mask_d],
                       c='royalblue', marker='x', s=25, alpha=0.7,
                       label=f'SIDM dwarf (σ={np.std(resid_sidm_all[mask_d]):.3f})')
        if np.any(mask_s):
            ax.scatter(np.log10(all_g_bar[mask_s]), resid_obs_all[mask_s],
                       c='firebrick', s=20, alpha=0.5, label=f'obs spiral (σ={np.std(resid_obs_all[mask_s]):.3f})')
            ax.scatter(np.log10(all_g_bar[mask_s]), resid_sidm_all[mask_s],
                       c='firebrick', marker='x', s=25, alpha=0.7,
                       label=f'SIDM spiral (σ={np.std(resid_sidm_all[mask_s]):.3f})')

        ax.axhline(0, color='k', ls='--', lw=0.5)
        ax.set_xlabel(r'$\log_{10}\, g_{\rm bar}$ [m/s$^2$]')
        ax.set_ylabel(r'$\log_{10}(g / g_{\rm McGaugh})$')
        ax.set_title(f'{bp_name} — Residuals from McGaugh RAR')
        ax.legend(fontsize=7)
        ax.set_ylim(-1.0, 1.0)

        fig.tight_layout()
        fig_path = os.path.join(OUT_DIR, f'rar_v2_{bp_name}.png')
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)
        print(f"  Plot saved → {fig_path}")

        # ── Fit summary CSV ──
        fit_csv = os.path.join(OUT_DIR, f'rar_v2_fits_{bp_name}.csv')
        with open(fit_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['galaxy', 'category', 'f_gas', 'upsilon', 'chi2_dof', 'r_1_kpc'])
            for fr in fit_results:
                w.writerow([fr['name'], fr['category'], fr['f_gas'],
                            f"{fr['upsilon']:.4f}", f"{fr['chi2_dof']:.3f}", f"{fr['r_1']:.4f}"])
        print(f"  Fit results → {fit_csv}")

    print(f"\n✓ §7.3 RAR v2 (gas-fraction fix) DONE.")


if __name__ == '__main__':
    main()
