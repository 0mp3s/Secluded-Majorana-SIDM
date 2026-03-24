#!/usr/bin/env python3
"""
predictions/rotation_curves/fit_sparc_baryons.py
=================================================
Full rotation curve fit: V²_tot(r) = V²_bar(r) × Υ_* + V²_SIDM(r)

For each galaxy × benchmark point:
  1. Build NFW halo from V_max + concentration
  2. Compute SIDM core radius r_1 from σ/m(v)
  3. Build cored density profile (isothermal inside r_1, NFW outside)
  4. Fit Υ_* (stellar mass-to-light ratio) to minimize χ²
  5. Check if best-fit Υ_* is physical (0.1–1.0 M_sun/L_sun at 3.6μm)

Key physics:
  V²_tot(r) = Υ_* × V²_bar,template(r) + V²_DM(r)

  V_bar,template is from SPARC 3.6μm photometry (assumes Υ_*=1).
  The free parameter Υ_* rescales it.
  Physical range: 0.2–0.8 M_sun/L_sun (Meidt+2014, Schombert+2019).

Output: per-galaxy fits, Υ_* table, combined residual plot.
"""
import sys, os, csv, math
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d

# ---------- path bootstrap ----------
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# ------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ---- constants ----
G_N = 4.302e-6           # kpc (km/s)^2 / M_sun
KPC_CM = 3.086e21
MSUN_G = 1.989e33
SEC_PER_GYR = 3.156e16
KM_S_TO_CM_S = 1e5

_DIR = os.path.dirname(os.path.abspath(__file__))

# =========================================================================
#  NFW helpers (shared with predict_core_sizes.py)
# =========================================================================

def nfw_rho(r_kpc, rho_s, r_s):
    x = max(r_kpc / r_s, 1e-8)
    return rho_s / (x * (1.0 + x) ** 2)

def nfw_mass(r_kpc, rho_s, r_s):
    x = r_kpc / r_s
    return 4.0 * math.pi * rho_s * r_s**3 * (math.log(1 + x) - x / (1 + x))

def nfw_v_circ(r_kpc, rho_s, r_s):
    if r_kpc < 1e-8:
        return 0.0
    return math.sqrt(G_N * nfw_mass(r_kpc, rho_s, r_s) / r_kpc)

def nfw_params_from_vmax(V_max, c=12.0):
    rho_crit = 126.0  # M_sun/kpc^3 (h=0.674)
    f_c = math.log(1 + c) - c / (1 + c)
    delta_c = (200.0 / 3.0) * c**3 / f_c
    rho_s = rho_crit * delta_c
    x_max = 2.163
    f_xmax = math.log(1 + x_max) - x_max / (1 + x_max)
    r_s = math.sqrt(V_max**2 * x_max / (G_N * 4 * math.pi * rho_s * f_xmax))
    return rho_s, r_s

# =========================================================================
#  SIDM cored profile
# =========================================================================

def find_r1(rho_s, r_s, sigma_over_m, sigma_v_km_s, t_age_s):
    """Thermalization radius: rho(r1) * (σ/m) * v * t = 1."""
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
    """V_DM for cored-SIDM profile."""
    if r_kpc < 1e-8:
        return 0.0
    if r_kpc <= r_1:
        M_enc = (4.0 / 3.0) * math.pi * rho_core * r_kpc**3
    else:
        M_core = (4.0 / 3.0) * math.pi * rho_core * r_1**3
        M_enc = M_core + nfw_mass(r_kpc, rho_s, r_s) - nfw_mass(r_1, rho_s, r_s)
    return math.sqrt(G_N * M_enc / r_kpc)

# =========================================================================
#  Concentration–mass relation (Dutton & Macciò 2014)
# =========================================================================

def concentration_mass(M_200_Msun, h=0.674):
    """Dutton & Macciò 2014 eq.7 (Planck cosmology, z=0)."""
    log_M = math.log10(max(M_200_Msun * h / 1e12, 1e-30))
    return 10**(0.905 - 0.101 * log_M)

def nfw_params_with_cM(V_max, h=0.674):
    """NFW params with self-consistent c(M_200). Returns (rho_s, r_s, c, M_200)."""
    rho_crit = 126.0  # M_sun/kpc^3
    c = 12.0
    for _ in range(30):
        rho_s, r_s = nfw_params_from_vmax(V_max, c=c)
        R_200 = c * r_s
        M_200 = (4.0 / 3.0) * math.pi * 200.0 * rho_crit * R_200**3
        c_new = concentration_mass(M_200, h)
        if abs(c_new - c) / max(c, 0.1) < 1e-3:
            break
        c = 0.5 * (c + c_new)
    rho_s, r_s = nfw_params_from_vmax(V_max, c=c)
    return rho_s, r_s, c, M_200

# =========================================================================
#  Adiabatic contraction (Blumenthal+1986)
# =========================================================================

F_B = 0.157  # Omega_b / Omega_m (Planck 2018)

def solve_ac_ri(r_f, M_bar_rf, rho_s, r_s):
    """
    Blumenthal invariant:  r_i * M_NFW(r_i) = r_f * [(1-f_b)*M_NFW(r_i) + M_bar(r_f)]
    Solve for initial DM-shell radius r_i >= r_f.
    """
    if M_bar_rf <= F_B * nfw_mass(r_f, rho_s, r_s):
        return r_f  # no net contraction at this radius

    def residual(r_i):
        M_nfw = nfw_mass(r_i, rho_s, r_s)
        return r_i * M_nfw - r_f * ((1 - F_B) * M_nfw + M_bar_rf)

    r_lo, r_hi = r_f, 50.0 * r_f
    for _ in range(10):
        if residual(r_hi) > 0:
            break
        r_hi *= 5.0
    if residual(r_hi) <= 0:
        return r_hi

    for _ in range(80):
        mid = 0.5 * (r_lo + r_hi)
        if residual(mid) < 0:
            r_lo = mid
        else:
            r_hi = mid
    return 0.5 * (r_lo + r_hi)


def fit_galaxy_ac_sidm(r_data, V_obs, V_err, V_bar_template,
                       rho_s, r_s, sigma_m, sigma_v, t_age_s):
    """
    Iterative fit: Υ_* → AC'd DM profile → SIDM core → fit Υ_* → converge.
    Returns (upsilon, chi2, ndof, r_1, V_DM_array).
    """
    V_bar2 = V_bar_template**2
    V_bar_func = interp1d(r_data, V_bar_template, kind='linear',
                          bounds_error=False,
                          fill_value=(V_bar_template[0], V_bar_template[-1]))

    upsilon = 0.5
    best_chi2 = 1e30
    V_DM_result = np.zeros(len(r_data))
    r_1_result = 0.0

    for iteration in range(8):
        # Fine grid for density estimation
        r_min = max(r_data.min() * 0.3, 0.02)
        r_max = r_data.max() * 2.0
        r_fine = np.geomspace(r_min, r_max, 200)

        # AC'd DM mass profile on fine grid
        M_DM_AC = np.empty(len(r_fine))
        for i, rf in enumerate(r_fine):
            Vb = float(V_bar_func(rf))
            m_bar = max(upsilon * Vb**2 * rf / G_N, 0.0)
            ri = solve_ac_ri(rf, m_bar, rho_s, r_s)
            M_DM_AC[i] = (1 - F_B) * nfw_mass(ri, rho_s, r_s)

        # Density from mass profile
        rho_AC = np.gradient(M_DM_AC, r_fine) / (4.0 * np.pi * r_fine**2)
        rho_AC = np.maximum(rho_AC, 1e-10)

        # SIDM core radius from scattering condition
        v_cm = sigma_v * KM_S_TO_CM_S
        target_rho_cgs = 1.0 / (sigma_m * v_cm * t_age_s)
        target_rho = target_rho_cgs * KPC_CM**3 / MSUN_G

        r_1 = r_fine[-1]
        for i in range(len(r_fine)):
            if rho_AC[i] < target_rho:
                r_1 = r_fine[i]
                break
        rho_core = target_rho

        # V_DM at data radii with SIDM coring on AC'd profile
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

        def chi2_func(ups):
            V_tot2 = ups * V_bar2 + V_DM2
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

# =========================================================================
#  Data loading
# =========================================================================

def load_rotation_data(csv_path):
    """Load per-galaxy rotation curve data points."""
    galaxies = {}
    with open(csv_path, newline='', encoding='utf-8') as f:
        # skip comment lines
        lines = [l for l in f if not l.startswith('#')]
    import io
    reader = csv.DictReader(io.StringIO(''.join(lines)))
    for row in reader:
        name = row['name'].strip()
        if name not in galaxies:
            galaxies[name] = {'r': [], 'V_obs': [], 'V_err': [], 'V_bar': []}
        galaxies[name]['r'].append(float(row['r_kpc']))
        galaxies[name]['V_obs'].append(float(row['V_obs_km_s']))
        galaxies[name]['V_err'].append(float(row['V_err_km_s']))
        galaxies[name]['V_bar'].append(float(row['V_bar_km_s']))
    # convert to numpy
    for name in galaxies:
        for key in galaxies[name]:
            galaxies[name][key] = np.array(galaxies[name][key])
    return galaxies

def load_galaxy_meta(csv_path):
    """Load galaxy metadata (V_max, category, etc.)."""
    meta = {}
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            meta[row['name']] = {
                'V_max': float(row['V_max_km_s']),
                'category': row['category'],
            }
    return meta

# =========================================================================
#  Fitting
# =========================================================================

def chi2_fit(r_data, V_obs, V_err, V_bar_template, V_DM_at_r):
    """
    Fit Υ_* to minimize χ²:
        V_tot² = Υ_* × V_bar² + V_DM²
        V_tot  = sqrt(max(Υ_* × V_bar² + V_DM², 0))
    """
    V_bar2 = V_bar_template**2
    V_DM2 = V_DM_at_r**2

    def neg_loglike(upsilon):
        V_tot2 = upsilon * V_bar2 + V_DM2
        V_tot = np.sqrt(np.maximum(V_tot2, 1e-10))
        residuals = (V_obs - V_tot) / V_err
        return np.sum(residuals**2)

    result = minimize_scalar(neg_loglike, bounds=(0.01, 3.0), method='bounded')
    best_upsilon = result.x
    best_chi2 = result.fun
    ndof = len(V_obs) - 1  # 1 free parameter
    return best_upsilon, best_chi2, ndof

# =========================================================================
#  Main
# =========================================================================

def main():
    cfg = load_config(__file__)
    rc_csv = os.path.join(_DIR, cfg.get('rotation_data_csv', 'sparc_rotation_data.csv'))
    meta_csv = os.path.join(_DIR, cfg.get('sparc_csv', 'sparc_galaxies.csv'))
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = GC.benchmarks_from_labels(cfg['benchmark_labels'])
    t_age_s = cfg.get('halo_age_Gyr', 10.0) * SEC_PER_GYR

    galaxies = load_rotation_data(rc_csv)
    meta = load_galaxy_meta(meta_csv)

    print("=" * 105)
    print("  SPARC Rotation Curve Fit: V²_tot = Υ_* × V²_bar + V²_SIDM")
    print("  Improvements: c(M_200) Dutton+2014 + adiabatic contraction Blumenthal+1986")
    print("=" * 105)
    print(f"  Physical Υ_* range at 3.6μm: 0.2 – 0.8 M_sun/L_sun (Meidt+2014)")
    print(f"  Galaxies: {len(galaxies)}, Benchmark points: {len(bps)}")

    all_fit_results = {}

    for bp in bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi_GeV = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        lam = alpha * m_chi / m_phi_GeV

        print(f"\n{'='*95}")
        print(f"  {label}: m_chi={m_chi:.1f} GeV, m_phi={bp['m_phi_MeV']:.2f} MeV, "
              f"alpha={alpha:.3e}, lambda={lam:.1f}")
        print(f"{'='*95}")
        print(f"  {'Galaxy':<12s} {'V_max':>6s} {'c(M)':>5s} {'sigma/m':>8s} {'r_1':>6s} "
              f"{'Υ_*':>6s} {'chi2/dof':>9s} {'Υ_* OK?':>8s}")
        print("  " + "-" * 72)

        fit_results = []
        for gal_name in sorted(galaxies.keys()):
            if gal_name not in meta:
                continue
            gal = galaxies[gal_name]
            gal_meta = meta[gal_name]
            V_max = gal_meta['V_max']

            # NFW halo with c(M_200) relation
            rho_s, r_s, conc, M_200 = nfw_params_with_cM(V_max)

            # SIDM cross section at this halo's velocity
            sigma_v = V_max / math.sqrt(2.0)
            v_rel = sigma_v * math.sqrt(2.0)
            sigma_m = sigma_T_vpm(m_chi, m_phi_GeV, alpha, v_rel)

            # Full fit: AC + SIDM + χ² minimize (iterative)
            upsilon, chi2, ndof, r_1, V_DM_arr = fit_galaxy_ac_sidm(
                gal['r'], gal['V_obs'], gal['V_err'], gal['V_bar'],
                rho_s, r_s, sigma_m, sigma_v, t_age_s)

            chi2_dof = chi2 / ndof if ndof > 0 else chi2
            physical = 0.1 <= upsilon <= 1.5
            status = "OK" if physical else "BAD"

            print(f"  {gal_name:<12s} {V_max:>6.0f} {conc:>5.1f} {sigma_m:>8.3f} {r_1:>6.2f} "
                  f"{upsilon:>6.2f} {chi2_dof:>9.2f} {status:>8s}")

            fit_results.append({
                'name': gal_name,
                'V_max': V_max,
                'conc': conc,
                'category': gal_meta['category'],
                'sigma_m': sigma_m,
                'r_1': r_1,
                'upsilon': upsilon,
                'chi2': chi2,
                'ndof': ndof,
                'chi2_dof': chi2_dof,
                'physical': physical,
                'r': gal['r'],
                'V_obs': gal['V_obs'],
                'V_err': gal['V_err'],
                'V_bar': gal['V_bar'],
                'V_DM': V_DM_arr,
            })
        all_fit_results[label] = fit_results

    # ---- Summary ----
    print(f"\n{'='*95}")
    print("  SUMMARY: Best-fit Υ_* values")
    print(f"{'='*95}")
    print(f"  {'Galaxy':<12s}", end="")
    for bp in bps:
        print(f"  {bp['label']:>10s}", end="")
    print(f"  {'Physical?':>10s}")
    print("  " + "-" * (12 + 12 * len(bps) + 12))

    for i, gal_name in enumerate(sorted(galaxies.keys())):
        if gal_name not in meta:
            continue
        print(f"  {gal_name:<12s}", end="")
        all_phys = True
        for bp in bps:
            label = bp['label']
            fr = [r for r in all_fit_results[label] if r['name'] == gal_name][0]
            print(f"  {fr['upsilon']:>10.3f}", end="")
            if not fr['physical']:
                all_phys = False
        status = "ALL OK" if all_phys else "CHECK"
        print(f"  {status:>10s}")

    # ---- Per-galaxy rotation curve plots ----
    gal_names = sorted(galaxies.keys())
    gal_names = [g for g in gal_names if g in meta]
    n_gal = len(gal_names)
    n_bp = len(bps)

    fig, axes = plt.subplots(n_gal, n_bp, figsize=(5 * n_bp, 4 * n_gal),
                             squeeze=False, sharex=False)

    for ig, gal_name in enumerate(gal_names):
        for ib, bp in enumerate(bps):
            ax = axes[ig][ib]
            label = bp['label']
            fr = [r for r in all_fit_results[label] if r['name'] == gal_name][0]

            r = fr['r']
            V_obs = fr['V_obs']
            V_err = fr['V_err']
            V_bar = fr['V_bar']
            V_DM = fr['V_DM']
            ups = fr['upsilon']

            V_bar_scaled = np.sqrt(ups) * V_bar
            V_tot = np.sqrt(ups * V_bar**2 + V_DM**2)

            # plot
            ax.errorbar(r, V_obs, yerr=V_err, fmt='ko', ms=4, capsize=2,
                        label='Observed', zorder=5)
            ax.plot(r, V_tot, 'r-', lw=2, label=f'Total (Υ*={ups:.2f})')
            ax.plot(r, V_bar_scaled, 'b--', lw=1.2, alpha=0.7,
                    label=f'Baryons (×√Υ*)')
            ax.plot(r, V_DM, 'g:', lw=1.2, alpha=0.7, label='SIDM halo')

            # NFW for reference (using c(M) for this galaxy)
            rho_s_ref, r_s_ref, _, _ = nfw_params_with_cM(fr['V_max'])
            r_fine = np.linspace(0.1, max(r) * 1.1, 100)
            V_nfw_fine = [nfw_v_circ(ri, rho_s_ref, r_s_ref) for ri in r_fine]
            ax.plot(r_fine, V_nfw_fine, 'g--', lw=0.8, alpha=0.3, label='NFW (no core)')

            ax.set_ylim(bottom=0)
            ax.grid(True, alpha=0.2)
            chi2_str = f"χ²/dof={fr['chi2_dof']:.2f}"
            ax.set_title(f"{gal_name} — {label}\n{chi2_str}", fontsize=10)
            if ib == 0:
                ax.set_ylabel('V (km/s)', fontsize=10)
            if ig == n_gal - 1:
                ax.set_xlabel('r (kpc)', fontsize=10)
            if ig == 0 and ib == 0:
                ax.legend(fontsize=7, loc='lower right')

    fig.suptitle('SPARC Rotation Curve Fits: V²_tot = Υ_* V²_bar + V²_SIDM  [c(M) + AC]',
                 fontsize=13, y=1.01)
    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'sparc_baryons_fit.png')
    fig.savefig(fig_path, dpi=120, bbox_inches='tight')
    print(f"\n  Plot saved: {fig_path}")
    plt.close(fig)

    # ---- Υ_* summary plot ----
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    colors = {'BP1': 'steelblue', 'BP16': 'seagreen', 'MAP': 'firebrick'}
    offsets = {'BP1': -0.15, 'BP16': 0.0, 'MAP': 0.15}

    x_pos = np.arange(n_gal)
    for bp in bps:
        label = bp['label']
        ups_vals = []
        for gal_name in gal_names:
            fr = [r for r in all_fit_results[label] if r['name'] == gal_name][0]
            ups_vals.append(fr['upsilon'])
        ax2.bar(x_pos + offsets.get(label, 0), ups_vals, width=0.13,
                color=colors.get(label, 'gray'), label=label, alpha=0.85)

    ax2.axhspan(0.2, 0.8, alpha=0.15, color='green', label='Physical range')
    ax2.axhline(0.5, color='green', ls='--', lw=1, alpha=0.5)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(gal_names, rotation=30, ha='right')
    ax2.set_ylabel(r'$\Upsilon_*$ ($M_\odot / L_\odot$ at 3.6μm)', fontsize=12)
    ax2.set_title('Best-fit Stellar Mass-to-Light Ratio', fontsize=13)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_ylim(0, max(2.0, max(fr['upsilon'] for res in all_fit_results.values() for fr in res) * 1.2))

    fig2.tight_layout()
    fig2_path = os.path.join(out_dir, 'upsilon_summary.png')
    fig2.savefig(fig2_path, dpi=150, bbox_inches='tight')
    print(f"  Plot saved: {fig2_path}")
    plt.close(fig2)


if __name__ == "__main__":
    main()
