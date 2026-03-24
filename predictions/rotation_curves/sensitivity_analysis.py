#!/usr/bin/env python3
"""
predictions/rotation_curves/sensitivity_analysis.py
====================================================
Systematic sensitivity analysis for SPARC rotation curve fits.

Scans 4 uncertain astrophysical parameters on a grid:
  1. c_factor   — multiplicative scatter on c(M_200) relation   [0.6, 1.5]
  2. nu_AC      — adiabatic contraction strength (0=none, 1=full Blumenthal) [0, 1]
  3. t_age      — halo age in Gyr                               [5, 13]
  4. f_Vbar     — V_bar template systematic rescaling            [0.85, 1.15]

For each grid point × galaxy × BP:
  - Fit Υ_* (single free parameter)
  - Record (Υ_*, χ²/dof, r_1)

Output:
  - CSV results table (every grid point)
  - "Physical region" maps: which parameter combos give 0.2 ≤ Υ_* ≤ 0.8
  - Per-galaxy sensitivity heat maps
  - Summary: which galaxies are parameter-robust vs fragile

Uses multiprocessing on 14 CPU cores.
"""

import sys, os, csv, math, time, itertools
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---------- path bootstrap ----------
_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
# ------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm

# JIT warmup
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ---- constants ----
G_N          = 4.302e-6       # kpc (km/s)^2 / M_sun
KPC_CM       = 3.086e21
MSUN_G       = 1.989e33
SEC_PER_GYR  = 3.156e16
KM_S_TO_CM_S = 1e5
F_B_COSMO    = 0.157          # Omega_b / Omega_m  (Planck 2018)

_DIR = os.path.dirname(os.path.abspath(__file__))

# =====================================================================
#  NFW / SIDM / AC  physics  (copied from fit_sparc_baryons.py)
# =====================================================================

def nfw_mass(r, rho_s, r_s):
    x = r / r_s
    return 4 * math.pi * rho_s * r_s**3 * (math.log(1 + x) - x / (1 + x))

def nfw_rho(r, rho_s, r_s):
    x = max(r / r_s, 1e-8)
    return rho_s / (x * (1 + x)**2)

def nfw_params_from_vmax(V_max, c):
    rho_crit = 126.0
    f_c = math.log(1 + c) - c / (1 + c)
    delta_c = (200 / 3) * c**3 / f_c
    rho_s = rho_crit * delta_c
    x_max = 2.163
    f_xmax = math.log(1 + x_max) - x_max / (1 + x_max)
    r_s = math.sqrt(V_max**2 * x_max / (G_N * 4 * math.pi * rho_s * f_xmax))
    return rho_s, r_s

def concentration_mass(M_200, h=0.674):
    log_M = math.log10(max(M_200 * h / 1e12, 1e-30))
    return 10**(0.905 - 0.101 * log_M)

def self_consistent_c(V_max, c_factor=1.0):
    """Iterate c(M) relation with multiplicative scatter factor."""
    c = 12.0
    rho_crit = 126.0
    for _ in range(30):
        rho_s, r_s = nfw_params_from_vmax(V_max, c)
        R_200 = c * r_s
        M_200 = (4 / 3) * math.pi * 200 * rho_crit * R_200**3
        c_new = concentration_mass(M_200) * c_factor
        if abs(c_new - c) / max(c, 0.1) < 1e-3:
            break
        c = 0.5 * (c + c_new)
    return c

def solve_ac_ri(r_f, M_bar_rf, rho_s, r_s, nu):
    """AC with partial contraction: nu=0 → no AC, nu=1 → full Blumenthal."""
    if nu < 1e-6 or M_bar_rf <= F_B_COSMO * nfw_mass(r_f, rho_s, r_s):
        return r_f

    def residual(r_i):
        M_nfw = nfw_mass(r_i, rho_s, r_s)
        return r_i * M_nfw - r_f * ((1 - F_B_COSMO) * M_nfw + M_bar_rf)

    r_lo, r_hi = r_f, 50 * r_f
    for _ in range(10):
        if residual(r_hi) > 0:
            break
        r_hi *= 5
    if residual(r_hi) <= 0:
        return r_hi

    for _ in range(80):
        mid = 0.5 * (r_lo + r_hi)
        if residual(mid) < 0:
            r_lo = mid
        else:
            r_hi = mid
    r_i_full = 0.5 * (r_lo + r_hi)
    # Partial AC: interpolate between r_f (no AC) and r_i_full
    return r_f + nu * (r_i_full - r_f)

# =====================================================================
#  Single grid-point evaluation
# =====================================================================

def evaluate_one(args):
    """
    Evaluate a single (galaxy, BP, c_factor, nu_AC, t_age_Gyr, f_Vbar) point.
    Returns dict with all results.
    """
    (gal_name, r_data, V_obs, V_err, V_bar_template, V_max,
     bp_label, m_chi, m_phi_GeV, alpha,
     c_factor, nu_AC, t_age_Gyr, f_Vbar) = args

    t_age_s = t_age_Gyr * SEC_PER_GYR

    # Concentration
    c = self_consistent_c(V_max, c_factor)
    rho_s, r_s = nfw_params_from_vmax(V_max, c)

    # SIDM cross section
    sigma_v = V_max / math.sqrt(2)
    v_rel = sigma_v * math.sqrt(2)
    sigma_m = sigma_T_vpm(m_chi, m_phi_GeV, alpha, v_rel)

    # Rescale V_bar template
    V_bar = V_bar_template * f_Vbar
    V_bar2 = V_bar**2

    V_bar_func = interp1d(r_data, V_bar, kind='linear',
                          bounds_error=False,
                          fill_value=(V_bar[0], V_bar[-1]))

    upsilon = 0.5
    V_DM = np.zeros(len(r_data))
    r_1 = 0.0

    for iteration in range(8):
        r_min = max(r_data.min() * 0.3, 0.02)
        r_max = r_data.max() * 2.0
        r_fine = np.geomspace(r_min, r_max, 200)

        M_DM_AC = np.empty(len(r_fine))
        for i, rf in enumerate(r_fine):
            Vb = float(V_bar_func(rf))
            m_bar = max(upsilon * Vb**2 * rf / G_N, 0.0)
            ri = solve_ac_ri(rf, m_bar, rho_s, r_s, nu_AC)
            M_DM_AC[i] = (1 - F_B_COSMO) * nfw_mass(ri, rho_s, r_s)

        rho_AC = np.gradient(M_DM_AC, r_fine) / (4 * np.pi * r_fine**2)
        rho_AC = np.maximum(rho_AC, 1e-10)

        v_cm = sigma_v * KM_S_TO_CM_S
        target_rho_cgs = 1.0 / (sigma_m * v_cm * t_age_s)
        target_rho = target_rho_cgs * KPC_CM**3 / MSUN_G

        r_1 = r_fine[-1]
        for i in range(len(r_fine)):
            if rho_AC[i] < target_rho:
                r_1 = r_fine[i]
                break
        rho_core = target_rho

        M_DM_AC_at_r1 = float(np.interp(r_1, r_fine, M_DM_AC))
        M_DM_AC_data = np.interp(r_data, r_fine, M_DM_AC)

        for i, r in enumerate(r_data):
            if r <= r_1:
                M_enc = (4 / 3) * math.pi * rho_core * r**3
            else:
                M_core = (4 / 3) * math.pi * rho_core * r_1**3
                M_enc = M_core + (M_DM_AC_data[i] - M_DM_AC_at_r1)
            V_DM[i] = math.sqrt(max(G_N * M_enc / r, 0))

        V_DM2 = V_DM**2

        def chi2_func(ups):
            V_tot2 = ups * V_bar2 + V_DM2
            V_tot = np.sqrt(np.maximum(V_tot2, 1e-10))
            return np.sum(((V_obs - V_tot) / V_err)**2)

        result = minimize_scalar(chi2_func, bounds=(0.01, 3.0), method='bounded')
        new_upsilon = result.x
        chi2 = result.fun

        if abs(new_upsilon - upsilon) < 0.01:
            upsilon = new_upsilon
            break
        upsilon = new_upsilon

    ndof = max(len(V_obs) - 1, 1)
    chi2_dof = chi2 / ndof
    physical = 0.2 <= upsilon <= 0.8

    return {
        'galaxy': gal_name,
        'bp': bp_label,
        'c_factor': c_factor,
        'nu_AC': nu_AC,
        't_age_Gyr': t_age_Gyr,
        'f_Vbar': f_Vbar,
        'c_actual': c,
        'sigma_m': sigma_m,
        'r_1': r_1,
        'upsilon': upsilon,
        'chi2': chi2,
        'chi2_dof': chi2_dof,
        'physical': physical,
    }

# =====================================================================
#  Data loading  (same as fit_sparc_baryons.py)
# =====================================================================

def load_rotation_data(csv_path):
    galaxies = {}
    with open(csv_path, newline='', encoding='utf-8') as f:
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
    for name in galaxies:
        for key in galaxies[name]:
            galaxies[name][key] = np.array(galaxies[name][key])
    return galaxies

def load_galaxy_meta(csv_path):
    meta = {}
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            meta[row['name']] = {'V_max': float(row['V_max_km_s']),
                                 'category': row['category']}
    return meta

# =====================================================================
#  Plotting
# =====================================================================

def plot_sensitivity_heatmaps(results, gal_names, bp_labels, out_dir):
    """
    For each galaxy × BP: 2D heat maps of Υ_* in (c_factor, nu_AC) space,
    marginalized over t_age and f_Vbar at their central values.
    """
    # Identify central t_age and f_Vbar
    t_ages = sorted(set(r['t_age_Gyr'] for r in results))
    f_vbars = sorted(set(r['f_Vbar'] for r in results))
    t_mid = t_ages[len(t_ages) // 2]
    f_mid = f_vbars[len(f_vbars) // 2]

    c_factors = sorted(set(r['c_factor'] for r in results))
    nu_ACs = sorted(set(r['nu_AC'] for r in results))

    n_gal = len(gal_names)
    n_bp = len(bp_labels)

    fig, axes = plt.subplots(n_gal, n_bp, figsize=(5 * n_bp, 3.5 * n_gal),
                             squeeze=False)

    for ig, gal in enumerate(gal_names):
        for ib, bp in enumerate(bp_labels):
            ax = axes[ig][ib]
            subset = [r for r in results
                      if r['galaxy'] == gal and r['bp'] == bp
                      and abs(r['t_age_Gyr'] - t_mid) < 0.1
                      and abs(r['f_Vbar'] - f_mid) < 0.01]

            Z = np.full((len(nu_ACs), len(c_factors)), np.nan)
            for r in subset:
                ic = c_factors.index(r['c_factor'])
                inu = nu_ACs.index(r['nu_AC'])
                Z[inu, ic] = r['upsilon']

            im = ax.pcolormesh(c_factors, nu_ACs, Z,
                               cmap='RdYlGn_r', vmin=0, vmax=1.5, shading='nearest')
            # Physical region contour
            ax.contour(c_factors, nu_ACs, Z, levels=[0.2, 0.8],
                       colors=['blue', 'blue'], linewidths=1.5, linestyles=['--', '--'])
            ax.contourf(c_factors, nu_ACs, Z, levels=[0.2, 0.8],
                        colors=['blue'], alpha=0.08)

            ax.set_title(f'{gal} — {bp}', fontsize=10)
            if ig == n_gal - 1:
                ax.set_xlabel('c_factor', fontsize=9)
            if ib == 0:
                ax.set_ylabel(r'$\nu_{\rm AC}$', fontsize=10)

    fig.suptitle(r'$\Upsilon_*$ sensitivity  (t_age={:.0f} Gyr, f_Vbar={:.2f})'
                 .format(t_mid, f_mid), fontsize=13, y=1.01)
    cbar = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02,
                        label=r'$\Upsilon_*$')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'sensitivity_heatmaps.png'),
                dpi=130, bbox_inches='tight')
    plt.close(fig)

def plot_physical_fraction(results, gal_names, bp_labels, out_dir):
    """
    Bar chart: per galaxy × BP, fraction of grid points with physical Υ_*.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(gal_names))
    width = 0.25
    colors_bp = {'BP1': 'steelblue', 'BP16': 'seagreen', 'MAP': 'firebrick'}

    for ib, bp in enumerate(bp_labels):
        fracs = []
        for gal in gal_names:
            subset = [r for r in results if r['galaxy'] == gal and r['bp'] == bp]
            if subset:
                fracs.append(sum(1 for r in subset if r['physical']) / len(subset))
            else:
                fracs.append(0)
        ax.bar(x + (ib - 1) * width, fracs, width=width,
               color=colors_bp.get(bp, 'gray'), label=bp, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(gal_names, rotation=30, ha='right')
    ax.set_ylabel('Fraction with physical $\\Upsilon_*$  (0.2–0.8)', fontsize=11)
    ax.set_title('Parameter Space Coverage: Physical $\\Upsilon_*$', fontsize=13)
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color='gray', ls=':', lw=1)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, axis='y')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'physical_fraction.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)

def plot_upsilon_vs_param(results, gal_names, bp_labels, out_dir):
    """1D slices: Υ_* vs each parameter, other params at central values."""
    params = [
        ('c_factor', 'c_factor', 'c / c(M) factor'),
        ('nu_AC',    'nu_AC',    'AC strength $\\nu$'),
        ('t_age_Gyr','t_age_Gyr','Halo age (Gyr)'),
        ('f_Vbar',   'f_Vbar',   '$f_{V_{\\rm bar}}$'),
    ]
    # central values
    all_vals = {p[0]: sorted(set(r[p[0]] for r in results)) for p in params}
    central = {p[0]: all_vals[p[0]][len(all_vals[p[0]]) // 2] for p in params}

    colors_gal = plt.cm.tab10(np.linspace(0, 1, len(gal_names)))

    for bp in bp_labels:
        fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
        for ip, (pkey, _, plabel) in enumerate(params):
            ax = axes[ip]
            for ig, gal in enumerate(gal_names):
                # Fix all other params at central, vary this one
                subset = [r for r in results
                          if r['galaxy'] == gal and r['bp'] == bp
                          and all(abs(r[pk] - central[pk]) < 0.01
                                  for pk, _, _ in params if pk != pkey)]
                if not subset:
                    continue
                subset.sort(key=lambda r: r[pkey])
                xs = [r[pkey] for r in subset]
                ys = [r['upsilon'] for r in subset]
                ax.plot(xs, ys, 'o-', color=colors_gal[ig], ms=4, lw=1.3,
                        label=gal if ip == 0 else None)

            ax.axhspan(0.2, 0.8, alpha=0.1, color='green')
            ax.axhline(0.5, color='green', ls='--', lw=0.8, alpha=0.5)
            ax.set_xlabel(plabel, fontsize=11)
            if ip == 0:
                ax.set_ylabel(r'$\Upsilon_*$', fontsize=12)
            ax.set_title(plabel, fontsize=10)
            ax.grid(True, alpha=0.2)

        axes[0].legend(fontsize=7, loc='upper left', ncol=2)
        fig.suptitle(f'{bp}: $\\Upsilon_*$ sensitivity to each parameter '
                     f'(others at central)', fontsize=13, y=1.02)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f'sensitivity_1D_{bp}.png'),
                    dpi=130, bbox_inches='tight')
        plt.close(fig)

# =====================================================================
#  Main
# =====================================================================

def main():
    t0 = time.time()
    cfg = load_config(__file__)
    rc_csv = os.path.join(_DIR, cfg.get('rotation_data_csv', 'sparc_rotation_data.csv'))
    meta_csv = os.path.join(_DIR, cfg.get('sparc_csv', 'sparc_galaxies.csv'))
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = GC.benchmarks_from_labels(cfg['benchmark_labels'])
    galaxies = load_rotation_data(rc_csv)
    meta = load_galaxy_meta(meta_csv)
    gal_names = sorted(g for g in galaxies if g in meta)

    # ---- Parameter grid ----
    c_factors  = np.round(np.linspace(0.6, 1.5, 10), 3)     # 10 values
    nu_ACs     = np.round(np.linspace(0.0, 1.0, 6), 3)      #  6 values
    t_ages     = np.round(np.linspace(5.0, 13.0, 5), 1)     #  5 values
    f_Vbars    = np.round(np.linspace(0.85, 1.15, 4), 3)    #  4 values

    grid = list(itertools.product(c_factors, nu_ACs, t_ages, f_Vbars))
    total_jobs = len(gal_names) * len(bps) * len(grid)

    print("=" * 90)
    print("  SPARC Rotation Curve — Sensitivity Analysis")
    print("=" * 90)
    print(f"  Galaxies: {len(gal_names)}  |  BPs: {len(bps)}  |  Grid: {len(grid)} combos")
    print(f"  Total evaluations: {total_jobs:,}")
    print(f"  c_factor:  {c_factors[0]}  →  {c_factors[-1]}  ({len(c_factors)} pts)")
    print(f"  nu_AC:     {nu_ACs[0]}  →  {nu_ACs[-1]}  ({len(nu_ACs)} pts)")
    print(f"  t_age:     {t_ages[0]}  →  {t_ages[-1]}  ({len(t_ages)} pts)")
    print(f"  f_Vbar:    {f_Vbars[0]}  →  {f_Vbars[-1]}  ({len(f_Vbars)} pts)")
    print()

    # Build job list
    jobs = []
    for gal_name in gal_names:
        gal = galaxies[gal_name]
        V_max = meta[gal_name]['V_max']
        for bp in bps:
            m_chi = bp['m_chi_GeV']
            m_phi_GeV = bp['m_phi_MeV'] / 1000.0
            alpha = bp['alpha']
            for (cf, nu, ta, fv) in grid:
                jobs.append((
                    gal_name, gal['r'], gal['V_obs'], gal['V_err'],
                    gal['V_bar'], V_max,
                    bp['label'], m_chi, m_phi_GeV, alpha,
                    cf, nu, ta, fv
                ))

    # Parallel execution
    n_workers = min(12, os.cpu_count() or 1)
    print(f"  Launching {n_workers} workers for {total_jobs:,} jobs ...")
    results = []
    done = 0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(evaluate_one, job) for job in jobs]
        for future in as_completed(futures):
            results.append(future.result())
            done += 1
            if done % 500 == 0 or done == total_jobs:
                elapsed = time.time() - t0
                rate = done / elapsed
                eta = (total_jobs - done) / max(rate, 1)
                print(f"    [{done:>6,}/{total_jobs:,}]  "
                      f"{elapsed:.0f}s elapsed  |  ETA {eta:.0f}s  |  "
                      f"{rate:.1f} evals/s")

    elapsed = time.time() - t0
    print(f"\n  Done: {total_jobs:,} evaluations in {elapsed:.1f}s "
          f"({total_jobs / elapsed:.1f} evals/s)")

    # ---- Save CSV ----
    csv_path = os.path.join(out_dir, 'sensitivity_results.csv')
    fields = ['galaxy', 'bp', 'c_factor', 'nu_AC', 't_age_Gyr', 'f_Vbar',
              'c_actual', 'sigma_m', 'r_1', 'upsilon', 'chi2', 'chi2_dof', 'physical']
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(results, key=lambda x: (x['galaxy'], x['bp'],
                                                 x['c_factor'], x['nu_AC'])):
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r[k], float) else r[k])
                        for k in fields})
    print(f"  CSV saved: {csv_path}  ({len(results)} rows)")

    # ---- Summary table ----
    bp_labels = [bp['label'] for bp in bps]
    print(f"\n{'='*90}")
    print("  PHYSICAL Υ_* FRACTION  (0.2 ≤ Υ_* ≤ 0.8)")
    print(f"{'='*90}")
    print(f"  {'Galaxy':<12s}", end="")
    for bp in bp_labels:
        print(f"  {bp:>10s}", end="")
    print(f"  {'Overall':>10s}")
    print("  " + "-" * (12 + 12 * len(bp_labels) + 12))

    for gal in gal_names:
        print(f"  {gal:<12s}", end="")
        total_phys = 0
        total_count = 0
        for bp in bp_labels:
            sub = [r for r in results if r['galaxy'] == gal and r['bp'] == bp]
            n_phys = sum(1 for r in sub if r['physical'])
            frac = n_phys / len(sub) if sub else 0
            total_phys += n_phys
            total_count += len(sub)
            print(f"  {frac:>9.1%}  ", end="")
        overall = total_phys / total_count if total_count else 0
        print(f"  {overall:>9.1%}")

    # ---- "Sweet spot" analysis: best (c_factor, nu_AC) per galaxy ----
    print(f"\n{'='*90}")
    print("  BEST (c_factor, nu_AC) per galaxy × BP  (lowest |Υ_* - 0.5|)")
    print(f"{'='*90}")
    print(f"  {'Galaxy':<12s}  {'BP':<6s}  {'c_fac':>6s}  {'nu_AC':>6s}  "
          f"{'t_age':>6s}  {'f_Vb':>5s}  {'Υ_*':>6s}  {'χ²/dof':>7s}")
    print("  " + "-" * 70)

    for gal in gal_names:
        for bp in bp_labels:
            sub = [r for r in results if r['galaxy'] == gal and r['bp'] == bp]
            if not sub:
                continue
            best = min(sub, key=lambda r: abs(r['upsilon'] - 0.5))
            if best['physical']:
                marker = "✓"
            else:
                marker = "✗"
            print(f"  {gal:<12s}  {bp:<6s}  {best['c_factor']:>6.2f}  "
                  f"{best['nu_AC']:>6.2f}  {best['t_age_Gyr']:>6.1f}  "
                  f"{best['f_Vbar']:>5.2f}  {best['upsilon']:>6.3f}  "
                  f"{best['chi2_dof']:>7.2f}  {marker}")

    # ---- Chi2 quality: best chi2/dof among physical fits ----
    print(f"\n{'='*90}")
    print("  BEST χ²/dof AMONG PHYSICAL FITS")
    print(f"{'='*90}")
    for gal in gal_names:
        for bp in bp_labels:
            phys = [r for r in results
                    if r['galaxy'] == gal and r['bp'] == bp and r['physical']]
            if not phys:
                print(f"  {gal:<12s} {bp:<6s}  — NO PHYSICAL FITS —")
                continue
            best = min(phys, key=lambda r: r['chi2_dof'])
            print(f"  {gal:<12s} {bp:<6s}  Υ*={best['upsilon']:.3f}  "
                  f"χ²/dof={best['chi2_dof']:.2f}  "
                  f"c_fac={best['c_factor']:.2f}  nu={best['nu_AC']:.2f}  "
                  f"t={best['t_age_Gyr']:.0f}Gyr  fV={best['f_Vbar']:.2f}")

    # ---- Plots ----
    print("\n  Generating plots ...")
    plot_sensitivity_heatmaps(results, gal_names, bp_labels, out_dir)
    print("    → sensitivity_heatmaps.png")
    plot_physical_fraction(results, gal_names, bp_labels, out_dir)
    print("    → physical_fraction.png")
    plot_upsilon_vs_param(results, gal_names, bp_labels, out_dir)
    for bp in bp_labels:
        print(f"    → sensitivity_1D_{bp}.png")

    print(f"\n  Total wall time: {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()
