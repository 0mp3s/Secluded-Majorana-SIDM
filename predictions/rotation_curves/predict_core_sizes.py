#!/usr/bin/env python3
"""
predictions/rotation_curves/predict_core_sizes.py
==================================================
Predict SIDM core sizes for SPARC-like galaxies and compare
the resulting V(2 kpc) to observed values.

Physics (Kaplinghat+2016, Kamada+2017):
  SIDM thermalizes the inner halo within a radius r_1 where
      rho_NFW(r_1) × (sigma/m) × sigma_v × t_age ~ 1

  Inside r_1 → isothermal core with constant density rho_core.
  Outside r_1 → NFW profile is unchanged.

  The "diversity problem" (Oman+2015): at fixed V_max, observed
  galaxies show wide scatter in V(2 kpc).  NFW predicts a tight
  relation.  Velocity-dependent SIDM can produce diversity
  because sigma/m varies with halo velocity.

Output: V(2 kpc) vs V_max plot + table.
"""
import sys, os, csv, math
import numpy as np

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
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ---- constants ----
G_N = 4.302e-6          # kpc (km/s)^2 / M_sun
KPC_CM = 3.086e21        # cm
MSUN_G = 1.989e33        # g
SEC_PER_GYR = 3.156e16
KM_S_TO_CM_S = 1e5

_DIR = os.path.dirname(os.path.abspath(__file__))

# =========================================================================
#  NFW helpers
# =========================================================================

def nfw_rho(r_kpc, rho_s, r_s):
    """NFW density [M_sun/kpc^3]."""
    x = max(r_kpc / r_s, 1e-8)
    return rho_s / (x * (1.0 + x) ** 2)


def nfw_mass(r_kpc, rho_s, r_s):
    """NFW enclosed mass [M_sun]."""
    x = r_kpc / r_s
    return 4.0 * math.pi * rho_s * r_s ** 3 * (math.log(1.0 + x) - x / (1.0 + x))


def nfw_v_circ(r_kpc, rho_s, r_s):
    """NFW circular velocity [km/s]."""
    if r_kpc < 1e-8:
        return 0.0
    return math.sqrt(G_N * nfw_mass(r_kpc, rho_s, r_s) / r_kpc)


def nfw_params_from_vmax(V_max, c=12.0):
    """
    Return (rho_s [M_sun/kpc^3], r_s [kpc]) for an NFW halo with the
    given V_max and concentration c.
    """
    rho_crit = 126.0  # h=0.674 → M_sun/kpc^3
    f_c = math.log(1.0 + c) - c / (1.0 + c)
    delta_c = (200.0 / 3.0) * c ** 3 / f_c
    rho_s = rho_crit * delta_c

    x_max = 2.163
    f_xmax = math.log(1.0 + x_max) - x_max / (1.0 + x_max)
    r_s = math.sqrt(V_max ** 2 * x_max
                    / (G_N * 4.0 * math.pi * rho_s * f_xmax))
    return rho_s, r_s


# =========================================================================
#  SIDM core
# =========================================================================

def find_r1(rho_s, r_s, sigma_over_m, sigma_v_km_s, t_age_s):
    """
    Thermalization radius r_1 [kpc] defined by
        rho(r_1) * (σ/m) * v * t_age = 1.
    """
    v_cm = sigma_v_km_s * KM_S_TO_CM_S
    target_rho_cgs = 1.0 / (sigma_over_m * v_cm * t_age_s)
    target_rho = target_rho_cgs * KPC_CM ** 3 / MSUN_G   # → M_sun/kpc^3

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
    """Circular velocity of a cored-SIDM + NFW profile [km/s]."""
    if r_kpc < 1e-8:
        return 0.0
    if r_kpc <= r_1:
        M_enc = (4.0 / 3.0) * math.pi * rho_core * r_kpc ** 3
    else:
        M_core = (4.0 / 3.0) * math.pi * rho_core * r_1 ** 3
        M_enc = M_core + nfw_mass(r_kpc, rho_s, r_s) - nfw_mass(r_1, rho_s, r_s)
    return math.sqrt(G_N * M_enc / r_kpc)


# =========================================================================
#  I/O
# =========================================================================

def load_sparc(csv_path):
    gals = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gals.append({
                'name':       row['name'],
                'V_max':      float(row['V_max_km_s']),
                'V_2kpc_obs': float(row['V_2kpc_km_s']),
                'V_2kpc_err': float(row['V_2kpc_err']),
                'r_core_obs': float(row['r_core_kpc']),
                'r_core_err': float(row['r_core_err_kpc']),
                'sigma_v':    float(row['sigma_v_km_s']),
                'category':   row['category'],
            })
    return gals


# =========================================================================
#  main
# =========================================================================

def main():
    cfg = load_config(__file__)
    csv_path  = os.path.join(_DIR, cfg.get('sparc_csv', 'sparc_subset.csv'))
    out_dir   = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps  = cfg['benchmark_points']
    conc = cfg.get('nfw_concentration', 12.0)
    t_age_s = cfg.get('halo_age_Gyr', 10.0) * SEC_PER_GYR

    galaxies = load_sparc(csv_path)

    print("=" * 90)
    print("  SIDM Core-Size Predictions  vs  SPARC Rotation Curves")
    print("=" * 90)

    all_results = {}

    for bp in bps:
        label     = bp['label']
        m_chi     = bp['m_chi_GeV']
        m_phi_GeV = bp['m_phi_MeV'] / 1000.0
        alpha     = bp['alpha']
        lam       = alpha * m_chi / m_phi_GeV

        print(f"\n  --- {label}: m_chi={m_chi:.1f} GeV, m_phi={bp['m_phi_MeV']:.2f} MeV, "
              f"alpha={alpha:.3e}, lambda={lam:.1f} ---")
        print(f"  {'Galaxy':<12s} {'V_max':>6s} {'sigma/m':>9s} {'r_1':>7s} "
              f"{'r_c_obs':>7s} {'V2_pred':>8s} {'V2_obs':>7s} {'Match':>6s}")
        print("  " + "-" * 70)

        results = []
        for gal in galaxies:
            rho_s, r_s = nfw_params_from_vmax(gal['V_max'], c=conc)

            v_rel = gal['sigma_v'] * math.sqrt(2.0)
            sigma_m = sigma_T_vpm(m_chi, m_phi_GeV, alpha, v_rel)

            r_1 = find_r1(rho_s, r_s, sigma_m, gal['sigma_v'], t_age_s)
            rho_core = nfw_rho(r_1, rho_s, r_s)

            V2_nfw  = nfw_v_circ(2.0, rho_s, r_s)
            V2_sidm = sidm_v_circ(2.0, rho_s, r_s, r_1, rho_core)

            diff  = abs(V2_sidm - gal['V_2kpc_obs'])
            match = "OK" if diff < 2.0 * gal['V_2kpc_err'] else "MISS"

            print(f"  {gal['name']:<12s} {gal['V_max']:>6.0f} {sigma_m:>9.3f} "
                  f"{r_1:>7.2f} {gal['r_core_obs']:>7.1f} "
                  f"{V2_sidm:>8.1f} {gal['V_2kpc_obs']:>7.1f} {match:>6s}")

            results.append({
                'name':      gal['name'],
                'V_max':     gal['V_max'],
                'sigma_m':   sigma_m,
                'r_1':       r_1,
                'V2_nfw':    V2_nfw,
                'V2_sidm':   V2_sidm,
                'V2_obs':    gal['V_2kpc_obs'],
                'V2_err':    gal['V_2kpc_err'],
                'category':  gal['category'],
            })
        all_results[label] = results

    # ------------------------------------------------------------------
    #  Plot: V(2 kpc) vs V_max
    # ------------------------------------------------------------------
    n_bp = len(bps)
    fig, axes = plt.subplots(1, n_bp, figsize=(6 * n_bp, 5), sharey=True)
    if n_bp == 1:
        axes = [axes]

    for ax, bp in zip(axes, bps):
        label   = bp['label']
        results = all_results[label]

        vm   = [r['V_max']  for r in results]
        v2o  = [r['V2_obs'] for r in results]
        v2e  = [r['V2_err'] for r in results]
        v2s  = [r['V2_sidm'] for r in results]

        # NFW curve
        vm_line = np.linspace(30, 260, 120)
        v2_nfw_line = []
        for v in vm_line:
            rs, rss = nfw_params_from_vmax(v, c=conc)
            v2_nfw_line.append(nfw_v_circ(2.0, rs, rss))

        ax.plot(vm_line, v2_nfw_line, 'k--', lw=1.5, alpha=0.5, label='NFW (no SIDM)')
        ax.errorbar(vm, v2o, yerr=v2e, fmt='ko', ms=5, capsize=3,
                    label='Observed', zorder=3)
        ax.scatter(vm, v2s, c='crimson', s=60, marker='D', edgecolors='k',
                   linewidths=0.5, label=f'SIDM ({label})', zorder=4)
        for i in range(len(results)):
            ax.plot([vm[i], vm[i]], [v2o[i], v2s[i]], 'r-', alpha=0.3, lw=1)

        lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000.0)
        ax.set_title(f'{label}  (λ = {lam:.0f})', fontsize=12)
        ax.set_xlabel(r'$V_{\rm max}$ (km s$^{-1}$)', fontsize=11)
        if ax is axes[0]:
            ax.set_ylabel(r'$V(2\,{\rm kpc})$ (km s$^{-1}$)', fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Rotation-Curve Diversity: SIDM Prediction vs SPARC Data',
                 fontsize=13, y=1.02)
    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'rotation_curve_diversity.png')
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved: {fig_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
