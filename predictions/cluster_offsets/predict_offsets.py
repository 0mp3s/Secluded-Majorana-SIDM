#!/usr/bin/env python3
"""
predictions/cluster_offsets/predict_offsets.py
==============================================
Compute sigma/m at galaxy-cluster merger velocities and compare
to observational upper bounds (Harvey+2015 and individual systems).

Physics:
  During a cluster merger the DM halos pass through each other.
  Self-interactions create a drag force that separates DM from
  the collisionless galaxies:
      a_drag ~ -(sigma/m) * rho_DM * v^2

  The resulting DM-galaxy offset is proportional to sigma/m at
  the infall velocity v ~ 1000-5000 km/s.
  Harvey+2015 stacked 72 merging clusters and derived
      sigma/m < 0.47 cm^2/g  (95 % CL).

  Our velocity-dependent cross section MUST satisfy this bound.
  Because sigma/m drops at high v (resonant regime), this is
  a natural outcome of the model.

Output: sigma/m(v) curve for each BP + cluster data points.
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
from global_config import GC
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

_DIR = os.path.dirname(os.path.abspath(__file__))


def load_clusters(csv_path):
    clusters = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            clusters.append({
                'name':     row['system'],
                'v_infall': float(row['v_infall_km_s']),
                'v_err':    float(row['v_infall_err']),
                'sigma_m_upper': float(row['sigma_m_upper_cm2_g']),
                'ref':      row['ref'],
            })
    return clusters


def main():
    cfg = load_config(__file__)
    csv_path = os.path.join(_DIR, cfg.get('clusters_csv', 'clusters_data.csv'))
    out_dir  = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = GC.benchmarks_from_labels(cfg['benchmark_labels'])
    v_lo, v_hi = cfg.get('velocity_range_km_s', [500, 5000])
    n_v = cfg.get('n_velocity_points', 200)

    clusters = load_clusters(csv_path)

    print("=" * 80)
    print("  Cluster-Merger Offsets: sigma/m at High Velocities")
    print("=" * 80)

    # ---- Table per BP ----
    for bp in bps:
        label     = bp['label']
        m_chi     = bp['m_chi_GeV']
        m_phi_GeV = bp['m_phi_MeV'] / 1000.0
        alpha     = bp['alpha']
        lam       = alpha * m_chi / m_phi_GeV

        print(f"\n  --- {label}: m_chi={m_chi:.1f} GeV, "
              f"m_phi={bp['m_phi_MeV']:.2f} MeV, "
              f"alpha={alpha:.3e}, lambda={lam:.1f} ---")
        print(f"  {'System':<22s} {'v_inf':>7s} {'sigma/m_pred':>13s} "
              f"{'sigma/m_upper':>14s} {'Status':>8s}")
        print("  " + "-" * 70)

        all_ok = True
        for cl in clusters:
            sm = sigma_T_vpm(m_chi, m_phi_GeV, alpha, cl['v_infall'])
            ok = sm < cl['sigma_m_upper']
            if not ok:
                all_ok = False
            print(f"  {cl['name']:<22s} {cl['v_infall']:>7.0f} "
                  f"{sm:>13.4f} {cl['sigma_m_upper']:>14.2f} "
                  f"{'PASS' if ok else 'FAIL':>8s}")

        print(f"\n  => {label}: {'ALL PASS' if all_ok else 'SOME FAIL'}")

    # ---- Plot: sigma/m(v) + upper bounds ----
    v_arr = np.linspace(v_lo, v_hi, n_v)

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = {'BP1': 'steelblue', 'BP16': 'seagreen', 'MAP': 'firebrick'}

    for bp in bps:
        label     = bp['label']
        m_chi     = bp['m_chi_GeV']
        m_phi_GeV = bp['m_phi_MeV'] / 1000.0
        alpha     = bp['alpha']

        sm_arr = [sigma_T_vpm(m_chi, m_phi_GeV, alpha, float(v)) for v in v_arr]
        ax.plot(v_arr, sm_arr, lw=2, color=colors.get(label, 'gray'),
                label=label)

    # Cluster upper bounds as markers
    for cl in clusters:
        ax.errorbar(cl['v_infall'], cl['sigma_m_upper'],
                    xerr=cl['v_err'], fmt='v', ms=8, color='black',
                    capsize=3, zorder=5)
        ax.annotate(cl['name'].replace('_', ' '),
                    (cl['v_infall'], cl['sigma_m_upper']),
                    textcoords='offset points', xytext=(5, 5),
                    fontsize=7, rotation=15)

    ax.axhline(0.47, color='gray', ls=':', lw=1, alpha=0.6)
    ax.text(v_hi * 0.85, 0.50, 'Harvey+15 (0.47)', fontsize=8,
            color='gray', ha='right')

    ax.set_xlabel(r'$v_{\rm rel}$ (km s$^{-1}$)', fontsize=12)
    ax.set_ylabel(r'$\sigma/m$ (cm$^2$ g$^{-1}$)', fontsize=12)
    ax.set_title('SIDM Cross Section at Cluster Merger Velocities', fontsize=13)
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-4)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'cluster_sigma_m.png')
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved: {fig_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
