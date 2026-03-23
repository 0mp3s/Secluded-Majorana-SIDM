#!/usr/bin/env python3
"""
predictions/gravothermal/predict_gravothermal.py
=================================================
Compute gravothermal evolution timescales for classical dSphs
and check whether observed cores are consistent with our model.

Physics:
  The SIDM scattering rate per particle is:
      Gamma = (sigma/m) * rho_s * sigma_v
  where rho_s is the NFW scale density and sigma_v the 1D velocity dispersion.

  The number of scatterings per Hubble time is:
      N_scatter = Gamma * t_age

  Regimes (Essig+2019, Nishikawa+2020, Kaplinghat+2016):
      N_scatter < 1    --> too few scatterings, NFW cusp survives
      1 < N_scatter < ~100  --> core formation phase (observed cored dSphs)
      N_scatter > ~100-300  --> gravothermal collapse begins --> re-cusping

  Prediction: if lambda is large (MAP, lambda~333), N_scatter may exceed
  the collapse threshold in old dSphs --> inconsistent with observed cores.
  If lambda ~ 4 (BP1), N_scatter ~ 1-10 --> stable cores.

Output: table + bar chart comparing BP1/BP16/MAP predictions.
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

# Constants
MSUN_PER_KPC3_TO_GEV_PER_CM3 = 3.796e-2  # 1 M_sun/kpc^3 ~ 3.796e-2 GeV/cm^3
SEC_PER_GYR = 3.156e16
KM_S_TO_CM_S = 1e5
KPC_TO_CM = 3.086e21

_DIR = os.path.dirname(os.path.abspath(__file__))


def load_dsphs(csv_path):
    """Load dSph data from CSV."""
    dsphs = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            dsphs.append({
                'name': row['name'],
                'sigma_v': float(row['sigma_v_km_s']),
                'r_half_pc': float(row['r_half_pc']),
                'rho_s': float(row['rho_s_Msun_kpc3']),
                'r_s_kpc': float(row['r_s_kpc']),
                'core_observed': row['core_observed'],
                't_age_Gyr': float(row['t_age_Gyr']),
            })
    return dsphs


def compute_predictions(dsphs, benchmark_points):
    """
    For each dSph × benchmark, compute:
      - sigma/m at the dSph velocity
      - Gamma = (sigma/m) * rho_s * sigma_v
      - N_scatter = Gamma * t_age
    """
    results = []
    for dsph in dsphs:
        v = dsph['sigma_v']  # km/s (characteristic DM velocity ~ 1D dispersion)
        # relative velocity ~ sqrt(2) * sigma_v for MB distribution
        v_rel = v * math.sqrt(2)
        rho_s = dsph['rho_s']  # M_sun / kpc^3
        r_s = dsph['r_s_kpc']  # kpc
        t_age = dsph['t_age_Gyr']

        row = {'dsph': dsph['name'], 'v_rel_km_s': v_rel,
               'core_observed': dsph['core_observed']}

        for bp in benchmark_points:
            label = bp['label']
            m_chi = bp['m_chi_GeV']
            m_phi_GeV = bp['m_phi_MeV'] / 1000.0
            alpha = bp['alpha']

            # sigma/m in cm^2/g
            sigma_over_m = sigma_T_vpm(m_chi, m_phi_GeV, alpha, v_rel)

            # Convert units for rate:
            # Gamma = (sigma/m) [cm^2/g] * rho_DM [g/cm^3] * v_rel [cm/s]
            # rho_DM in g/cm^3:
            # 1 M_sun = 1.989e33 g, 1 kpc = 3.086e21 cm
            rho_g_cm3 = rho_s * 1.989e33 / (3.086e21)**3  # g/cm^3

            v_rel_cm_s = v_rel * KM_S_TO_CM_S

            # Gamma in 1/s
            gamma = sigma_over_m * rho_g_cm3 * v_rel_cm_s

            # N_scatter = Gamma * t_age
            t_age_s = t_age * SEC_PER_GYR
            n_scatter = gamma * t_age_s

            # lambda = 2*alpha*m_chi/m_phi_GeV
            lam = alpha * m_chi / m_phi_GeV

            row[f'{label}_sigma_m'] = sigma_over_m
            row[f'{label}_gamma_Gyr'] = gamma * SEC_PER_GYR
            row[f'{label}_N_scatter'] = n_scatter
            row[f'{label}_lambda'] = lam

        results.append(row)
    return results


def classify(n_scatter, threshold_core=1.0, threshold_collapse=100.0):
    if n_scatter < threshold_core:
        return "CUSPY (too few)"
    elif n_scatter < threshold_collapse:
        return "CORED (stable)"
    else:
        return "COLLAPSE (re-cusp)"


def main():
    cfg = load_config(__file__)
    csv_path = os.path.join(_DIR, cfg.get('dsphs_csv', 'dsphs_data.csv'))
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = cfg.get('benchmark_points', [
        {"label": "BP1",  "m_chi_GeV": 20.69, "m_phi_MeV": 11.34, "alpha": 1.048e-3},
        {"label": "BP16", "m_chi_GeV": 20.7,  "m_phi_MeV": 9.91,  "alpha": 1.048e-3},
        {"label": "MAP",  "m_chi_GeV": 90.64, "m_phi_MeV": 13.85, "alpha": 2.546e-2},
    ])
    threshold_core = cfg.get('core_formation_threshold', 1.0)
    threshold_collapse = cfg.get('collapse_threshold', 100.0)

    print("=" * 80)
    print("  Gravothermal Collapse Predictions for Classical dSphs")
    print("=" * 80)

    dsphs = load_dsphs(csv_path)
    results = compute_predictions(dsphs, bps)

    # Print results table per benchmark
    for bp in bps:
        label = bp['label']
        lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000.0)
        print(f"\n  --- {label}: m_chi={bp['m_chi_GeV']:.1f} GeV, "
              f"m_phi={bp['m_phi_MeV']:.2f} MeV, alpha={bp['alpha']:.3e}, lambda={lam:.1f} ---")
        print(f"  {'dSph':<14s} {'v_rel':>6s} {'sigma/m':>10s} {'N_scatter':>10s} "
              f"{'Prediction':>18s} {'Observed':>10s} {'Match':>6s}")
        print("  " + "-" * 78)

        for r in results:
            sigma_m = r[f'{label}_sigma_m']
            n_scat = r[f'{label}_N_scatter']
            pred = classify(n_scat, threshold_core, threshold_collapse)
            obs = r['core_observed']

            # Does prediction match observation?
            if obs == "YES":
                match = "OK" if "CORED" in pred else "FAIL"
            elif obs == "AMBIGUOUS":
                match = "~"
            else:
                match = "OK" if "CUSPY" in pred or "COLLAPSE" in pred else "?"

            print(f"  {r['dsph']:<14s} {r['v_rel_km_s']:>6.1f} {sigma_m:>10.4f} "
                  f"{n_scat:>10.2f} {pred:>18s} {obs:>10s} {match:>6s}")

    # Summary
    print(f"\n{'=' * 80}")
    print("  SUMMARY")
    print(f"{'=' * 80}")
    for bp in bps:
        label = bp['label']
        n_ok = 0
        n_fail = 0
        n_ambig = 0
        for r in results:
            n_scat = r[f'{label}_N_scatter']
            pred = classify(n_scat, threshold_core, threshold_collapse)
            obs = r['core_observed']
            if obs == "YES":
                if "CORED" in pred:
                    n_ok += 1
                else:
                    n_fail += 1
            else:
                n_ambig += 1
        print(f"  {label}: {n_ok} OK, {n_fail} FAIL, {n_ambig} ambiguous "
              f"(out of {len(results)} dSphs)")

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(12, 6))

    dsph_names = [r['dsph'] for r in results]
    x = np.arange(len(dsph_names))
    width = 0.25
    colors = {'BP1': 'steelblue', 'BP16': 'seagreen', 'MAP': 'firebrick'}

    for i, bp in enumerate(bps):
        label = bp['label']
        n_scats = [r[f'{label}_N_scatter'] for r in results]
        bars = ax.bar(x + i * width, n_scats, width, label=label,
                      color=colors.get(label, 'gray'), alpha=0.8)

    ax.axhline(threshold_core, color='green', ls='--', lw=1.5,
               label=f'Core formation (N={threshold_core})')
    ax.axhline(threshold_collapse, color='red', ls='--', lw=2,
               label=f'Collapse threshold (N={threshold_collapse})')

    ax.set_yscale('log')
    ax.set_ylabel('N_scatter = Γ × t_age', fontsize=12)
    ax.set_xlabel('Dwarf Spheroidal', fontsize=12)
    ax.set_title('Gravothermal Evolution: Scatterings per Hubble Time', fontsize=13)
    ax.set_xticks(x + width)
    ax.set_xticklabels(dsph_names, rotation=30, ha='right')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    # Shade regions
    ylim = ax.get_ylim()
    ax.axhspan(ylim[0], threshold_core, alpha=0.05, color='blue', label='_nolegend_')
    ax.axhspan(threshold_core, threshold_collapse, alpha=0.05, color='green', label='_nolegend_')
    ax.axhspan(threshold_collapse, ylim[1], alpha=0.08, color='red', label='_nolegend_')
    ax.set_ylim(ylim)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'gravothermal_prediction.png')
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved: {fig_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
