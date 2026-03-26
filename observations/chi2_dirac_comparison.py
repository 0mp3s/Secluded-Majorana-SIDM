#!/usr/bin/env python3
"""
observations/chi2_dirac_comparison.py
=======================================
Majorana vs Dirac χ² — is the NGC 1560 tension from odd-l triplet weights?

Physics:
  Majorana scattering weights partial waves:
    w_even = 1 (singlet, CP-even)
    w_odd  = 3 (triplet, CP-odd)
  Dirac scattering treats all partial waves equally (w = 1).

  If the tension at v ≈ 50-60 km/s comes from the triplet-enhanced
  odd-l waves, then replacing Majorana → Dirac should REDUCE the pull.

  This script:
  1. Computes σ_T^{Maj}(v₀) and σ_T^{Dir}(v₀) at each observational velocity
  2. Computes Δχ² decomposition for both
  3. Compares pulls: is Majorana pull > Dirac pull at 50-60 km/s?
  4. Reports the "Majorana penalty" per observation

  Key diagnostic: if |pull_Maj| > |pull_Dir| at NGC 1560 / diverse RC,
  the triplet channel is the source of tension, pointing to real
  CP-odd physics rather than a generic fitting artifact.

Produces:
  - Console: side-by-side Majorana vs Dirac pulls & χ²
  - output/chi2_dirac_comparison.png  (comparison bar chart)
  - output/chi2_dirac_comparison.pdf
  - output/chi2_dirac_comparison.csv
"""
import sys, os, math, time, csv
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
sys.path.insert(0, os.path.join(_ROOT, 'cross_checks'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm
from numba import jit

# Import Dirac σ_T from sommerfeld_scattering_ratio module
from sommerfeld_scattering_ratio import sigma_T_dirac

# ── JIT warmup ──
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)
sigma_T_dirac(20.0, 10e-3, 1e-3, 100.0)

cfg = load_config(__file__)
_bps = {bp['label']: bp for bp in
        GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1', 'MAP', 'MAP_relic']))}
OBSERVATIONS = GC.observations_as_tuples()


def _compute_row(sigma_fn, label_type, name, v, central, lo, hi, ref):
    """Compute single-observation chi2 contribution."""
    theory = sigma_fn(v)
    one_sided = (lo == 0.0)
    if one_sided and theory <= hi:
        return dict(name=name, v=v, central=central, lo=lo, hi=hi, ref=ref,
                    theory=theory, pull=0.0, delta_chi2=0.0, type=label_type)
    if theory >= central:
        sigma_err = hi - central if hi > central else 0.5 * central
    else:
        sigma_err = central - lo if central > lo else 0.5 * central
    if sigma_err <= 0:
        sigma_err = 0.5 * max(central, 0.01)
    pull = (theory - central) / sigma_err
    return dict(name=name, v=v, central=central, lo=lo, hi=hi, ref=ref,
                theory=theory, pull=pull, delta_chi2=pull ** 2, type=label_type)


def main():
    t0 = time.time()
    hdr = "=" * 80
    print(hdr)
    print("  Majorana vs Dirac χ² — triplet-weight diagnostic")
    print(hdr)

    all_results = {}

    for label, bp in _bps.items():
        mc = bp['m_chi_GeV']
        mp_MeV = bp['m_phi_MeV']
        al = bp['alpha']
        mp_GeV = mp_MeV * 1e-3
        lam = al * mc / mp_GeV

        print(f"\n  {label}:  m_χ={mc:.2f} GeV, m_φ={mp_MeV:.2f} MeV, "
              f"α={al:.4e}, λ={lam:.2f}")

        def sigma_maj(v, _mc=mc, _mp=mp_GeV, _al=al):
            return sigma_T_vpm(_mc, _mp, _al, float(v))

        def sigma_dir(v, _mc=mc, _mp=mp_GeV, _al=al):
            return sigma_T_dirac(_mc, _mp, _al, float(v))

        rows_maj = []
        rows_dir = []

        for name, v, central, lo, hi, ref in OBSERVATIONS:
            rows_maj.append(
                _compute_row(sigma_maj, 'Majorana', name, v, central, lo, hi, ref))
            rows_dir.append(
                _compute_row(sigma_dir, 'Dirac', name, v, central, lo, hi, ref))

        chi2_maj = sum(r['delta_chi2'] for r in rows_maj)
        chi2_dir = sum(r['delta_chi2'] for r in rows_dir)
        penalty = chi2_maj - chi2_dir

        print(f"\n    {'System':<24s} {'v₀':>5s}  │  {'σ_Maj':>8s}  {'pull_M':>8s}  "
              f"{'Δχ²_M':>7s}  │  {'σ_Dir':>8s}  {'pull_D':>8s}  {'Δχ²_D':>7s}  │  penalty")
        print("    " + "-" * 105)

        for rm, rd in zip(rows_maj, rows_dir):
            pen = rm['delta_chi2'] - rd['delta_chi2']
            arrow = "⚠" if pen > 0.1 else " "
            print(f"    {rm['name']:<24s} {rm['v']:>5.0f}  │  {rm['theory']:>8.4f}  "
                  f"{rm['pull']:>+8.3f}  {rm['delta_chi2']:>7.4f}  │  "
                  f"{rd['theory']:>8.4f}  {rd['pull']:>+8.3f}  "
                  f"{rd['delta_chi2']:>7.4f}  │  {pen:>+7.4f} {arrow}")

        sign = "WORSE" if penalty > 0 else "BETTER"
        print(f"\n    Total χ²:  Majorana={chi2_maj:.4f}   Dirac={chi2_dir:.4f}    "
              f"Majorana penalty={penalty:+.4f}  ({sign})")
        print(f"    σ_Maj/σ_Dir ratio range: "
              f"[{min(rm['theory']/rd['theory'] if rd['theory'] > 0 else float('inf') for rm, rd in zip(rows_maj, rows_dir)):.3f}, "
              f"{max(rm['theory']/rd['theory'] if rd['theory'] > 0 else float('inf') for rm, rd in zip(rows_maj, rows_dir)):.3f}]")

        all_results[label] = (rows_maj, rows_dir, penalty)

    # ──────────────────────────────────────────────────────────
    #  Figure: Majorana vs Dirac comparison
    # ──────────────────────────────────────────────────────────
    n_bp = len(all_results)
    fig, axes = plt.subplots(n_bp, 2, figsize=(14, 4.5 * n_bp))
    if n_bp == 1:
        axes = axes.reshape(1, -1)

    for idx, (label, (rows_m, rows_d, pen)) in enumerate(all_results.items()):
        names = [r['name'] for r in rows_m]
        y = np.arange(len(names))

        # Panel 1: Δχ² comparison
        dchi2_m = [r['delta_chi2'] for r in rows_m]
        dchi2_d = [r['delta_chi2'] for r in rows_d]
        ax1 = axes[idx, 0]
        ax1.barh(y - 0.2, dchi2_m, 0.35, label='Majorana', color='#d62728', alpha=0.8)
        ax1.barh(y + 0.2, dchi2_d, 0.35, label='Dirac', color='#1f77b4', alpha=0.8)
        ax1.set_yticks(y)
        ax1.set_yticklabels(names, fontsize=8)
        ax1.set_xlabel('Δχ²', fontsize=11)
        sign_str = f"pen={pen:+.3f}" if abs(pen) > 0.001 else "≈0"
        ax1.set_title(f'{label} — Δχ² ({sign_str})', fontsize=12)
        ax1.legend(fontsize=9)
        ax1.invert_yaxis()

        # Panel 2: pull comparison
        pulls_m = [r['pull'] for r in rows_m]
        pulls_d = [r['pull'] for r in rows_d]
        ax2 = axes[idx, 1]
        ax2.barh(y - 0.2, pulls_m, 0.35, label='Majorana', color='#d62728', alpha=0.8)
        ax2.barh(y + 0.2, pulls_d, 0.35, label='Dirac', color='#1f77b4', alpha=0.8)
        ax2.set_yticks(y)
        ax2.set_yticklabels(names, fontsize=8)
        ax2.set_xlabel('Pull', fontsize=11)
        ax2.set_title(f'{label} — pulls: Majorana vs Dirac', fontsize=12)
        ax2.axvline(0, color='k', lw=0.8)
        ax2.axvline(-1, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax2.axvline(+1, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax2.legend(fontsize=9)
        ax2.invert_yaxis()

    plt.tight_layout()
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'chi2_dirac_comparison.{ext}'), dpi=200)
    plt.close(fig)
    print(f"\n  Figures saved to output/chi2_dirac_comparison.{{png,pdf}}")

    # ──────────────────────────────────────────────────────────
    #  CSV
    # ──────────────────────────────────────────────────────────
    csv_path = os.path.join(out_dir, 'chi2_dirac_comparison.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['benchmark', 'system', 'v_km_s',
                     'sigma_maj', 'pull_maj', 'dchi2_maj',
                     'sigma_dirac', 'pull_dirac', 'dchi2_dirac',
                     'penalty'])
        for label, (rows_m, rows_d, _) in all_results.items():
            for rm, rd in zip(rows_m, rows_d):
                w.writerow([label, rm['name'], rm['v'],
                            f"{rm['theory']:.6e}", f"{rm['pull']:.6f}", f"{rm['delta_chi2']:.6f}",
                            f"{rd['theory']:.6e}", f"{rd['pull']:.6f}", f"{rd['delta_chi2']:.6f}",
                            f"{rm['delta_chi2'] - rd['delta_chi2']:.6f}"])
    print(f"  CSV saved to {csv_path}")

    dt = time.time() - t0
    print(f"\n  Total runtime: {dt:.1f}s")
    print(hdr)

    return all_results


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("✅ chi2_dirac_comparison done!")
    except Exception:
        pass
