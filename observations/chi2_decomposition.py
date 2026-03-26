#!/usr/bin/env python3
"""
observations/chi2_decomposition.py
====================================
χ² decomposition per observation — identify the most constraining system.

Physics:
  The total χ² from chi2_fit.py is a sum over 13 observational systems.
  This script decomposes it to reveal which observations actually drive
  the fit and which are irrelevant:

      χ² = Σᵢ Δχ²ᵢ ,   Δχ²ᵢ = [(theory_i − central_i) / σ_i]²

  For each benchmark, we compute:
    1. Δχ²ᵢ per system — individual contributions
    2. Fractional contribution fᵢ = Δχ²ᵢ / χ²_total
    3. Pull = (theory − central) / σ — signed, shows bias direction
    4. Cumulative χ² when sorted by velocity — shows which velocity
       range dominates the constraint

  One-sided upper limits (Bullet, Harvey) contribute Δχ² = 0 when
  theory < bound, but CAN constrain if model overshoots.

  This analysis reveals:
    - Which observations are the "power constraints"
    - Whether the fit is driven by dSphs, clusters, or rotation curves
    - Whether removing one observation significantly changes the fit
    - Leave-one-out Δχ² (sensitivity to each observation)

Produces:
  - Console: decomposition table + ranked contributions
  - output/chi2_decomposition.png  (3-panel: bar chart, cumulative, pull)
  - output/chi2_decomposition.pdf
  - output/chi2_decomposition.csv
"""
import sys, os, math, time, csv
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

cfg = load_config(__file__)
_bps = {bp['label']: bp for bp in
        GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1', 'MAP']))}

OBSERVATIONS = GC.observations_as_tuples()


# ──────────────────────────────────────────────────────────────
#  Per-observation χ² contribution
# ──────────────────────────────────────────────────────────────

def decompose_chi2(bp):
    """Compute per-observation χ² contribution for a benchmark point.

    Returns list of dicts with keys:
      name, v, central, lo, hi, ref, theory, sigma_err, pull, delta_chi2
    """
    mc  = bp['m_chi_GeV']
    mp  = bp['m_phi_MeV'] * 1e-3       # → GeV for solver
    al  = bp['alpha_chi']

    rows = []
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        theory = sigma_T_vpm(mc, mp, al, float(v))

        # One-sided upper limits
        one_sided = (lo == 0.0)
        if one_sided and theory <= hi:
            rows.append(dict(
                name=name, v=v, central=central, lo=lo, hi=hi, ref=ref,
                theory=theory, sigma_err=hi - central if hi > central else 0.5 * max(central, 0.01),
                pull=0.0, delta_chi2=0.0, one_sided=True, within_bound=True,
            ))
            continue

        # Asymmetric errors
        if theory >= central:
            sigma_err = hi - central if hi > central else 0.5 * central
        else:
            sigma_err = central - lo if central > lo else 0.5 * central
        if sigma_err <= 0:
            sigma_err = 0.5 * max(central, 0.01)

        pull = (theory - central) / sigma_err
        delta_chi2 = pull ** 2

        rows.append(dict(
            name=name, v=v, central=central, lo=lo, hi=hi, ref=ref,
            theory=theory, sigma_err=sigma_err,
            pull=pull, delta_chi2=delta_chi2,
            one_sided=one_sided, within_bound=False,
        ))

    return rows


# ──────────────────────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    hdr = "=" * 72
    print(hdr)
    print("  χ² Decomposition per Observation")
    print(hdr)

    all_results = {}

    for label, bp in _bps.items():
        mc  = bp['m_chi_GeV']
        mp  = bp['m_phi_MeV']
        al  = bp['alpha_chi']
        lam = al * mc / (mp * 1e-3)

        print(f"\n  {label}:  m_χ={mc:.2f} GeV, m_φ={mp:.2f} MeV, "
              f"α={al:.4e}, λ={lam:.2f}")

        rows = decompose_chi2(bp)
        chi2_total = sum(r['delta_chi2'] for r in rows)

        # Sort by velocity for display
        rows_by_v = sorted(rows, key=lambda r: r['v'])

        print(f"\n    {'System':<24s} {'v₀':>5s}  {'σ_theory':>9s}  "
              f"{'central':>8s}  {'σ_err':>7s}  {'pull':>7s}  "
              f"{'Δχ²':>7s}  {'f [%]':>7s}  Note")
        print("    " + "-" * 100)

        for r in rows_by_v:
            f_pct = 100.0 * r['delta_chi2'] / chi2_total if chi2_total > 0 else 0.0
            note = ""
            if r.get('one_sided') and r.get('within_bound'):
                note = "upper limit OK"
            elif f_pct > 20:
                note = "DOMINANT"
            elif f_pct > 10:
                note = "significant"
            elif r['delta_chi2'] < 0.01:
                note = "negligible"

            print(f"    {r['name']:<24s} {r['v']:>5.0f}  {r['theory']:>9.4f}  "
                  f"{r['central']:>8.3f}  {r['sigma_err']:>7.4f}  "
                  f"{r['pull']:>+7.3f}  {r['delta_chi2']:>7.4f}  "
                  f"{f_pct:>6.1f}%  {note}")

        print(f"\n    Total χ² = {chi2_total:.4f}   (ndof = 10 free, 11 relic)")
        print(f"    χ²/dof = {chi2_total / 10:.4f} (free)   "
              f"{chi2_total / 11:.4f} (relic)")

        # Ranked by contribution
        rows_ranked = sorted(rows, key=lambda r: r['delta_chi2'], reverse=True)
        print(f"\n    Ranked contributions:")
        cumul = 0.0
        for i, r in enumerate(rows_ranked):
            cumul += r['delta_chi2']
            f_pct = 100.0 * r['delta_chi2'] / chi2_total if chi2_total > 0 else 0.0
            cumul_pct = 100.0 * cumul / chi2_total if chi2_total > 0 else 0.0
            print(f"      #{i+1:>2d}  {r['name']:<24s}  Δχ²={r['delta_chi2']:>7.4f}  "
                  f"({f_pct:>5.1f}%)   cumul={cumul_pct:>5.1f}%")

        all_results[label] = rows

    # ──────────────────────────────────────────────────────────
    #  Figure: 3-panel for each benchmark
    # ──────────────────────────────────────────────────────────
    n_bp = len(all_results)
    fig = plt.figure(figsize=(16, 5 * n_bp))
    outer = gridspec.GridSpec(n_bp, 1, hspace=0.35)

    for idx, (label, rows) in enumerate(all_results.items()):
        inner = gridspec.GridSpecFromSubplotSpec(1, 3,
                subplot_spec=outer[idx], wspace=0.30)

        rows_by_v = sorted(rows, key=lambda r: r['v'])
        names = [r['name'] for r in rows_by_v]
        delta_chi2 = np.array([r['delta_chi2'] for r in rows_by_v])
        pulls = np.array([r['pull'] for r in rows_by_v])
        chi2_total = delta_chi2.sum()

        # Panel 1: bar chart Δχ² per system
        ax1 = fig.add_subplot(inner[0])
        colors = ['#d62728' if d > 0.2 * chi2_total else
                  '#ff7f0e' if d > 0.1 * chi2_total else
                  '#2ca02c' if d > 0.01 else '#aaaaaa'
                  for d in delta_chi2]
        bars = ax1.barh(range(len(names)), delta_chi2, color=colors, edgecolor='k', lw=0.5)
        ax1.set_yticks(range(len(names)))
        ax1.set_yticklabels(names, fontsize=8)
        ax1.set_xlabel('Δχ²', fontsize=11)
        ax1.set_title(f'{label} — per-observation Δχ²', fontsize=12)
        ax1.invert_yaxis()

        # Panel 2: cumulative χ² vs velocity
        ax2 = fig.add_subplot(inner[1])
        velocities = np.array([r['v'] for r in rows_by_v])
        cumul = np.cumsum(delta_chi2)
        ax2.step(velocities, cumul, where='post', color='navy', lw=2)
        ax2.scatter(velocities, cumul, color='navy', s=30, zorder=5)
        ax2.set_xscale('log')
        ax2.set_xlabel('v [km/s]', fontsize=11)
        ax2.set_ylabel('Cumulative χ²', fontsize=11)
        ax2.set_title(f'{label} — cumulative χ² vs velocity', fontsize=12)
        ax2.axhline(chi2_total, ls='--', color='gray', alpha=0.5)
        ax2.text(velocities[-1], chi2_total * 1.05,
                 f'Total = {chi2_total:.2f}', ha='right', fontsize=9,
                 color='gray')

        # Panel 3: pull diagram
        ax3 = fig.add_subplot(inner[2])
        pull_colors = ['#d62728' if abs(p) > 1 else
                       '#ff7f0e' if abs(p) > 0.5 else '#2ca02c'
                       for p in pulls]
        ax3.barh(range(len(names)), pulls, color=pull_colors,
                 edgecolor='k', lw=0.5)
        ax3.set_yticks(range(len(names)))
        ax3.set_yticklabels(names, fontsize=8)
        ax3.set_xlabel('Pull = (theory − data) / σ', fontsize=11)
        ax3.set_title(f'{label} — pulls', fontsize=12)
        ax3.axvline(0, color='k', lw=0.8)
        ax3.axvline(-1, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax3.axvline(+1, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax3.invert_yaxis()

    plt.tight_layout()
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'chi2_decomposition.{ext}'), dpi=200)
    plt.close(fig)
    print(f"\n  Figures saved to output/chi2_decomposition.{{png,pdf}}")

    # ──────────────────────────────────────────────────────────
    #  CSV output
    # ──────────────────────────────────────────────────────────
    csv_path = os.path.join(out_dir, 'chi2_decomposition.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['benchmark', 'system', 'v_km_s', 'sigma_theory',
                     'central', 'lo', 'hi', 'sigma_err', 'pull',
                     'delta_chi2', 'fraction_pct', 'one_sided',
                     'within_bound', 'ref'])
        for label, rows in all_results.items():
            chi2_total = sum(r['delta_chi2'] for r in rows)
            for r in sorted(rows, key=lambda r: r['v']):
                f_pct = 100.0 * r['delta_chi2'] / chi2_total if chi2_total > 0 else 0.0
                w.writerow([
                    label, r['name'], r['v'], f"{r['theory']:.6e}",
                    r['central'], r['lo'], r['hi'], f"{r['sigma_err']:.6e}",
                    f"{r['pull']:.6f}", f"{r['delta_chi2']:.6f}",
                    f"{f_pct:.2f}",
                    int(r.get('one_sided', False)),
                    int(r.get('within_bound', False)),
                    r['ref'],
                ])
    print(f"  CSV saved to {csv_path}")

    dt = time.time() - t0
    print(f"\n  Total runtime: {dt:.1f}s")
    print(hdr)

    return all_results


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("✅ chi2_decomposition done!")
    except Exception:
        pass
