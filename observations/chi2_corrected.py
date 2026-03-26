#!/usr/bin/env python3
"""
observations/chi2_corrected.py
================================
MB-corrected χ² decomposition — does tension survive velocity averaging?

Physics:
  chi2_decomposition.py uses σ(v₀) — the point-estimate cross section.
  averaging_error.py showed this introduces 10–16% systematic at clusters.
  This script replaces σ(v₀) with ⟨σ/m⟩_MB in the χ² calculation:

      Δχ²_corrected = [(⟨σ/m⟩_MB(v₀) − central) / σ_err]²

  If tension disappears → it was a systematic artifact.
  If tension persists → it's real physics requiring new ingredients.

  Key test: NGC 1560 pull in MAP_relic was −0.70.
  After MB correction, does it become |pull| < 0.5 ?

Produces:
  - Console: side-by-side comparison (point vs averaged)
  - output/chi2_corrected.png  (comparison bar chart)
  - output/chi2_corrected.pdf
  - output/chi2_corrected.csv
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

from concurrent.futures import ThreadPoolExecutor
from scipy.interpolate import interp1d

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm

sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

cfg = load_config(__file__)
_bps = {bp['label']: bp for bp in
        GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1', 'MAP']))}

OBSERVATIONS = GC.observations_as_tuples()

N_WORKERS = os.cpu_count() or 4
N_DENSE = 500
V_DENSE = np.logspace(0, np.log10(25000.0), N_DENSE)


def _build_sigma_interp(m_chi, m_phi_GeV, alpha):
    """Dense grid + cubic log-log interpolator for σ(v)."""
    def _eval(v):
        return sigma_T_vpm(m_chi, m_phi_GeV, alpha, v)
    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        sig_dense = np.array(list(pool.map(_eval, V_DENSE)))
    valid = sig_dense > 0
    if np.sum(valid) < 10:
        return lambda v: 0.0
    log_interp = interp1d(np.log(V_DENSE[valid]), np.log(sig_dense[valid]),
                          kind='cubic', fill_value='extrapolate')
    return lambda v: np.exp(log_interp(np.log(v)))


def sigma_m_averaged(sigma_func, v_char_km_s, n_gauss=30):
    """⟨σ/m⟩ over MB distribution."""
    v0 = v_char_km_s / np.sqrt(2.0)
    v_max = max(5.0 * v_char_km_s, 100.0)
    a, b = 1.0, v_max
    nodes, weights_gl = np.polynomial.legendre.leggauss(n_gauss)
    v_pts = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    jac = 0.5 * (b - a)
    sigma_vals = sigma_func(v_pts)
    w_vals = v_pts**3 * np.exp(-(v_pts / v0)**2 / 2.0)
    num = jac * np.sum(weights_gl * sigma_vals * w_vals)
    den = jac * np.sum(weights_gl * w_vals)
    return num / den if den > 1e-300 else 0.0


def _compute_chi2_row(theory, name, v, central, lo, hi, ref):
    """Compute single-observation chi2 contribution."""
    one_sided = (lo == 0.0)
    if one_sided and theory <= hi:
        return dict(name=name, v=v, central=central, lo=lo, hi=hi, ref=ref,
                    theory=theory, sigma_err=hi - central if hi > central else 0.5 * max(central, 0.01),
                    pull=0.0, delta_chi2=0.0, one_sided=True)
    if theory >= central:
        sigma_err = hi - central if hi > central else 0.5 * central
    else:
        sigma_err = central - lo if central > lo else 0.5 * central
    if sigma_err <= 0:
        sigma_err = 0.5 * max(central, 0.01)
    pull = (theory - central) / sigma_err
    return dict(name=name, v=v, central=central, lo=lo, hi=hi, ref=ref,
                theory=theory, sigma_err=sigma_err,
                pull=pull, delta_chi2=pull ** 2, one_sided=False)


def main():
    t0 = time.time()
    hdr = "=" * 80
    print(hdr)
    print("  MB-Corrected χ² Decomposition")
    print("  σ(v₀)  vs  ⟨σ/m⟩_MB  — does tension survive averaging?")
    print(hdr)

    all_results = {}

    for label, bp in _bps.items():
        mc = bp['m_chi_GeV']
        mp = bp['m_phi_MeV']
        al = bp['alpha']
        mp_GeV = mp * 1e-3
        lam = al * mc / mp_GeV

        print(f"\n  {label}:  m_χ={mc:.2f} GeV, m_φ={mp:.2f} MeV, "
              f"α={al:.4e}, λ={lam:.2f}")

        # Build interpolator
        t1 = time.time()
        sigma_func = _build_sigma_interp(mc, mp_GeV, al)
        dt_interp = time.time() - t1

        rows_point = []
        rows_avg = []

        for name, v, central, lo, hi, ref in OBSERVATIONS:
            # Point estimate
            theory_pt = sigma_T_vpm(mc, mp_GeV, al, float(v))
            rows_point.append(
                _compute_chi2_row(theory_pt, name, v, central, lo, hi, ref))

            # MB averaged
            theory_avg = sigma_m_averaged(sigma_func, float(v))
            rows_avg.append(
                _compute_chi2_row(theory_avg, name, v, central, lo, hi, ref))

        chi2_pt = sum(r['delta_chi2'] for r in rows_point)
        chi2_avg = sum(r['delta_chi2'] for r in rows_avg)

        print(f"    Dense grid: {dt_interp:.1f}s")
        print(f"\n    {'System':<24s} {'v₀':>5s}  {'σ_pt':>8s}  {'pull_pt':>8s}  "
              f"{'Δχ²_pt':>7s}  │  {'⟨σ⟩_MB':>8s}  {'pull_MB':>8s}  {'Δχ²_MB':>7s}  Δ")
        print("    " + "-" * 100)

        for rp, ra in zip(rows_point, rows_avg):
            delta = ra['delta_chi2'] - rp['delta_chi2']
            arrow = "↑" if delta > 0.01 else "↓" if delta < -0.01 else "≈"
            print(f"    {rp['name']:<24s} {rp['v']:>5.0f}  {rp['theory']:>8.4f}  "
                  f"{rp['pull']:>+8.3f}  {rp['delta_chi2']:>7.4f}  │  "
                  f"{ra['theory']:>8.4f}  {ra['pull']:>+8.3f}  "
                  f"{ra['delta_chi2']:>7.4f}  {arrow}{abs(delta):>6.3f}")

        print(f"\n    Total χ²:  point={chi2_pt:.4f}   MB-corrected={chi2_avg:.4f}   "
              f"Δ={chi2_avg - chi2_pt:+.4f}")

        all_results[label] = (rows_point, rows_avg)

    # ──────────────────────────────────────────────────────────
    #  Figure: comparison bar chart
    # ──────────────────────────────────────────────────────────
    n_bp = len(all_results)
    fig, axes = plt.subplots(n_bp, 2, figsize=(14, 4.5 * n_bp))
    if n_bp == 1:
        axes = axes.reshape(1, -1)

    for idx, (label, (rows_pt, rows_avg)) in enumerate(all_results.items()):
        names = [r['name'] for r in rows_pt]
        pulls_pt = [r['pull'] for r in rows_pt]
        pulls_avg = [r['pull'] for r in rows_avg]
        dchi2_pt = [r['delta_chi2'] for r in rows_pt]
        dchi2_avg = [r['delta_chi2'] for r in rows_avg]
        y = np.arange(len(names))

        # Panel 1: Δχ² comparison
        ax1 = axes[idx, 0]
        ax1.barh(y - 0.2, dchi2_pt, 0.35, label='σ(v₀)', color='#d62728', alpha=0.8)
        ax1.barh(y + 0.2, dchi2_avg, 0.35, label='⟨σ/m⟩_MB', color='#2ca02c', alpha=0.8)
        ax1.set_yticks(y)
        ax1.set_yticklabels(names, fontsize=8)
        ax1.set_xlabel('Δχ²', fontsize=11)
        ax1.set_title(f'{label} — Δχ² comparison', fontsize=12)
        ax1.legend(fontsize=9)
        ax1.invert_yaxis()

        # Panel 2: pull comparison
        ax2 = axes[idx, 1]
        ax2.barh(y - 0.2, pulls_pt, 0.35, label='σ(v₀)', color='#d62728', alpha=0.8)
        ax2.barh(y + 0.2, pulls_avg, 0.35, label='⟨σ/m⟩_MB', color='#2ca02c', alpha=0.8)
        ax2.set_yticks(y)
        ax2.set_yticklabels(names, fontsize=8)
        ax2.set_xlabel('Pull', fontsize=11)
        ax2.set_title(f'{label} — pull comparison', fontsize=12)
        ax2.axvline(0, color='k', lw=0.8)
        ax2.axvline(-1, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax2.axvline(+1, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax2.legend(fontsize=9)
        ax2.invert_yaxis()

    plt.tight_layout()
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'chi2_corrected.{ext}'), dpi=200)
    plt.close(fig)
    print(f"\n  Figures saved to output/chi2_corrected.{{png,pdf}}")

    # ──────────────────────────────────────────────────────────
    #  CSV
    # ──────────────────────────────────────────────────────────
    csv_path = os.path.join(out_dir, 'chi2_corrected.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['benchmark', 'system', 'v_km_s',
                     'sigma_point', 'pull_point', 'dchi2_point',
                     'sigma_avg', 'pull_avg', 'dchi2_avg',
                     'dchi2_change'])
        for label, (rows_pt, rows_avg) in all_results.items():
            for rp, ra in zip(rows_pt, rows_avg):
                w.writerow([label, rp['name'], rp['v'],
                            f"{rp['theory']:.6e}", f"{rp['pull']:.6f}", f"{rp['delta_chi2']:.6f}",
                            f"{ra['theory']:.6e}", f"{ra['pull']:.6f}", f"{ra['delta_chi2']:.6f}",
                            f"{ra['delta_chi2'] - rp['delta_chi2']:.6f}"])
    print(f"  CSV saved to {csv_path}")

    dt = time.time() - t0
    print(f"\n  Total runtime: {dt:.1f}s")
    print(hdr)

    return all_results


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("✅ chi2_corrected done!")
    except Exception:
        pass
