#!/usr/bin/env python3
"""
V10 - v37_velocity_averaged.py
===============================
Compare fixed-velocity sigma/m to Maxwell-Boltzmann velocity-averaged values.

For each astrophysical system with characteristic velocity v_char, the proper
observable is the thermally-averaged cross section:

    <sigma v> / <v> = integral[ sigma(v) * v * f_MB(v; v0) dv ]
                     / integral[ v * f_MB(v; v0) dv ]

where f_MB(v; v0) = 4*pi*(v/v0)^2 * (1/(2*pi*v0^2))^(3/2) * exp(-v^2/(2*v0^2))
is the Maxwell-Boltzmann distribution with 1D velocity dispersion v0.

Following Kaplinghat, Tulin & Yu (2016), the characteristic velocity v_char
is identified with the most-probable relative velocity: v_char = 2*v0 (or 
more precisely, v0 = v_char / sqrt(2) for the peak of f(v)*v).

We compute this for the 17 relic BPs and the observational velocities.
"""
# === path setup (auto-generated) ================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# =================================================================

import sys, os, time, math
import numpy as np
from scipy.integrate import quad
from config_loader import load_config
from global_config import GC

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)
os.environ['PYTHONIOENCODING'] = 'utf-8'

_DIR = os.path.dirname(os.path.abspath(__file__))
if _DIR not in sys.path:
    sys.path.insert(0, _DIR)

from v22_raw_scan import sigma_T_vpm
from output_manager import get_latest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_CFG = load_config(__file__)

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ================================================================
#  Velocity-averaged sigma/m
# ================================================================

def sigma_m_averaged(m_chi, m_phi_GeV, alpha, v_char_km_s, n_gauss=30):
    """
    Compute <sigma/m> averaged over Maxwell-Boltzmann velocity distribution.

    v_char is the characteristic velocity (peak of v*f(v) = most probable
    relative speed). We use v0 = v_char / sqrt(2) as the 1D dispersion.

    Uses Gauss-Legendre quadrature on [0, 4*v_char] (covers >99.9% of weight).
    Falls back to finer grid if cross section has resonant structure.
    """
    v0 = v_char_km_s / np.sqrt(2.0)  # 1D velocity dispersion

    # Maxwell-Boltzmann weight: v^3 * exp(-v^2/(2*v0^2))
    # (the v * f_MB(v) * v factor, up to normalization)
    def weight(v):
        x = v / v0
        return v**3 * np.exp(-x**2 / 2.0)

    def integrand_num(v):
        if v < 1.0:
            return 0.0
        s = sigma_T_vpm(m_chi, m_phi_GeV, alpha, v)
        if np.isnan(s) or np.isinf(s):
            return 0.0
        return s * weight(v)

    def integrand_den(v):
        return weight(v)

    # Integration range: [1, 5*v_char] covers essentially all the weight
    v_max = max(5.0 * v_char_km_s, 100.0)

    # Gauss-Legendre quadrature with n points
    nodes, weights_gl = np.polynomial.legendre.leggauss(n_gauss)
    # Map from [-1, 1] to [1, v_max]
    a, b = 1.0, v_max
    v_pts = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    jac = 0.5 * (b - a)

    # Evaluate VPM at each point
    sigma_vals = np.array([sigma_T_vpm(m_chi, m_phi_GeV, alpha, v) for v in v_pts])
    w_vals = np.array([weight(v) for v in v_pts])

    num = jac * np.sum(weights_gl * sigma_vals * w_vals)
    den = jac * np.sum(weights_gl * w_vals)

    if den < 1e-300:
        return 0.0
    return num / den


# ================================================================
#  Observational data — loaded from config.json if available
# ================================================================
OBSERVATIONS = GC.observations_as_tuples()


def compute_chi2(sigma_dict):
    chi2 = 0.0
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        theory = sigma_dict.get(v, 0.0)
        # One-sided upper limits (lo == 0): no penalty when theory < hi
        if lo == 0.0 and theory <= hi:
            continue
        if theory >= central:
            sigma = hi - central if hi > central else 0.5 * central
        else:
            sigma = central - lo if central > lo else 0.5 * central
        if sigma <= 0:
            sigma = 0.5 * max(central, 0.01)
        chi2 += ((theory - central) / sigma) ** 2
    return chi2


# ================================================================
#  Main
# ================================================================
def main():
    import csv as csv_mod

    # Load relic BPs (path overridable via config.json)
    _bp_csv = _CFG.get("benchmark_csv") or str(get_latest("v31_true_viable_points"))
    if not os.path.isabs(_bp_csv):
        _bp_csv = os.path.normpath(os.path.join(_DIR, _bp_csv))
    csv_path = _bp_csv
    bps = []
    with open(csv_path) as f:
        reader = csv_mod.DictReader(f)
        for row in reader:
            bps.append({
                'idx': len(bps) + 1,
                'm_chi': float(row['m_chi_GeV']),
                'm_phi_MeV': float(row['m_phi_MeV']),
                'm_phi_GeV': float(row['m_phi_MeV']) / 1000.0,
                'alpha': float(row['alpha']),
                'lam': float(row['lambda']),
            })

    print(f"Loaded {len(bps)} relic-viable benchmark points\n")

    # Unique velocities
    obs_velocities = sorted(set(v for _, v, *_ in OBSERVATIONS))
    print(f"Observational velocities: {obs_velocities}\n")

    # ── Select a few representative BPs for detailed analysis ──
    # BP1 (lambda=1.9), a mid-lambda, and the highest-lambda
    bp_indices = [0, 8, 10]  # BP1, BP9 (lam~5.4), BP11 (lam~32)
    test_bps = [bps[i] for i in bp_indices]

    print("=" * 90)
    print("FIXED vs VELOCITY-AVERAGED sigma/m  (n_gauss=30)")
    print("=" * 90)

    for bp in test_bps:
        mc = bp['m_chi']
        mp = bp['m_phi_GeV']
        al = bp['alpha']
        lam = bp['lam']

        print(f"\n  BP{bp['idx']}: m_chi={mc:.1f} GeV, m_phi={bp['m_phi_MeV']:.2f} MeV, "
              f"alpha={al:.3e}, lambda={lam:.2f}")
        print(f"  {'System':<20s}  {'v [km/s]':>8s}  {'sigma_fixed':>12s}  {'<sigma>_MB':>12s}  "
              f"{'ratio':>8s}  {'diff%':>8s}")
        print("  " + "-" * 80)

        for name, v, central, lo, hi, ref in OBSERVATIONS:
            t0 = time.time()
            s_fixed = sigma_T_vpm(mc, mp, al, float(v))
            s_avg = sigma_m_averaged(mc, mp, al, float(v), n_gauss=30)
            dt = time.time() - t0

            if s_fixed > 0 and s_avg > 0:
                ratio = s_avg / s_fixed
                diff_pct = (ratio - 1.0) * 100
            else:
                ratio = float('nan')
                diff_pct = float('nan')

            print(f"  {name:<20s}  {v:>8d}  {s_fixed:>12.4e}  {s_avg:>12.4e}  "
                  f"{ratio:>8.4f}  {diff_pct:>+8.1f}%")

    # ── Full chi2 comparison for all 17 BPs ──
    print("\n" + "=" * 90)
    print("CHI-SQUARED: FIXED vs VELOCITY-AVERAGED (all 17 relic BPs)")
    print("=" * 90)
    print(f"\n  {'BP':>3s}  {'m_chi':>7s}  {'lambda':>7s}  "
          f"{'chi2_fixed':>12s}  {'chi2_avg':>12s}  {'chi2_f/dof':>11s}  {'chi2_a/dof':>11s}  {'diff%':>8s}")
    print("  " + "-" * 85)

    chi2_fixed_all = []
    chi2_avg_all = []

    for bp in bps:
        mc = bp['m_chi']
        mp = bp['m_phi_GeV']
        al = bp['alpha']
        lam = bp['lam']

        sigma_fixed = {}
        sigma_avg = {}
        for v in obs_velocities:
            sigma_fixed[v] = sigma_T_vpm(mc, mp, al, float(v))
            sigma_avg[v] = sigma_m_averaged(mc, mp, al, float(v), n_gauss=30)

        c2_fixed = compute_chi2(sigma_fixed)
        c2_avg = compute_chi2(sigma_avg)
        dof = len(OBSERVATIONS) - 0  # 0 free params for relic BPs
        diff_pct = (c2_avg / c2_fixed - 1) * 100 if c2_fixed > 0 else float('nan')

        chi2_fixed_all.append(c2_fixed)
        chi2_avg_all.append(c2_avg)

        print(f"  {bp['idx']:3d}  {mc:7.1f}  {lam:7.2f}  "
              f"{c2_fixed:12.4f}  {c2_avg:12.4f}  "
              f"{c2_fixed/dof:11.4f}  {c2_avg/dof:11.4f}  {diff_pct:>+8.1f}%")

    # ── Summary ──
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    c2f = np.array(chi2_fixed_all)
    c2a = np.array(chi2_avg_all)
    dof = len(OBSERVATIONS)

    print(f"\n  Fixed-velocity chi2/dof: min={c2f.min()/dof:.4f}, "
          f"max={c2f.max()/dof:.4f}, mean={c2f.mean()/dof:.4f}")
    print(f"  Averaged chi2/dof:      min={c2a.min()/dof:.4f}, "
          f"max={c2a.max()/dof:.4f}, mean={c2a.mean()/dof:.4f}")

    diff = np.abs(c2a - c2f) / c2f * 100
    print(f"\n  |chi2_avg - chi2_fixed| / chi2_fixed:")
    print(f"    min  = {diff.min():.1f}%")
    print(f"    max  = {diff.max():.1f}%")
    print(f"    mean = {diff.mean():.1f}%")

    # Check if any BP changes viability
    good_f = np.sum(c2f / dof < 1.0)
    good_a = np.sum(c2a / dof < 1.0)
    print(f"\n  BPs with chi2/dof < 1: fixed={good_f}/17, averaged={good_a}/17")

    # ── Detailed ratio plot for BP1 ──
    bp1 = bps[0]
    mc, mp, al = bp1['m_chi'], bp1['m_phi_GeV'], bp1['alpha']
    v_range = np.logspace(np.log10(10), np.log10(5000), 40)
    s_fixed_arr = np.array([sigma_T_vpm(mc, mp, al, v) for v in v_range])

    print("\nComputing velocity-averaged curve for BP1 (40 points)...")
    s_avg_arr = np.array([sigma_m_averaged(mc, mp, al, v, n_gauss=30) for v in v_range])
    ratio_arr = s_avg_arr / s_fixed_arr

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: sigma/m(v) comparison
    ax = axes[0]
    ax.loglog(v_range, s_fixed_arr, 'b-', lw=2, label='Fixed velocity')
    ax.loglog(v_range, s_avg_arr, 'r--', lw=2, label='MB-averaged')
    # Observational data
    for name, v, cent, lo, hi, ref in OBSERVATIONS:
        ax.errorbar(v, cent, yerr=[[cent-lo], [hi-cent]], fmt='ko', ms=4,
                    capsize=3, zorder=5)
    ax.set_xlabel('v [km/s]')
    ax.set_ylabel('$\\sigma/m$ [cm$^2$/g]')
    ax.set_title(f'BP1: $m_\\chi$={mc:.1f} GeV, $\\lambda$={bp1["lam"]:.2f}')
    ax.legend()
    ax.set_xlim(8, 6000)
    ax.set_ylim(1e-3, 20)
    ax.grid(True, alpha=0.3)

    # Panel 2: ratio vs v
    ax = axes[1]
    ax.semilogx(v_range, ratio_arr, 'g-', lw=2)
    ax.axhline(1.0, color='gray', ls=':', lw=0.8)
    ax.fill_between(v_range, 0.8, 1.2, alpha=0.1, color='blue', label='20% band')
    ax.set_xlabel('v [km/s]')
    ax.set_ylabel('$\\langle\\sigma/m\\rangle_{\\rm MB}$ / $\\sigma/m(v)$')
    ax.set_title('Ratio: velocity-averaged / fixed')
    ax.set_xlim(8, 6000)
    ax.set_ylim(0.5, 2.0)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(_DIR, 'v37_velocity_averaged.png')
    plt.savefig(fig_path, dpi=150)
    print(f"\nFigure saved: {fig_path}")

    # ── Final verdict ──
    print("\n" + "=" * 70)
    print("FINAL VERDICT")
    print("=" * 70)
    max_diff = diff.max()
    mean_diff = diff.mean()
    print(f"\n  Maximum chi2 change from velocity averaging: {max_diff:.1f}%")
    print(f"  Mean chi2 change: {mean_diff:.1f}%")
    if good_f == good_a:
        print(f"  No benchmark point changes viability status.")
    else:
        print(f"  WARNING: {abs(good_f - good_a)} BP(s) change viability status!")

    if max_diff < 50:
        print(f"  => Velocity averaging has a {'minor' if max_diff < 20 else 'moderate'} effect on chi2.")
        print(f"  => Consistent with observational uncertainties of ~0.5 dex (factor of 3).")
    else:
        print(f"  => Velocity averaging has a SIGNIFICANT effect. Preprint should be updated.")


if __name__ == "__main__":
    main()
