#!/usr/bin/env python3
"""
observations/opusB_observational_comparison.py
================================================
FIXED COPY of observational_comparison.py — by Opus B.

Changes vs original:
  1. BUG FIX (Harvey+15 velocity): changed from v=1500 to v=1000 km/s.
     Harvey+15 (MNRAS 449, 3393) reports σ/m < 0.47 cm²/g from 72 cluster
     mergers.  The effective pairwise velocity is ~1000 km/s per Robertson+17
     and most recent SIDM literature.  KTY16 mapped it to ~1500, but this
     overestimates the collision speed.

Everything else is identical to the original.
"""
# === path setup (auto-generated) ================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# =================================================================

import sys, os, math, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

# ==============================================================
#  REAL OBSERVATIONAL DATA — from global_config.json
# ==============================================================
_CFG = load_config(__file__)

OBSERVATIONS = GC.observations_as_tuples()

# ==============================================================
#  Our benchmark points — from global_config.json
# ==============================================================

_BP_LABELS = _CFG.get("benchmark_labels", ["BP1", "MAP"])
BENCHMARKS = []
for _lbl in _BP_LABELS:
    _b = GC.benchmark(_lbl)
    BENCHMARKS.append((_lbl, _b["m_chi_GeV"], _b["m_phi_MeV"] * 1e-3, _b["alpha"]))

VELOCITIES = np.logspace(np.log10(5), np.log10(5000), 100)

# ==============================================================
#  Compute theoretical curves
# ==============================================================

def compute_curve(m_chi, m_phi, alpha, velocities):
    """Compute σ/m at each velocity."""
    sigma_m = []
    for v in velocities:
        s = sigma_T_vpm(m_chi, m_phi, alpha, float(v))
        sigma_m.append(s)
    return np.array(sigma_m)


def main():
    t0 = time.time()

    print("=" * 80)
    print("  V10 — v33 Observational Comparison")
    print("  Theoretical σ/m(v) vs Real Astrophysical Data")
    print("=" * 80)
    print()

    # --- Warm up Numba ---
    print("  Warming up Numba JIT...", end=" ", flush=True)
    _ = sigma_T_vpm(10.0, 0.01, 1e-3, 30.0)
    print("done.\n")

    # --- Compute theoretical curves ---
    curves = {}
    for name, mc, mp, alpha in BENCHMARKS:
        print(f"  Computing {name}: m_chi={mc}, m_phi={mp*1e3:.2f} MeV, alpha={alpha:.4e}")
        curve = compute_curve(mc, mp, alpha, VELOCITIES)
        curves[name] = curve
        # Report key velocities
        s30 = sigma_T_vpm(mc, mp, alpha, 30.0)
        s100 = sigma_T_vpm(mc, mp, alpha, 100.0)
        s1000 = sigma_T_vpm(mc, mp, alpha, 1000.0)
        print(f"    σ/m(30)={s30:.4f}, σ/m(100)={s100:.4f}, σ/m(1000)={s1000:.4f} cm²/g")
    print()

    # --- Print observational data ---
    print("  OBSERVATIONAL DATA POINTS:")
    print(f"  {'System':<25} {'v [km/s]':>10} {'σ/m':>8} {'lo':>8} {'hi':>8} {'Ref':<12}")
    print("  " + "-" * 75)
    for name, v, sm, lo, hi, ref in OBSERVATIONS:
        print(f"  {name:<25} {v:10.0f} {sm:8.2f} {lo:8.2f} {hi:8.2f} {ref:<12}")
    print()

    # --- Compatibility check (all benchmarks) ---
    for bp_name, bp_mc, bp_mp, bp_alpha in BENCHMARKS:
        print(f"  COMPATIBILITY CHECK ({bp_name} vs observations):")
        print("  " + "-" * 60)
        n_compat = 0
        n_total = len(OBSERVATIONS)
        for name, v, sm_obs, lo, hi, ref in OBSERVATIONS:
            sm_theory = sigma_T_vpm(bp_mc, bp_mp, bp_alpha, float(v))
            in_range = lo <= sm_theory <= hi
            if in_range:
                n_compat += 1
                status = "COMPATIBLE"
            else:
                if sm_theory < lo:
                    status = f"BELOW (theory={sm_theory:.3f} < {lo:.2f})"
                else:
                    status = f"ABOVE (theory={sm_theory:.3f} > {hi:.2f})"
            print(f"  {name:<25} v={v:5.0f}  theory={sm_theory:.4f}  obs=[{lo:.2f},{hi:.2f}]  {status}")

        print()
        print(f"  {bp_name} compatible with {n_compat}/{n_total} observations")
        print()

    # --- Generate Figure 3 ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    # Plot theoretical curves (dynamic from config benchmarks)
    _palette = ['#2563EB', '#DC2626', '#16A34A', '#F59E0B', '#8B5CF6', '#6B7280']
    for i, (bname, mc, mp, alpha) in enumerate(BENCHMARKS):
        color = _palette[i % len(_palette)]
        label = rf'{bname}: $m_\chi$={mc:.1f} GeV, $m_\phi$={mp*1e3:.1f} MeV'
        ax.plot(VELOCITIES, curves[bname], color=color,
                linewidth=2.5, label=label, zorder=5)

    # Plot observational data with error bars
    obs_colors = {
        'KTY16': '#F59E0B',
        'KKPY17': '#8B5CF6',
        'Randall+08': '#EF4444',
        'Harvey+15': '#EF4444',
        'Elbert+15': '#10B981',
    }
    obs_markers = {
        'KTY16': 'o',
        'KKPY17': 's',
        'Randall+08': 'v',
        'Harvey+15': '^',
        'Elbert+15': 'D',
    }
    plotted_refs = set()
    for name, v, sm, lo, hi, ref in OBSERVATIONS:
        yerr_lo = sm - lo
        yerr_hi = hi - sm
        color = obs_colors.get(ref, '#6B7280')
        marker = obs_markers.get(ref, 'o')
        label = ref if ref not in plotted_refs else None
        plotted_refs.add(ref)
        ax.errorbar(v, sm, yerr=[[yerr_lo], [yerr_hi]],
                    fmt=marker, color=color, markersize=8,
                    capsize=4, capthick=1.5, linewidth=1.5,
                    label=label, zorder=10, markeredgecolor='black',
                    markeredgewidth=0.5)

    # Shaded bands for key constraints
    ax.axhspan(0.5, 10.0, xmin=0, xmax=0.15, alpha=0.08, color='green',
               label=r'Dwarf requirement ($\sigma/m \geq 0.5$)')
    ax.axhspan(0, 0.1, xmin=0.75, xmax=1.0, alpha=0.08, color='red',
               label=r'Cluster bound ($\sigma/m < 0.1$)')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(5, 5000)
    ax.set_ylim(1e-3, 50)
    ax.set_xlabel('Relative velocity $v$ [km/s]', fontsize=14)
    ax.set_ylabel(r'$\sigma_T/m_\chi$ [cm$^2$/g]', fontsize=14)
    ax.set_title('Theoretical Prediction vs Observational Data\n'
                 'Secluded Majorana SIDM — Scalar Mediator', fontsize=14)
    ax.legend(fontsize=8.5, loc='upper right', ncol=2, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    # Reference lines
    ax.axhline(y=0.1, color='red', linestyle='--', alpha=0.4, linewidth=0.8)
    ax.axhline(y=0.5, color='green', linestyle='--', alpha=0.4, linewidth=0.8)

    plt.tight_layout()
    outpath = os.path.join(_DIR, 'v33_observational_comparison.png')
    fig.savefig(outpath, dpi=200)
    print(f"  Saved figure → {outpath}")
    plt.close()

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")
    print("=" * 80)


if __name__ == "__main__":
    main()
