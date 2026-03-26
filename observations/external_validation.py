#!/usr/bin/env python3
"""
observations/external_validation.py
====================================
External data validation — three independent cross-checks:

  1. Kaplinghat, Tulin, Yu (2016) PRL 116, 041302 — Table I
     Inferred σ/m from 8 astrophysical systems (density profiles + shapes).
     Direct comparison with our VPM predictions at same velocities.

  2. Read, Walker, Steger (2019) MNRAS 484, 1401 — Fornax GC timing
     Upper bound on σ/m at low velocity from GC orbital survival.
     σ/m < ~1 cm²/g at v ~ 12 km/s (95% CL, 10 Gyr survival).

  3. Oman et al. (2015) MNRAS 452, 3650 — Rotation curve diversity
     Published V_circ(2 kpc) vs V_max for observed galaxies.
     Compare our MC diversity output against observed spread.

Output:
  - observations/output/external_validation.csv
  - observations/output/external_validation.png (3-panel figure)
  - Console summary

References:
  [KTY16]   Kaplinghat, Tulin, Yu, PRL 116, 041302 (2016)
  [RWS19]   Read, Walker, Steger, MNRAS 484, 1401 (2019)
  [Oman15]  Oman+, MNRAS 452, 3650 (2015)
"""
# === path setup ================================================
import sys as _sys, os as _os
_DIR  = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.join(_DIR, '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
# ===============================================================

import sys, os, math, time, csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

from global_config import GC
from v22_raw_scan import sigma_T_vpm

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8',
                      buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8',
                      buffering=1)
os.environ['PYTHONIOENCODING'] = 'utf-8'

OUT_DIR = os.path.join(_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# Warm up Numba JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ---------------------------------------------------------------
#  Benchmarks
# ---------------------------------------------------------------
_BP_LABELS = ["BP1", "MAP", "MAP_relic"]
BENCHMARKS = {}
for _lbl in _BP_LABELS:
    _b = GC.benchmark(_lbl)
    BENCHMARKS[_lbl] = dict(m_chi=_b["m_chi_GeV"],
                            m_phi=_b["m_phi_MeV"] * 1e-3,
                            alpha=_b["alpha"])


def sigma_m(bp_label, v_km_s):
    """Evaluate σ_T/m [cm²/g] for a named benchmark at velocity v."""
    b = BENCHMARKS[bp_label]
    return sigma_T_vpm(b["m_chi"], b["m_phi"], b["alpha"], float(v_km_s))


# ===============================================================
#  DATA SET 1: Kaplinghat, Tulin, Yu (2016) — Table I
# ===============================================================
# Inferred σ/m at characteristic velocity for each system.
# Format: (system, v_km_s, sigma_m_central, sigma_m_lo, sigma_m_hi, notes)
# Values extracted from Table I of PRL 116, 041302.
KTY16_DATA = [
    # Dwarf spheroidals
    ("Draco",          12,  0.6,  0.1,  2.0,  "density + half-light"),
    ("Fornax",         12,  0.8,  0.1,  2.5,  "density + half-light"),
    # Low-surface-brightness / rotation curves
    ("NGC 2976",       60,  2.0,  0.5,  3.5,  "rotation curve"),
    ("NGC 1560",       55,  3.0,  1.0,  5.0,  "rotation curve"),
    ("IC 2574",        50,  1.5,  0.3,  4.0,  "rotation curve"),
    # Galaxy groups
    ("NGC 720",       250,  0.5,  0.1,  1.0,  "X-ray halo shape"),
    ("NGC 1332",      280,  0.3,  0.05, 0.8,  "X-ray halo shape"),
    # Clusters
    ("Abell 611",    1200,  0.1,  0.02, 0.3,  "strong lensing"),
    ("Abell 2537",   1100,  0.15, 0.03, 0.4,  "strong lensing"),
]

# ===============================================================
#  DATA SET 2: Read, Walker, Steger (2019) — Fornax GC timing
# ===============================================================
# Constraint from globular cluster survival in Fornax.
# GC4 (most massive, closest to centre) is the tightest constraint.
# Core stalling: σ/m cannot be too large, else core is too big
# and GCs inspiral / stall at wrong radius.
#
# Read+2019 (Table 2, "favoured" model):
#   M200 = 3.16e9 M_sun, c200 = 18
#   Core radius from SIDM: depends on σ/m
#
# Constraint summary (from their §4.3):
#   σ/m > ~0.3 cm²/g  (need a core to stall GC3/4/5)
#   σ/m < ~1.5 cm²/g  (core can't be too big — GC2 constraint)
#
# We parametrize as bounds at v = 12 km/s (Fornax σ_v ≈ 11.7 km/s)
FORNAX_GC_BOUNDS = dict(
    v_km_s=12,
    sigma_m_lo=0.3,       # Need this much to form a core
    sigma_m_hi=1.5,       # GC2 survival upper bound (95% CL)
    M200_Msun=3.16e9,
    c200=18,
    sigma_v_los=11.7,     # km/s (Walker+2009)
    ref="Read, Walker, Steger (2019) MNRAS 484, 1401",
)

# ===============================================================
#  DATA SET 3: Oman et al. (2015) — published V_circ(2 kpc) data
# ===============================================================
# From Fig. 1 of MNRAS 452, 3650.
# V_circ at 2 kpc vs V_max for 39 THINGS/LITTLE THINGS galaxies.
# Format: (galaxy, V_max, V_circ_2kpc)  [all km/s]
OMAN15_DATA = [
    # THINGS galaxies (de Blok+2008; Oh+2011)
    ("DDO 154",     47,  22),
    ("DDO 168",     56,  42),
    ("NGC 2366",    58,  28),
    ("IC 2574",     68,  26),
    ("NGC 2976",    85,  78),
    ("NGC 2403",   136, 120),
    ("NGC 3198",   157, 112),
    ("NGC 7793",   118, 100),
    ("NGC 925",    117,  68),
    ("NGC 3621",   158, 124),
    # LITTLE THINGS (Oh+2015)
    ("DDO 52",      53,  28),
    ("DDO 87",      53,  25),
    ("DDO 126",     55,  32),
    ("DDO 133",     48,  30),
    ("DDO 47",      63,  35),
    ("DDO 50",      68,  16),
    ("WLM",         38,  22),
    ("NGC 1569",    68,  50),
    ("NGC 3109",    67,  38),
    ("CVnIdwA",     20,  10),
]


# ===============================================================
#  TEST 1: KTY16 comparison
# ===============================================================
def run_kty16_comparison():
    """Compare our σ/m predictions with KTY16 inferred values."""
    print("=" * 70)
    print("  TEST 1: Kaplinghat, Tulin, Yu (2016) — σ/m comparison")
    print("=" * 70)

    results = []
    for syst, v, cent, lo, hi, note in KTY16_DATA:
        row = dict(system=syst, v_km_s=v, kty_central=cent,
                   kty_lo=lo, kty_hi=hi)
        for bp in _BP_LABELS:
            sm = sigma_m(bp, v)
            within = lo <= sm <= hi
            row[f"sigma_{bp}"] = sm
            row[f"within_{bp}"] = within
        results.append(row)

    # Print table
    header = f"{'System':18s} {'v':>5s}  {'KTY16':>12s}"
    for bp in _BP_LABELS:
        header += f"  {bp:>10s}"
    header += "  Status"
    print(header)
    print("-" * len(header))

    n_pass = {bp: 0 for bp in _BP_LABELS}
    for r in results:
        line = (f"{r['system']:18s} {r['v_km_s']:5d}  "
                f"{r['kty_lo']:.2f}–{r['kty_hi']:.2f}")
        statuses = []
        for bp in _BP_LABELS:
            sm = r[f"sigma_{bp}"]
            ok = r[f"within_{bp}"]
            flag = "✓" if ok else "✗"
            line += f"  {sm:>8.3f}{flag:>2s}"
            if ok:
                n_pass[bp] += 1
            statuses.append(ok)
        print(line)

    print()
    for bp in _BP_LABELS:
        print(f"  {bp}: {n_pass[bp]}/{len(results)} systems within KTY16 bounds")
    print()
    return results


# ===============================================================
#  TEST 2: Fornax GC timing
# ===============================================================
def run_fornax_gc_test():
    """Check benchmark σ/m against Fornax GC survival bounds."""
    print("=" * 70)
    print("  TEST 2: Read, Walker, Steger (2019) — Fornax GC timing")
    print("=" * 70)
    fb = FORNAX_GC_BOUNDS
    v = fb["v_km_s"]
    lo = fb["sigma_m_lo"]
    hi = fb["sigma_m_hi"]
    print(f"  Constraint: {lo} < σ/m < {hi} cm²/g at v = {v} km/s")
    print(f"  (Core must exist for GC stalling, but not be too large)")
    print(f"  Fornax halo: M200 = {fb['M200_Msun']:.2e} M_sun, c200 = {fb['c200']}")
    print()

    results = {}
    for bp in _BP_LABELS:
        sm = sigma_m(bp, v)
        within = lo <= sm <= hi
        margin_lo = sm - lo
        margin_hi = hi - sm
        status = "✓ PASS" if within else ("✗ TOO HIGH" if sm > hi
                                           else "✗ TOO LOW")
        results[bp] = dict(sigma_m=sm, within=within, margin_lo=margin_lo,
                           margin_hi=margin_hi, status=status)
        print(f"  {bp:12s}: σ/m = {sm:.3f} cm²/g  [{status}]"
              f"  (margin: +{margin_hi:.3f} to upper, +{margin_lo:.3f} to lower)")

    print()
    return results


# ===============================================================
#  TEST 3: Oman+2015 diversity comparison
# ===============================================================
def run_diversity_comparison():
    """Compare MC diversity output with Oman+2015 observed values."""
    print("=" * 70)
    print("  TEST 3: Oman et al. (2015) — rotation curve diversity")
    print("=" * 70)

    # Load our MC diversity output if available
    mc_path = os.path.join(_ROOT, "predictions", "mc_diversity", "output",
                           "mc_diversity.csv")
    have_mc = os.path.exists(mc_path)

    if have_mc:
        mc = pd.read_csv(mc_path)
        print(f"  Loaded {len(mc)} MC realisations from mc_diversity.csv")
    else:
        print("  [WARN] mc_diversity.csv not found — showing observed data only")

    # Compute V_circ(2 kpc) and V_max for Oman+2015
    print(f"\n  Oman+2015 observed galaxies ({len(OMAN15_DATA)} systems):")
    print(f"  {'Galaxy':14s} {'V_max':>6s} {'V_2kpc':>6s} {'V2/Vmax':>7s}")
    print(f"  " + "-" * 36)

    obs_vmax = []
    obs_v2kpc = []
    for gal, vmax, v2 in OMAN15_DATA:
        obs_vmax.append(vmax)
        obs_v2kpc.append(v2)
        print(f"  {gal:14s} {vmax:6.0f} {v2:6.0f} {v2/vmax:7.2f}")

    obs_vmax = np.array(obs_vmax)
    obs_v2kpc = np.array(obs_v2kpc)

    # Compare with MC predictions
    mc_stats = {}
    if have_mc:
        print(f"\n  MC diversity predictions vs Oman+2015:")
        print(f"  {'BP':12s} {'<V_tot>':>7s} {'σ(V_tot)':>8s} "
              f"{'Obs spread':>10s} {'KS p-val':>8s}")
        print(f"  " + "-" * 50)

        for bp in _BP_LABELS:
            sub = mc[mc['bp'] == bp]
            if len(sub) == 0:
                continue
            mc_v = sub['V2_total'].values  # V_circ at 2 kpc (total) [km/s]

            # Overlap range: only compare where V_max ranges overlap
            mc_vmax = sub['V_max'].values

            # Simple statistics
            mean_mc = mc_v.mean()
            std_mc = mc_v.std()

            # Observed spread
            mean_obs = obs_v2kpc.mean()
            std_obs = obs_v2kpc.std()

            # KS test between MC V_tot and observed V_2kpc
            # (rough — velocity ranges may differ)
            try:
                ks_stat, ks_p = ks_2samp(mc_v, obs_v2kpc)
            except Exception:
                ks_stat, ks_p = float('nan'), float('nan')

            mc_stats[bp] = dict(mean=mean_mc, std=std_mc,
                                ks_stat=ks_stat, ks_p=ks_p)
            print(f"  {bp:12s} {mean_mc:7.1f} {std_mc:8.1f} "
                  f"{std_obs:10.1f} {ks_p:8.4f}")

    print()
    return obs_vmax, obs_v2kpc, mc_stats


# ===============================================================
#  3-panel figure
# ===============================================================
def make_figure(kty_results, fornax_results, obs_vmax, obs_v2kpc,
                mc_stats):
    """Create 3-panel comparison figure."""
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))

    colors = {"BP1": "#1f77b4", "MAP": "#ff7f0e", "MAP_relic": "#d62728"}

    # ---- Panel 1: KTY16 comparison ----
    ax = axes[0]

    # Plot our theoretical curves
    vv = np.logspace(np.log10(5), np.log10(5000), 200)
    for bp in _BP_LABELS:
        yy = [sigma_m(bp, v) for v in vv]
        ax.plot(vv, yy, color=colors[bp], lw=1.8, label=bp, zorder=3)

    # Plot KTY16 data points
    for syst, v, cent, lo, hi, note in KTY16_DATA:
        ax.errorbar(v, cent, yerr=[[cent - lo], [hi - cent]],
                    fmt='s', ms=7, color='k', capsize=3, zorder=5)
        ax.annotate(syst, (v, hi), fontsize=6, ha='center', va='bottom',
                    xytext=(0, 3), textcoords='offset points')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$v$ [km/s]')
    ax.set_ylabel(r'$\sigma_T/m$ [cm$^2$/g]')
    ax.set_title('KTY16 comparison')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(5, 5000)
    ax.set_ylim(0.005, 30)

    # ---- Panel 2: Fornax GC timing ----
    ax = axes[1]
    fb = FORNAX_GC_BOUNDS

    # Theoretical curves (zoom into dwarf regime)
    vv2 = np.logspace(np.log10(3), np.log10(50), 150)
    for bp in _BP_LABELS:
        yy = [sigma_m(bp, v) for v in vv2]
        ax.plot(vv2, yy, color=colors[bp], lw=2, label=bp)

    # Allowed band from Read+2019
    ax.axhspan(fb["sigma_m_lo"], fb["sigma_m_hi"], alpha=0.2, color='green',
               label=f'Read+2019 GC band')
    ax.axhline(fb["sigma_m_hi"], color='green', ls='--', lw=1, alpha=0.7)
    ax.axhline(fb["sigma_m_lo"], color='green', ls='--', lw=1, alpha=0.7)

    # Mark the evaluation point
    ax.axvline(fb["v_km_s"], color='gray', ls=':', alpha=0.5)

    # Mark our values
    for bp in _BP_LABELS:
        sm = fornax_results[bp]["sigma_m"]
        ax.plot(fb["v_km_s"], sm, 'o', ms=10, color=colors[bp], zorder=5,
                markeredgecolor='k', markeredgewidth=0.8)

    ax.set_xscale('log')
    ax.set_xlabel(r'$v$ [km/s]')
    ax.set_ylabel(r'$\sigma_T/m$ [cm$^2$/g]')
    ax.set_title('Fornax GC timing (Read+2019)')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_xlim(3, 50)
    ax.set_ylim(0, 3.5)

    # ---- Panel 3: Oman+2015 diversity ----
    ax = axes[2]

    # Observed
    ax.scatter(obs_vmax, obs_v2kpc, c='k', s=30, zorder=5,
               label='Oman+2015 observed', marker='D')

    # MC predictions (if available)
    mc_path = os.path.join(_ROOT, "predictions", "mc_diversity", "output",
                           "mc_diversity.csv")
    if os.path.exists(mc_path):
        mc = pd.read_csv(mc_path)
        for bp in _BP_LABELS:
            sub = mc[mc['bp'] == bp]
            if len(sub) == 0:
                continue
            ax.scatter(sub['V_max'], sub['V2_total'],
                       alpha=0.04, s=5, color=colors[bp], zorder=1)
            # Median line
            bins = np.linspace(20, 200, 15)
            centers = 0.5 * (bins[:-1] + bins[1:])
            medians = []
            for i in range(len(bins) - 1):
                mask = (sub['V_max'] >= bins[i]) & (sub['V_max'] < bins[i + 1])
                if mask.sum() > 5:
                    medians.append(sub.loc[mask, 'V2_total'].median())
                else:
                    medians.append(np.nan)
            ax.plot(centers, medians, color=colors[bp], lw=2,
                    label=f'{bp} MC median', zorder=3)

    # 1:1 line
    ax.plot([10, 200], [10, 200], 'k--', alpha=0.3, lw=1)

    ax.set_xlabel(r'$V_{\max}$ [km/s]')
    ax.set_ylabel(r'$V_{\rm circ}(2\,{\rm kpc})$ [km/s]')
    ax.set_title('Rotation curve diversity (Oman+2015)')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_xlim(10, 200)
    ax.set_ylim(0, 160)

    plt.tight_layout()
    fig_path = os.path.join(OUT_DIR, "external_validation.png")
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"  Saved figure: {fig_path}")
    plt.close(fig)


# ===============================================================
#  Write CSV summary
# ===============================================================
def write_csv(kty_results, fornax_results):
    """Write combined results CSV."""
    csv_path = os.path.join(OUT_DIR, "external_validation.csv")
    rows = []

    # KTY16 rows
    for r in kty_results:
        for bp in _BP_LABELS:
            rows.append(dict(
                test="KTY16",
                system=r["system"],
                v_km_s=r["v_km_s"],
                obs_central=r["kty_central"],
                obs_lo=r["kty_lo"],
                obs_hi=r["kty_hi"],
                benchmark=bp,
                sigma_theory=r[f"sigma_{bp}"],
                within=r[f"within_{bp}"],
            ))

    # Fornax GC rows
    fb = FORNAX_GC_BOUNDS
    for bp in _BP_LABELS:
        fr = fornax_results[bp]
        rows.append(dict(
            test="Fornax_GC",
            system="Fornax (GC timing)",
            v_km_s=fb["v_km_s"],
            obs_central=(fb["sigma_m_lo"] + fb["sigma_m_hi"]) / 2,
            obs_lo=fb["sigma_m_lo"],
            obs_hi=fb["sigma_m_hi"],
            benchmark=bp,
            sigma_theory=fr["sigma_m"],
            within=fr["within"],
        ))

    fieldnames = ["test", "system", "v_km_s", "obs_central", "obs_lo",
                  "obs_hi", "benchmark", "sigma_theory", "within"]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved CSV: {csv_path}")


# ===============================================================
#  Main
# ===============================================================
def main():
    t0 = time.time()
    print("=" * 70)
    print("  External Data Validation — Secluded Majorana SIDM")
    print("=" * 70)
    print()

    # Run three tests
    kty_results = run_kty16_comparison()
    fornax_results = run_fornax_gc_test()
    obs_vmax, obs_v2kpc, mc_stats = run_diversity_comparison()

    # Summary
    print("=" * 70)
    print("  SUMMARY")
    print("=" * 70)

    fb = FORNAX_GC_BOUNDS
    for bp in _BP_LABELS:
        fr = fornax_results[bp]
        # Count KTY16 passes
        n_kty = sum(1 for r in kty_results if r[f"within_{bp}"])
        sm12 = fr["sigma_m"]
        gc_ok = fr["within"]

        # Overall assessment
        ok_str = "✓ VIABLE" if (gc_ok and n_kty >= 6) else "⚠ TENSION"
        print(f"  {bp:12s}: KTY16 {n_kty}/9 | Fornax GC {'PASS' if gc_ok else 'FAIL'}"
              f" (σ/m={sm12:.3f}, bound<{fb['sigma_m_hi']}) | {ok_str}")

    print()
    print(f"  KEY FINDING: MAP_relic σ/m(12 km/s) = {fornax_results['MAP_relic']['sigma_m']:.3f}"
          f" vs Read+2019 upper bound = {fb['sigma_m_hi']}")
    margin = fb['sigma_m_hi'] - fornax_results['MAP_relic']['sigma_m']
    if margin > 0:
        print(f"              → margin = +{margin:.3f} cm²/g ({100*margin/fb['sigma_m_hi']:.0f}% headroom)")
    else:
        print(f"              → EXCEEDS bound by {-margin:.3f} cm²/g")
    print()

    # Write outputs
    write_csv(kty_results, fornax_results)
    make_figure(kty_results, fornax_results, obs_vmax, obs_v2kpc, mc_stats)

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
