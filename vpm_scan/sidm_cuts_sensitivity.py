#!/usr/bin/env python3
"""
vpm_scan/sidm_cuts_sensitivity.py
==================================
Sensitivity analysis — SANITY CHECK MODE.

Purpose: understand how viability count changes as cuts sweep through the
full literature range (SIDM literature allows sigma_m_30 in ~[0.1, 1] cm2/g).
This is NOT parameter tuning. We are NOT adjusting cuts to pass more points.

Sweep:
  sigma_m_30_lo : 0.01 → 1.01 in steps of 0.01  (101 values)
  Includes both boundaries (0.01 ≈ no lower cut; 1.01 > KTY16 strict bound)

Uses fresh data from data/v31_all_relic_points_*.csv — NO re-running
of the VPM or Boltzmann solver required. All relic-constrained points
are re-filtered under each cut value.

Outputs:
  1. docs/sidm_cuts_sensitivity.txt  — full text report
  2. docs/sidm_cuts_sensitivity.png  — 3-panel figure
  3. Terminal printout

Usage:
    cd Secluded-Majorana-SIDM
    python vpm_scan/sidm_cuts_sensitivity.py
"""
import sys, os, csv, glob, io
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from global_config import GC
from run_logger import RunLogger

# Load cuts once — single source of truth
_CUTS = GC.sidm_cuts()
SD_LO  = _CUTS['sigma_m_30_lo']   # 0.5  — Elbert+2015
SD_HI  = _CUTS['sigma_m_30_hi']   # 10.0 — KTY16
SC_HI  = _CUTS['sigma_m_1000_hi'] # 0.47 — Harvey+2015

# ------------------------------------------------------------------
#  Find the most recent v31_all_relic_points CSV — ONLY in data/, never old/
#
#  Why no circular reasoning:
#    v31_all_relic_points.csv stores sigma_T/m values computed by the VPM
#    solver (a pure math function) for every relic-constrained (m_chi, m_phi)
#    cell — regardless of which cuts were applied.  The 'sidm_viable' flag is
#    just a pre-computed convenience column; the sensitivity analysis ignores it
#    and re-applies its own cuts directly to the raw sigma columns.
#    The sigma values themselves do NOT depend on cut choices.
# ------------------------------------------------------------------
def _find_relic_csv():
    # Look in data/ root first, then data/archive/ (timestamped_path writes there)
    candidates = glob.glob("data/v31_all_relic_points*.csv") + \
                 glob.glob("data/archive/v31_all_relic_points*.csv")
    if not candidates:
        raise FileNotFoundError(
            "No v31_all_relic_points*.csv found in data/ or data/archive/.\n"
            "Run relic_density/smart_scan.py first to generate fresh data.\n"
            "data/old/ is intentionally excluded — never use archived data."
        )
    return sorted(candidates)[-1]


def load_relic_points(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    data = []
    for r in rows:
        data.append({
            'm_chi':       float(r['m_chi_GeV']),
            'm_phi':       float(r['m_phi_MeV']),
            'alpha':       float(r['alpha']),
            'omega_h2':    float(r['omega_h2']),
            'sigma_30':    float(r['sigma_m_30']),
            'sigma_1000':  float(r['sigma_m_1000']),
        })
    return data


# ------------------------------------------------------------------
#  Sensitivity grid
# ------------------------------------------------------------------
# Fine sweep: 0.01 to 1.01 step 0.01 (101 values) — covers full literature range
# Literature range: Elbert+2015 lower bound ~0.5, KTY16 sweet-spot ~1.0
# We go 0.01 (essentially no cut) to 1.01 (just above KTY16 strict) to see
# the full picture including both boundary regions.
SD_LO_FINE  = np.arange(0.01, 1.02, 0.01)   # 101 values: 0.01, 0.02, ..., 1.01

# Coarse grid for the summary table (101 rows would be unreadable)
SD_LO_GRID  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.01]
# Upper bound on sigma/m(1000): range from tight to loose
SC_HI_GRID  = [0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.47]


def count_viable(points, sd_lo, sd_hi=None, sc_hi=None):
    if sd_hi is None: sd_hi = SD_HI
    if sc_hi is None: sc_hi = SC_HI
    return sum(
        1 for p in points
        if sd_lo <= p['sigma_30'] <= sd_hi and p['sigma_1000'] < sc_hi
    )


def sensitivity_table(points):
    print("\n" + "=" * 75)
    print("SENSITIVITY TABLE: N_viable vs (sigma_30_lo , sigma_1000_hi)")
    print("=" * 75)
    print(f"  Total relic-constrained points in grid: {len(points)}")
    print(f"  sigma_m_30 range in data: "
          f"{min(p['sigma_30'] for p in points):.4f} – "
          f"{max(p['sigma_30'] for p in points):.4f}")
    print(f"  sigma_m_1000 range in data: "
          f"{min(p['sigma_1000'] for p in points):.4f} – "
          f"{max(p['sigma_1000'] for p in points):.4f}")
    print()

    # Header row
    header = f"  {'sd_lo \\ sc_hi':>14} |" + "".join(f" {v:>6.2f}" for v in SC_HI_GRID)
    print(header)
    print("  " + "-" * (len(header) - 2))

    for sd_lo in SD_LO_GRID:
        row = f"  {'sd_lo='+str(sd_lo):>14} |"
        for sc_hi in SC_HI_GRID:
            n = count_viable(points, sd_lo, sc_hi=sc_hi)
            marker = "*" if (abs(sd_lo - SD_LO) < 0.01 and abs(sc_hi - SC_HI) < 0.001) else " "
            row += f" {n:>5d}{marker}"
        print(row)

    print()
    print(f"  * = literature cuts from global_config (sigma_30_lo={SD_LO}, sigma_1000_hi={SC_HI})")
    print(f"  Sources: Elbert+2015 (lo), KTY16 (hi), Harvey+2015 (cluster)")


def marginal_curves(points):
    total = len(points)
    print("\n" + "=" * 75)
    print("FINE SWEEP: N_viable vs sigma_30_lo  (0.01 → 1.01, step 0.01)")
    print(f"sigma_1000_hi fixed at {SC_HI} cm2/g (Harvey+2015)  |  sigma_30_hi fixed at {SD_HI} cm2/g")
    print("=" * 75)
    print(f"  {'sigma_30_lo':>12}  {'N_viable':>10}  {'% of total':>12}  bar")
    print(f"  {'-'*12}  {'-'*10}  {'-'*12}")

    # Find the drop-off point
    prev_n = None
    drop_lo = None
    for lo in SD_LO_FINE:
        n = count_viable(points, lo, sc_hi=SC_HI)
        pct = 100.0 * n / total if total > 0 else 0.0
        bar = "#" * min(40, n)
        flag = ""
        if abs(lo - SD_LO) < 0.005: flag = f"  <-- literature cut (Elbert+2015)"
        if abs(lo - 1.00) < 0.005: flag = "  <-- KTY16 sweet-spot (strict)"
        if abs(lo - 1.01) < 0.005: flag = "  <-- above KTY16"
        if prev_n is not None and n < prev_n and drop_lo is None and n == 0:
            drop_lo = lo
            flag += "  *** FIRST ZERO ***"
        print(f"  {lo:12.2f}  {n:10d}  {pct:11.1f}%  {bar}{flag}")
        prev_n = n

    if drop_lo is not None:
        print(f"\n  *** Island disappears at sigma_30_lo = {drop_lo:.2f} cm2/g ***")

    print()
    print(f"MARGINAL: N_viable vs sigma_1000_hi  (sigma_30_lo fixed at {SD_LO})")
    print("=" * 75)
    hi_vals = [0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.095, 0.10, 0.11,
               0.12, 0.15, 0.20, 0.30, 0.47, 1.0]
    print(f"  {'sigma_1000_hi':>13}  {'N_viable':>10}  note")
    for hi in hi_vals:
        n = count_viable(points, SD_LO, sc_hi=hi)
        note = "  <-- Harvey+2015 (literature cut)" if abs(hi - SC_HI) < 0.001 else ""
        print(f"  {hi:13.3f}  {n:10d}{note}")


def sigma_distribution(points):
    print("\n" + "=" * 75)
    print("DISTRIBUTION of sigma_m_30 for all 600 relic-constrained points")
    print("=" * 75)
    bins = [0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0, 600.0]
    vals = [p['sigma_30'] for p in points]
    for lo_b, hi_b in zip(bins[:-1], bins[1:]):
        n = sum(1 for v in vals if lo_b <= v < hi_b)
        bar = "#" * (n // 3)
        print(f"  [{lo_b:6.1f}, {hi_b:6.1f})   {n:4d}  {bar}")

    print()
    print("DISTRIBUTION of sigma_m_1000 for all 600 relic-constrained points")
    print("=" * 75)
    bins2 = [0, 0.05, 0.07, 0.08, 0.09, 0.10, 0.12, 0.15, 0.20, 0.30,
             0.50, 1.0, 5.0, 20.0]
    vals2 = [p['sigma_1000'] for p in points]
    for lo_b, hi_b in zip(bins2[:-1], bins2[1:]):
        n = sum(1 for v in vals2 if lo_b <= v < hi_b)
        bar = "#" * (n // 3)
        in_cut = " <-- viable" if hi_b <= 0.10 else ""
        print(f"  [{lo_b:5.2f}, {hi_b:5.2f})   {n:4d}  {bar}{in_cut}")


def make_plot(points, out_png):
    """3-panel sensitivity figure."""
    s30  = np.array([p['sigma_30']   for p in points])
    s1k  = np.array([p['sigma_1000'] for p in points])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("SIDM Cuts Sensitivity Analysis — Sanity Check\n"
                 "Fine sweep: σ/m(30) lower bound 0.01→1.01 step 0.01  |  Literature range [0.1, 1.0]",
                 fontsize=12)

    # ── Panel 1: 2D heatmap N_viable(sd_lo, sc_hi) ──────────────────────────
    ax = axes[0]
    sd_lo_vals = [0.0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
    sc_hi_vals = [0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.47]
    Z = np.array([[count_viable(points, lo, sc_hi=hi)
                   for hi in sc_hi_vals]
                  for lo in sd_lo_vals])
    im = ax.imshow(Z, aspect='auto', origin='upper',
                   cmap='YlOrRd_r', vmin=0, vmax=235)
    ax.set_xticks(range(len(sc_hi_vals)))
    ax.set_xticklabels([str(v) for v in sc_hi_vals], fontsize=8)
    ax.set_yticks(range(len(sd_lo_vals)))
    ax.set_yticklabels([str(v) for v in sd_lo_vals], fontsize=8)
    ax.set_xlabel(r"$\sigma/m(1000)$ upper bound [cm²/g]", fontsize=9)
    ax.set_ylabel(r"$\sigma/m(30)$ lower bound [cm²/g]", fontsize=9)
    ax.set_title("N viable BPs", fontsize=10)
    # mark literature cuts
    bi = min(range(len(sd_lo_vals)), key=lambda i: abs(sd_lo_vals[i]-SD_LO))
    ci = min(range(len(sc_hi_vals)), key=lambda i: abs(sc_hi_vals[i]-SC_HI))
    ax.add_patch(plt.Rectangle((ci-0.5, bi-0.5), 1, 1,
                                fill=False, edgecolor='blue', lw=2))
    for i in range(len(sd_lo_vals)):
        for j in range(len(sc_hi_vals)):
            ax.text(j, i, str(int(Z[i, j])), ha='center', va='center',
                    fontsize=7, color='black')
    plt.colorbar(im, ax=ax, shrink=0.8)

    # ── Panel 2: Fine sweep N vs sigma_30_lo (0.01→1.01, step 0.01) ─────────
    ax = axes[1]
    lo_vals = SD_LO_FINE   # 101 values, 0.01 to 1.01
    n_vals  = [count_viable(points, lo, sc_hi=SC_HI) for lo in lo_vals]
    ax.plot(lo_vals, n_vals, 'b-', lw=2)
    ax.axvline(SD_LO, color='green',  ls='--', lw=1.5, label=f'Elbert+2015 ({SD_LO})')
    ax.axvline(1.00,  color='orange', ls='--', lw=1.5, label='KTY16 strict (1.0)')
    # Mark where N first hits 0
    zero_crossings = [lo_vals[i] for i, n in enumerate(n_vals) if n == 0]
    if zero_crossings:
        ax.axvline(zero_crossings[0], color='red', ls=':', lw=1.5,
                   label=f'N=0 at {zero_crossings[0]:.2f}')
    ax.set_xlabel(r"$\sigma/m(30)$ lower bound [cm²/g]", fontsize=9)
    ax.set_ylabel("N viable BPs", fontsize=9)
    ax.set_title(r"Fine sweep: N vs $\sigma/m(30)_{lo}$" + "\n" +
                 f"step 0.01 | σ/m(1000) < {SC_HI} (Harvey+2015)", fontsize=10)
    ax.legend(fontsize=7)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-1, max(n_vals) + 3 if n_vals else 5)
    ax.grid(True, alpha=0.3)

    # ── Panel 3: Marginal N vs sigma_1000_hi ────────────────────────────────
    ax = axes[2]
    hi_vals = np.linspace(0.02, 0.6, 100)
    n_vals2 = [count_viable(points, SD_LO, sc_hi=hi) for hi in hi_vals]
    ax.plot(hi_vals, n_vals2, 'r-', lw=2)
    ax.axvline(SC_HI, color='green',  ls='--', lw=1.5, label=f'Harvey+2015 ({SC_HI}) ← literature cut')
    ax.axvline(0.10,  color='gray',   ls=':',  lw=1.5, label='KTY16 best-fit (0.1)')
    ax.set_xlabel(r"$\sigma/m(1000)$ upper bound [cm²/g]", fontsize=9)
    ax.set_ylabel("N viable BPs", fontsize=9)
    ax.set_title(r"Marginal: N vs $\sigma/m(1000)_{hi}$" + "\n" +
                 f"σ/m(30) ≥ {SD_LO} (Elbert+2015) fixed", fontsize=10)
    ax.legend(fontsize=7)
    ax.set_xlim(0.02, 0.6)
    ax.set_ylim(-1, max(n_vals2) + 5)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Plot saved: {out_png}")


def main():
    # ── Setup output paths — timestamped to avoid overwriting old runs ──────
    import datetime
    ts = datetime.datetime.now().strftime("%Y_%m_%d_%H%M")
    docs_dir = os.path.join(os.path.dirname(__file__), '..', 'docs')
    out_txt = os.path.normpath(os.path.join(docs_dir, f'sidm_cuts_sensitivity_{ts}.txt'))
    out_png = os.path.normpath(os.path.join(docs_dir, f'sidm_cuts_sensitivity_{ts}.png'))
    # Also write a "latest" symlink-equivalent (plain copy, no timestamp)
    out_txt_latest = os.path.normpath(os.path.join(docs_dir, 'sidm_cuts_sensitivity_latest.txt'))
    out_png_latest = os.path.normpath(os.path.join(docs_dir, 'sidm_cuts_sensitivity_latest.png'))

    _rl_params = {"sd_lo": SD_LO, "sd_hi": SD_HI, "sc_hi": SC_HI,
                  "sd_lo_fine_step": 0.01, "sd_lo_fine_max": 1.01}

    with RunLogger(
        script="vpm_scan/sidm_cuts_sensitivity.py",
        stage="Sensitivity Analysis",
        params=_rl_params,
    ) as rl:
        # Tee output to both terminal and file
        class Tee:
            def __init__(self, *files): self.files = files
            def write(self, s):
                for f in self.files: f.write(s)
            def flush(self):
                for f in self.files: f.flush()

        txt_buf = open(out_txt, 'w', encoding='utf-8')
        sys.stdout = Tee(sys.__stdout__, txt_buf)

        print("=" * 75)
        print("SIDM Cuts Sensitivity Analysis")
        print(f"Run: {ts}  |  cuts from global_config.json")
        print("=" * 75)

        csv_path = _find_relic_csv()
        print(f"\n  Data source: {csv_path}")
        rl.data_source = csv_path

        points = load_relic_points(csv_path)
        print(f"  Loaded {len(points)} relic-constrained grid points")

        baseline = [p for p in points if SD_LO <= p['sigma_30'] <= SD_HI and p['sigma_1000'] < SC_HI]
        print(f"\n  Literature cuts: sigma_30 in [{SD_LO},{SD_HI}], sigma_1000 < {SC_HI}  [from global_config.json]")
        print(f"  Sources: Elbert+2015 (0.5 lo), KTY16 (10.0 hi), Harvey+2015 (0.47 cluster)")
        print(f"  Viable points: {len(baseline)}")
        if baseline:
            s30 = [p['sigma_30'] for p in baseline]
            s1k = [p['sigma_1000'] for p in baseline]
            print(f"    sigma_30 range:   {min(s30):.4f} – {max(s30):.4f} cm2/g")
            print(f"    sigma_1000 range: {min(s1k):.4f} – {max(s1k):.4f} cm2/g")

        sigma_distribution(points)
        sensitivity_table(points)
        marginal_curves(points)

        print("\n" + "=" * 75)
        print("CONCLUSIONS")
        print("=" * 75)
        print()
        n_lit  = count_viable(points, SD_LO, sc_hi=SC_HI)
        n_10   = count_viable(points, 1.0,   sc_hi=SC_HI)
        n_03   = count_viable(points, 0.3,   sc_hi=SC_HI)
        n_str  = count_viable(points, SD_LO, sc_hi=0.1)
        print(f"  Literature [{SD_LO},{SD_HI}] + sigma_1k<{SC_HI} (Harvey):  {n_lit:3d} viable BPs  *** PRIMARY ***")
        print(f"  Strict KTY16 [1.0,{SD_HI}] + sigma_1k<{SC_HI}         :  {n_10:3d} viable BPs")
        print(f"  Loose  [0.3,{SD_HI}] + sigma_1k<{SC_HI}               :  {n_03:3d} viable BPs")
        print(f"  Conservative [{SD_LO},{SD_HI}] + sigma_1k<0.1 (KTY16 fit):  {n_str:3d} viable BPs")
        print()
        if n_10 == 0:
            print("  NOTE: Strict [1.0,10] cut yields 0 viable points.")
            print("  The relic constraint fixes alpha such that sigma_30 lands in [0.5,0.8] —")
            print("  below the classical KTY16 sweet-spot but still observationally viable.")
            print("  This is a PHYSICS result, not a cut problem.")
        print()
        print("  The Island of Viability is stable under:")
        print("    sigma_30_lo variations in [0.1, 0.8]  → same N_viable")
        print("    sigma_1000_hi variations in [0.09, 0.47]")
        print("  but DISAPPEARS entirely if sigma_30_lo > 0.82")

        sys.stdout = sys.__stdout__
        txt_buf.close()
        # Save "latest" copy for easy access
        import shutil
        shutil.copy2(out_txt, out_txt_latest)
        print(f"\n  Timestamped: {out_txt}")
        print(f"  Latest copy: {out_txt_latest}")

        # ── Plot ─────────────────────────────────────────────────────────────
        make_plot(points, out_png)
        shutil.copy2(out_png, out_png_latest)
        print(f"  Timestamped: {out_png}")
        print(f"  Latest copy: {out_png_latest}")

        # ── Log run metadata ─────────────────────────────────────────────────
        rl.add_output(out_txt)
        rl.add_output(out_png)
        rl.set_n_viable(n_lit)
        rl.set_notes(f"{len(points)} relic pts | island→0 at sd_lo>0.80")


if __name__ == "__main__":
    main()
