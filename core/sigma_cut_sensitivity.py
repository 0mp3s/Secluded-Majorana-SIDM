#!/usr/bin/env python3
"""
sigma_cut_sensitivity.py
========================
Post-process a raw viable CSV to find the optimal sigma/m(1000) upper bound.

Usage:
  # Analyse existing CSV (instant — cuts up to max(sigma_1000) in file):
  py core/sigma_cut_sensitivity.py

  # After running with relaxed cut for full picture:
  py core/v22_raw_scan_fast.py --sigma-cluster 1.0
  py core/sigma_cut_sensitivity.py

The script:
  1. Finds the latest all_viable_raw_v8_*.csv in data/archive/
  2. Sweeps sigma/m(1000) cut from 0.05 to max in fine steps
  3. Prints table of N_viable, N_BBN_safe, N_representative per cut
  4. Saves plot: sigma_cut_sensitivity.png
"""
import sys, os, glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from output_manager import timestamped_path


def find_latest_raw_csv():
    base = os.path.join(os.path.dirname(__file__), '..', 'data', 'archive')
    pattern = os.path.join(base, 'all_viable_raw_v8_*.csv')
    files = sorted(glob.glob(pattern))
    if not files:
        print("ERROR: No all_viable_raw_v8_*.csv found in data/archive/")
        sys.exit(1)
    return files[-1]


def analyse(csv_path):
    print(f"  Loading: {csv_path}")
    df = pd.read_csv(csv_path)
    n_total = len(df)
    max_sigma = df['sigma_m_1000'].max()
    print(f"  Total raw points: {n_total:,}")
    print(f"  sigma_1000 range: [{df['sigma_m_1000'].min():.6f}, {max_sigma:.6f}]")
    print(f"  m_chi range: [{df['m_chi_GeV'].min():.3f}, {df['m_chi_GeV'].max():.3f}] GeV")
    print()

    M_E_MEV = 0.511
    df['BBN_safe'] = df['m_phi_MeV'] > 2 * M_E_MEV

    # --- Sweep cuts ---
    cuts = sorted(set(
        list(np.arange(0.05, min(max_sigma + 0.05, 1.05), 0.05)) +
        [0.1, 0.2, 0.3, 0.35, 0.40, 0.42, 0.44, 0.46, 0.47, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0]
    ))
    cuts = [c for c in cuts if c <= max_sigma + 0.01]

    results = []
    print(f"  {'cut':>8}  {'N_viable':>10}  {'N_BBN':>8}  {'N_unique_mchi':>14}  {'hit_rate':>9}")
    print(f"  {'─'*8}  {'─'*10}  {'─'*8}  {'─'*14}  {'─'*9}")

    for cut in cuts:
        mask = df['sigma_m_1000'] < cut
        sub = df[mask]
        n_viable = len(sub)
        n_bbn = int(sub['BBN_safe'].sum())
        n_mchi = sub['m_chi_GeV'].nunique()
        hit_pct = n_viable / n_total * 100 if n_total > 0 else 0

        results.append({
            'sigma_cut': cut,
            'N_viable': n_viable,
            'N_BBN_safe': n_bbn,
            'N_unique_mchi': n_mchi,
            'hit_rate_pct': hit_pct
        })
        print(f"  {cut:8.3f}  {n_viable:10,}  {n_bbn:8,}  {n_mchi:14d}  {hit_pct:8.1f}%")

    print()
    res_df = pd.DataFrame(results)

    # --- Save results table ---
    out_csv = timestamped_path("sigma_cut_sensitivity")
    res_df.to_csv(out_csv, index=False)
    print(f"  Table saved: {out_csv}")

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax1 = plt.subplots(figsize=(10, 6))

        ax1.plot(res_df['sigma_cut'], res_df['N_viable'], 'b-o', ms=4, label='N viable (all)')
        ax1.plot(res_df['sigma_cut'], res_df['N_BBN_safe'], 'g-s', ms=4, label='N BBN-safe')
        ax1.set_xlabel(r'$\sigma/m(1000)$ upper cut [cm$^2$/g]', fontsize=12)
        ax1.set_ylabel('Number of viable points', fontsize=12, color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        ax1.legend(loc='upper left')
        ax1.grid(True, alpha=0.3)

        # Harvey+15 line
        ax1.axvline(x=0.47, color='red', ls='--', lw=1.5, alpha=0.7, label='Harvey+15 (0.47)')

        # Gradient (dN/dcut) on second axis
        if len(res_df) > 2:
            ax2 = ax1.twinx()
            grad = np.gradient(res_df['N_viable'], res_df['sigma_cut'])
            ax2.plot(res_df['sigma_cut'], grad, 'r-', alpha=0.5, lw=1.5, label='dN/d(cut)')
            ax2.set_ylabel('dN/d(cut)', fontsize=12, color='r')
            ax2.tick_params(axis='y', labelcolor='r')
            ax2.legend(loc='upper right')

        ax1.set_title(r'SIDM Viable Points vs $\sigma/m(1000\,\mathrm{km/s})$ Upper Bound', fontsize=13)
        fig.tight_layout()

        plot_path = timestamped_path("sigma_cut_sensitivity", ext=".png")
        fig.savefig(plot_path, dpi=150)
        print(f"  Plot saved: {plot_path}")
        plt.close()
    except ImportError:
        print("  (matplotlib not available — skipping plot)")

    # --- Key insight ---
    print()
    print("  === KEY INSIGHT ===")
    for i, row in res_df.iterrows():
        if i > 0:
            prev = res_df.iloc[i-1]
            delta = row['N_viable'] - prev['N_viable']
            if delta > 0:
                marginal = delta / (row['sigma_cut'] - prev['sigma_cut'])
                if marginal < 100 and prev['N_viable'] > 100:
                    print(f"  Diminishing returns above sigma_cut = {prev['sigma_cut']:.3f}")
                    print(f"    N jumps from {prev['N_viable']:,} to {row['N_viable']:,} "
                          f"(+{delta:,}, marginal = {marginal:.0f}/0.05)")
                    break

    # Where is the "knee"?
    if len(res_df) > 3:
        grad = np.gradient(res_df['N_viable'].values, res_df['sigma_cut'].values)
        grad2 = np.gradient(grad, res_df['sigma_cut'].values)
        knee_idx = np.argmax(np.abs(grad2[1:-1])) + 1
        knee_cut = res_df.iloc[knee_idx]['sigma_cut']
        print(f"  Steepest curvature (knee) at sigma_cut ≈ {knee_cut:.3f}")

    return res_df


def analyse_lower_bound(csv_path):
    """Sweep sigma/m(30 km/s) lower bound while keeping upper bounds fixed."""
    print()
    print("=" * 70)
    print("Lower Bound Sweep: sigma/m(30 km/s) >= X")
    print("=" * 70)
    print()

    df = pd.read_csv(csv_path)
    n_total = len(df)
    M_E_MEV = 0.511
    df['BBN_safe'] = df['m_phi_MeV'] > 2 * M_E_MEV
    max_s30 = df['sigma_m_30'].max()
    print(f"  sigma_30 range: [{df['sigma_m_30'].min():.4f}, {max_s30:.4f}]")
    print()

    lo_cuts = sorted(set([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]))

    results = []
    print(f"  {'lo_cut':>8}  {'N_viable':>10}  {'N_BBN':>8}  {'N_unique_mchi':>14}  {'m_chi_range':>20}")
    print(f"  {'─'*8}  {'─'*10}  {'─'*8}  {'─'*14}  {'─'*20}")

    for lo in lo_cuts:
        mask = df['sigma_m_30'] >= lo
        sub = df[mask]
        n_viable = len(sub)
        n_bbn = int(sub['BBN_safe'].sum())
        n_mchi = sub['m_chi_GeV'].nunique()
        if n_viable > 0:
            mrange = f"[{sub['m_chi_GeV'].min():.1f}, {sub['m_chi_GeV'].max():.1f}]"
        else:
            mrange = "—"

        results.append({
            'sigma_30_lo': lo,
            'N_viable': n_viable,
            'N_BBN_safe': n_bbn,
            'N_unique_mchi': n_mchi,
        })
        print(f"  {lo:8.2f}  {n_viable:10,}  {n_bbn:8,}  {n_mchi:14d}  {mrange:>20}")

    print()
    res_df = pd.DataFrame(results)

    out_csv = timestamped_path("sigma30_lo_sensitivity")
    res_df.to_csv(out_csv, index=False)
    print(f"  Table saved: {out_csv}")

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(res_df['sigma_30_lo'], res_df['N_viable'], 'b-o', ms=5, label='N viable (all)')
        ax.plot(res_df['sigma_30_lo'], res_df['N_BBN_safe'], 'g-s', ms=5, label='N BBN-safe')
        ax.set_xlabel(r'$\sigma/m(30\,\mathrm{km/s})$ lower bound [cm$^2$/g]', fontsize=12)
        ax.set_ylabel('Number of viable points', fontsize=12)
        ax.axvline(x=0.5, color='red', ls='--', lw=1.5, alpha=0.7, label='Elbert+15 (0.5)')
        ax.axvline(x=1.0, color='orange', ls='--', lw=1.5, alpha=0.7, label='Tulin & Yu 18 (1.0)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title(r'SIDM Viable Points vs $\sigma/m(30\,\mathrm{km/s})$ Lower Bound', fontsize=13)
        fig.tight_layout()

        plot_path = timestamped_path("sigma30_lo_sensitivity", ext=".png")
        fig.savefig(plot_path, dpi=150)
        print(f"  Plot saved: {plot_path}")
        plt.close()
    except ImportError:
        print("  (matplotlib not available — skipping plot)")

    return res_df


def main():
    print("=" * 70)
    print("SIDM Cut Sensitivity Analysis (Upper + Lower Bounds)")
    print("=" * 70)
    print()

    csv_path = sys.argv[1] if len(sys.argv) > 1 else find_latest_raw_csv()
    analyse(csv_path)
    analyse_lower_bound(csv_path)


if __name__ == "__main__":
    main()


if __name__ == '__main__':
    try:
        from tg_notify import notify
        notify("\u2705 sigma_cut_sensitivity done!")
    except Exception:
        pass
