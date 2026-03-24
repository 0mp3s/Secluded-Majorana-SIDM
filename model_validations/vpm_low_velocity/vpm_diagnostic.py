#!/usr/bin/env python3
"""
model_validations/vpm_low_velocity/vpm_diagnostic.py
====================================================
Compute σ_T/m at low velocities (v = 1–1000 km/s) for benchmark points.

Purpose: diagnostic table + plot for §7.2 (Fornax), §7.5 (UFDs),
and verification of s-wave plateau below first resonance.

For λ = α m_χ / m_φ ≈ 1.9 < π, σ_T saturates to a constant
at v ≲ v_med = (m_φ/m_χ)c ≈ 164 km/s (Tulin & Yu 2018, eq. 30).
"""
import sys, os, math
import numpy as np

# ---------- path bootstrap ----------
_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
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


def main():
    cfg = load_config(__file__)
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1']))
    velocities = cfg.get('velocities_km_s', [1, 3, 5, 8, 10, 12, 15, 20, 30, 50, 100, 300, 1000])

    print("=" * 90)
    print("  VPM Low-Velocity Diagnostic: σ_T/m [cm²/g] vs v [km/s]")
    print("=" * 90)

    # Header
    header = f"  {'v [km/s]':>10}"
    for bp in bps:
        header += f"  {bp['label']:>12}"
    print(header)
    print("  " + "─" * (10 + 14 * len(bps)))

    # Compute
    results = {bp['label']: [] for bp in bps}

    for v in velocities:
        row = f"  {v:>10.1f}"
        for bp in bps:
            m_chi = bp['m_chi_GeV']
            m_phi = bp['m_phi_MeV'] / 1000.0
            alpha = bp['alpha']
            sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, float(v))
            results[bp['label']].append(sigma_m)
            row += f"  {sigma_m:>12.4f}"
        print(row)

    # Physical scales
    print()
    print("  Physical scales:")
    print(f"  {'':>10}  {'λ':>12}  {'v_med [km/s]':>12}  {'v_Bohr [km/s]':>13}  {'Regime':>12}")
    print("  " + "─" * 65)
    for bp in bps:
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        lam = alpha * m_chi / m_phi
        v_med = (m_phi / m_chi) * 3e5  # c in km/s
        v_bohr = alpha * 3e5
        regime = "below 1st res" if lam < math.pi else f"resonance ({lam/math.pi:.1f}π)"
        print(f"  {bp['label']:>10}  {lam:>12.3f}  {v_med:>12.1f}  {v_bohr:>13.1f}  {regime:>12}")

    # Plateau analysis
    print()
    print("  Plateau check (σ/m variation for v ≤ 12 km/s):")
    for bp in bps:
        vals_low = []
        for i, v in enumerate(velocities):
            if v <= 12:
                vals_low.append(results[bp['label']][i])
        if len(vals_low) > 1:
            mean_val = np.mean(vals_low)
            max_dev = max(abs(x - mean_val) / mean_val for x in vals_low) * 100
            print(f"    {bp['label']}: σ/m(1–12) = {min(vals_low):.4f} – {max(vals_low):.4f} cm²/g "
                  f"(mean {mean_val:.4f}, max deviation {max_dev:.1f}%)")

    # Fornax-relevant values
    print()
    print("  Fornax-relevant (v = 12 km/s, σ_v = 11.7 km/s):")
    for bp in bps:
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        sig12 = sigma_T_vpm(m_chi, m_phi, alpha, 12.0)
        sig_vrel = sigma_T_vpm(m_chi, m_phi, alpha, 11.7 * math.sqrt(2))
        print(f"    {bp['label']}: σ/m(12) = {sig12:.4f}, σ/m(v_rel=16.5) = {sig_vrel:.4f} cm²/g")

    # UFD-relevant values
    print()
    print("  UFD-relevant (v = 3, 5, 8 km/s):")
    for bp in bps:
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        vals = []
        for v in [3, 5, 8]:
            s = sigma_T_vpm(m_chi, m_phi, alpha, float(v))
            vals.append(f"σ/m({v})={s:.3f}")
        print(f"    {bp['label']}: {', '.join(vals)} cm²/g")

    # ---- Plot ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    colors = {'BP1': 'steelblue', 'BP9': 'seagreen', 'MAP': 'firebrick'}
    v_fine = np.logspace(0, 3.2, 200)

    for bp in bps:
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        label = bp['label']
        sigma_fine = [sigma_T_vpm(m_chi, m_phi, alpha, float(v)) for v in v_fine]

        ax1.plot(v_fine, sigma_fine, color=colors.get(label, 'gray'),
                 lw=2, label=f"{label} (λ={alpha*m_chi/m_phi:.2f})")
        ax2.plot(v_fine, sigma_fine, color=colors.get(label, 'gray'), lw=2, label=label)

    # Mark key velocities
    for v_mark, lbl in [(3, 'Crater II'), (8, 'UFD'), (12, 'Fornax'), (30, 'dSph'), (1000, 'Cluster')]:
        ax1.axvline(v_mark, color='gray', ls=':', alpha=0.5)
        ax1.text(v_mark, ax1.get_ylim()[1] if ax1.get_ylim()[1] > 0 else 1.5,
                 lbl, rotation=90, va='top', ha='right', fontsize=7, alpha=0.7)

    ax1.set_xscale('log')
    ax1.set_xlabel('v [km/s]')
    ax1.set_ylabel('σ_T/m [cm²/g]')
    ax1.set_title('Full velocity range')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Zoom: low-v
    ax2.set_xlim(0, 50)
    ax2.set_xlabel('v [km/s]')
    ax2.set_ylabel('σ_T/m [cm²/g]')
    ax2.set_title('Low velocity zoom (UFD/Fornax regime)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Observational bands
    ax2.axhspan(0.2, 3.0, alpha=0.1, color='blue', label='Fornax obs (0.2–3)')
    ax2.axvspan(1, 12, alpha=0.05, color='green')

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'vpm_low_velocity.png')
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved: {fig_path}")
    plt.close(fig)

    # Save CSV
    csv_path = os.path.join(out_dir, 'vpm_low_velocity.csv')
    with open(csv_path, 'w') as f:
        cols = ['v_km_s'] + [bp['label'] + '_sigma_m' for bp in bps]
        f.write(','.join(cols) + '\n')
        for i, v in enumerate(velocities):
            vals = [f"{v:.1f}"]
            for bp in bps:
                vals.append(f"{results[bp['label']][i]:.6f}")
            f.write(','.join(vals) + '\n')
    print(f"  CSV saved: {csv_path}")


if __name__ == "__main__":
    main()
