#!/usr/bin/env python3
"""
cross_checks/sommerfeld_scattering_ratio.py
============================================
Sommerfeld / Scattering ratio  Φ(v) = S₀(v) / [σ_T(v)/m]

Physics:
  S₀(v) is the l=0 Sommerfeld enhancement — governs annihilation rate.
  σ_T(v)/m is the scattering cross section — governs halo dynamics.

  Both depend on the same Yukawa potential V(r) = -α e^{-m_φ r}/r,
  but weight partial waves differently:

    S₀ depends ONLY on l=0 (s-wave, even parity).
    σ_T^{Maj} sums ALL l with weights w_even=1, w_odd=3.
    σ_T^{Dir} sums ALL l with weights w=1.

  The ratio Φ(v) reveals the interplay between indirect detection
  (annihilation, proportional to S₀) and direct scattering constraints.

  Key prediction:
    Φ^{Maj}(v) ≠ Φ^{Dir}(v) because the scattering denominators differ.
    This means identical indirect-detection signals imply DIFFERENT
    self-interaction cross sections for Majorana vs Dirac dark matter.

  The ratio Φ^{Maj}/Φ^{Dir} = σ_T^{Dir}/σ_T^{Maj} at each velocity.

  Observational consequence:
    If CMB/Fermi-LAT constrain S₀ at some velocity, the implied σ_T/m
    differs for Majorana vs Dirac — providing a joint constraint that
    can distinguish particle nature.

Produces:
  - Console: Φ(v) summary table at reference velocities
  - output/sommerfeld_scattering_ratio.png  (3-panel figure)
  - output/sommerfeld_scattering_ratio.pdf
  - output/sommerfeld_scattering_ratio.csv
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

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm, vpm_phase_shift, GEV2_TO_CM2, GEV_IN_G, C_KM_S
from numba import jit

# Import Sommerfeld solver
sys.path.insert(0, _DIR)
from sommerfeld import sommerfeld_yukawa

# ── JIT warmup ──
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

cfg = load_config(__file__)
_bps = {bp['label']: bp for bp in GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1', 'MAP']))}
BP1 = _bps['BP1']
MAP = _bps['MAP']


# ══════════════════════════════════════════════════════════════
#  Dirac σ_T  (copied from predict_maj_vs_dir.py — weight=1 for all l)
# ══════════════════════════════════════════════════════════════
@jit(nopython=True, cache=True)
def sigma_T_dirac(m_chi, m_phi, alpha, v_km_s):
    """σ_T/m [cm²/g] for distinguishable Dirac fermion (all weights = 1)."""
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15:
        return 0.0

    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000

    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)

    sigma_sum = 0.0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        contrib = (2*l + 1) * math.sin(delta)**2
        sigma_sum += contrib
        if l > int(kappa) + 1 and sigma_sum > 0:
            if contrib / sigma_sum < 1e-3:
                break

    sigma_GeV2 = 4.0 * math.pi * sigma_sum / (k * k)
    sigma_cm2 = sigma_GeV2 * GEV2_TO_CM2
    return sigma_cm2 / (m_chi * GEV_IN_G)


# Warm up Dirac JIT
sigma_T_dirac(20.0, 10e-3, 1e-3, 100.0)


# ══════════════════════════════════════════════════════════════
#  Velocity grid
# ══════════════════════════════════════════════════════════════
N_PTS = 300
V_MIN, V_MAX = 3.0, 5000.0
V_GRID = np.logspace(np.log10(V_MIN), np.log10(V_MAX), N_PTS)

# Reference velocities for summary table
V_REFS = {
    'dSph': 12,
    'TBTF': 30,
    'LSB': 50,
    'MW': 220,
    'Cluster': 1200,
    'Bullet': 4700,
}


def compute_all(bp):
    """Compute S₀, σ_T^Maj, σ_T^Dir, and ratios Φ on the velocity grid."""
    m_chi = bp['m_chi_GeV']
    m_phi_GeV = bp['m_phi_MeV'] / 1000.0
    alpha = bp['alpha']

    S0 = np.zeros(N_PTS)
    sig_maj = np.zeros(N_PTS)
    sig_dir = np.zeros(N_PTS)

    for i, v in enumerate(V_GRID):
        S0[i] = sommerfeld_yukawa(m_chi, m_phi_GeV, alpha, v)
        sig_maj[i] = sigma_T_vpm(m_chi, m_phi_GeV, alpha, v)
        sig_dir[i] = sigma_T_dirac(m_chi, m_phi_GeV, alpha, v)

    # Φ = S₀ / (σ_T/m)  — units: [dimensionless] / [cm²/g] = [g/cm²]
    # Use masked arrays to avoid division by zero
    phi_maj = np.where(sig_maj > 0, S0 / sig_maj, np.nan)
    phi_dir = np.where(sig_dir > 0, S0 / sig_dir, np.nan)

    return S0, sig_maj, sig_dir, phi_maj, phi_dir


def main():
    t0 = time.time()
    print("=" * 70)
    print("  Sommerfeld / Scattering Ratio  Φ(v) = S₀(v) / [σ_T(v)/m]")
    print("=" * 70)

    results = {}
    for label, bp in [('BP1', BP1), ('MAP', MAP)]:
        lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000.0)
        print(f"\n  {label}:  m_χ={bp['m_chi_GeV']:.2f} GeV, "
              f"m_φ={bp['m_phi_MeV']:.2f} MeV, α={bp['alpha']:.4e}, λ={lam:.2f}")
        t1 = time.time()
        S0, sig_maj, sig_dir, phi_maj, phi_dir = compute_all(bp)
        dt = time.time() - t1
        print(f"    Computed {N_PTS} velocities in {dt:.1f}s")
        results[label] = (S0, sig_maj, sig_dir, phi_maj, phi_dir)

        # Summary table at reference velocities
        print(f"\n    {'v [km/s]':>10s}  {'S₀':>10s}  {'σ/m Maj':>10s}  "
              f"{'σ/m Dir':>10s}  {'Φ_Maj':>12s}  {'Φ_Dir':>12s}  {'Φ_M/Φ_D':>8s}")
        print("    " + "-" * 80)
        for name, v_ref in V_REFS.items():
            # Interpolate to reference velocity
            idx = np.argmin(np.abs(V_GRID - v_ref))
            s0_v = S0[idx]
            sm_v = sig_maj[idx]
            sd_v = sig_dir[idx]
            pm_v = phi_maj[idx]
            pd_v = phi_dir[idx]
            ratio = pm_v / pd_v if (not np.isnan(pm_v) and not np.isnan(pd_v) and pd_v != 0) else np.nan
            print(f"    {v_ref:>8d}  {s0_v:>10.4f}  {sm_v:>10.4f}  "
                  f"{sd_v:>10.4f}  {pm_v:>12.4e}  {pd_v:>12.4e}  {ratio:>8.4f}  ({name})")

    # ──────────────────────────────────────────────────────────
    #  3-panel figure
    # ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    colors = {'BP1': 'C0', 'MAP': 'C1'}
    ref_colors = {'dSph': 'green', 'TBTF': 'olive', 'LSB': 'teal',
                  'MW': 'orange', 'Cluster': 'red', 'Bullet': 'purple'}

    # Panel 1: S₀(v) and σ_T/m on same axis (dual y)
    ax1 = axes[0]
    ax1r = ax1.twinx()
    for label in ['BP1', 'MAP']:
        S0, sig_maj, sig_dir, _, _ = results[label]
        ax1.plot(V_GRID, S0, color=colors[label], ls='-', lw=2,
                 label=f'S₀ {label}')
        ax1r.plot(V_GRID, sig_maj, color=colors[label], ls='--', lw=1.5,
                  label=f'σ/m Maj {label}')
        ax1r.plot(V_GRID, sig_dir, color=colors[label], ls=':', lw=1.5,
                  label=f'σ/m Dir {label}')
    ax1.set_ylabel('S₀ (Sommerfeld)', fontsize=12)
    ax1r.set_ylabel('σ_T/m [cm²/g]', fontsize=12)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1r.set_yscale('log')
    ax1.set_title('Sommerfeld Enhancement & Scattering Cross Section', fontsize=13)
    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1r.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper right')
    # Reference velocity lines
    for name, v_ref in V_REFS.items():
        ax1.axvline(v_ref, color=ref_colors[name], alpha=0.3, ls='--', lw=0.8)

    # Panel 2: Φ_Maj and Φ_Dir
    ax2 = axes[1]
    for label in ['BP1', 'MAP']:
        _, _, _, phi_maj, phi_dir = results[label]
        ax2.plot(V_GRID, phi_maj, color=colors[label], ls='-', lw=2,
                 label=f'Φ_Maj {label}')
        ax2.plot(V_GRID, phi_dir, color=colors[label], ls='--', lw=1.5,
                 label=f'Φ_Dir {label}')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel('Φ = S₀ / (σ_T/m)  [g/cm²]', fontsize=12)
    ax2.set_title('Annihilation-to-Scattering Ratio Φ(v)', fontsize=13)
    ax2.legend(fontsize=9)
    for name, v_ref in V_REFS.items():
        ax2.axvline(v_ref, color=ref_colors[name], alpha=0.3, ls='--', lw=0.8)

    # Panel 3: Φ_Maj / Φ_Dir = σ_Dir / σ_Maj
    ax3 = axes[2]
    for label in ['BP1', 'MAP']:
        _, sig_maj, sig_dir, phi_maj, phi_dir = results[label]
        ratio = np.where((phi_dir > 0) & ~np.isnan(phi_dir),
                         phi_maj / phi_dir, np.nan)
        ax3.plot(V_GRID, ratio, color=colors[label], lw=2, label=label)
    ax3.axhline(1.0, color='grey', ls=':', lw=1)
    ax3.set_xscale('log')
    ax3.set_xlabel('v [km/s]', fontsize=12)
    ax3.set_ylabel('Φ_Maj / Φ_Dir = σ_Dir / σ_Maj', fontsize=12)
    ax3.set_title('Majorana vs Dirac Joint-Detection Discriminant', fontsize=13)
    ax3.legend(fontsize=10)
    for name, v_ref in V_REFS.items():
        ax3.axvline(v_ref, color=ref_colors[name], alpha=0.3, ls='--', lw=0.8)
        ax3.annotate(name, (v_ref, ax3.get_ylim()[0]),
                     fontsize=7, rotation=90, va='bottom', ha='right',
                     color=ref_colors[name], alpha=0.7)

    plt.tight_layout()
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'sommerfeld_scattering_ratio.{ext}'), dpi=200)
    plt.close(fig)
    print(f"\n  Figures saved to output/sommerfeld_scattering_ratio.{{png,pdf}}")

    # ──────────────────────────────────────────────────────────
    #  CSV output
    # ──────────────────────────────────────────────────────────
    csv_path = os.path.join(out_dir, 'sommerfeld_scattering_ratio.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['v_km_s',
                     'BP1_S0', 'BP1_sigmaT_Maj', 'BP1_sigmaT_Dir', 'BP1_Phi_Maj', 'BP1_Phi_Dir',
                     'MAP_S0', 'MAP_sigmaT_Maj', 'MAP_sigmaT_Dir', 'MAP_Phi_Maj', 'MAP_Phi_Dir'])
        for i in range(N_PTS):
            row = [f'{V_GRID[i]:.4f}']
            for label in ['BP1', 'MAP']:
                S0, sig_maj, sig_dir, phi_maj, phi_dir = results[label]
                row.extend([
                    f'{S0[i]:.8e}', f'{sig_maj[i]:.8e}', f'{sig_dir[i]:.8e}',
                    f'{phi_maj[i]:.8e}', f'{phi_dir[i]:.8e}'
                ])
            w.writerow(row)
    print(f"  CSV saved to {csv_path}")

    # ──────────────────────────────────────────────────────────
    #  Physics summary
    # ──────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  PHYSICS SUMMARY")
    print("=" * 70)

    for label in ['BP1', 'MAP']:
        S0, sig_maj, sig_dir, phi_maj, phi_dir = results[label]
        ratio = np.where((phi_dir > 0) & ~np.isnan(phi_dir),
                         phi_maj / phi_dir, np.nan)

        # Find where ratio deviates most from 1
        valid = ~np.isnan(ratio)
        if np.any(valid):
            max_dev_idx = np.argmax(np.abs(ratio[valid] - 1.0))
            v_max_dev = V_GRID[valid][max_dev_idx]
            max_ratio = ratio[valid][max_dev_idx]
        else:
            v_max_dev, max_ratio = np.nan, np.nan

        # Φ dynamic range
        phi_valid = phi_maj[~np.isnan(phi_maj)]
        if len(phi_valid) > 0:
            phi_range = phi_valid.max() / phi_valid.min()
        else:
            phi_range = np.nan

        print(f"\n  {label}:")
        print(f"    Φ_Maj dynamic range:  {phi_range:.1f}×  (across v = {V_MIN}–{V_MAX} km/s)")
        print(f"    Max Maj/Dir deviation: {max_ratio:.4f}  at v = {v_max_dev:.1f} km/s")
        print(f"      → {abs(max_ratio - 1)*100:.1f}% difference in indirect-detection implications")

        # Check at dSph velocities where CMB/Fermi constraints apply
        idx_dsph = np.argmin(np.abs(V_GRID - 12.0))
        phi_dsph_m = phi_maj[idx_dsph]
        phi_dsph_d = phi_dir[idx_dsph]
        print(f"    At dSph (v=12 km/s): Φ_Maj = {phi_dsph_m:.4e}, Φ_Dir = {phi_dsph_d:.4e}")
        print(f"      → Given same S₀, Majorana needs {phi_dsph_m/phi_dsph_d:.3f}× the σ/m of Dirac")

    dt_total = time.time() - t0
    print(f"\n  Total runtime: {dt_total:.1f}s")
    print("=" * 70)


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("✅ sommerfeld_scattering_ratio done!")
    except Exception:
        pass
