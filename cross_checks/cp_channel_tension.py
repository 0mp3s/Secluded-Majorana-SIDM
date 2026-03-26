#!/usr/bin/env python3
"""
cross_checks/cp_channel_tension.py
=====================================
CP-channel decomposition at tension velocities — is CP-odd peaked at 50-60 km/s?

Physics:
  For a Majorana fermion scattering via Yukawa mediator, the transport
  cross section weights even-l (singlet, CP-even) with w=1 and
  odd-l (triplet, CP-odd) with w=3:

    σ_T^{Maj} = (4π/k²) Σ_l (w_l)(2l+1) sin²δ_l

  where w_even=1, w_odd=3, and (2l+1) is the standard degeneracy factor.

  We decompose σ_T into:
    σ_even(v) = (4π/k²) Σ_{l even} (2l+1) sin²δ_l         [weight=1]
    σ_odd(v)  = (4π/k²) Σ_{l odd}  (2l+1) sin²δ_l          [weight=1, before triplet factor]

  so that σ_T^{Maj} = σ_even + 3·σ_odd
  and     σ_T^{Dir} = σ_even +   σ_odd

  The "triplet enhancement" at each velocity is:
    R_triplet(v) = 3·σ_odd / σ_T^{Maj}

  Key diagnostic:
    If R_triplet peaks at v ≈ 50-60 km/s where MAP_relic has tension,
    the tension is a direct signature of CP-odd physics —
    a smoking gun for the Majorana nature of dark matter.

    If R_triplet is flat, the tension is from generic VPM resonance
    structure (affects both Majorana and Dirac equally).

Produces:
  - Console: channel breakdown at observational velocities
  - output/cp_channel_tension.png  (3-panel figure)
  - output/cp_channel_tension.pdf
  - output/cp_channel_tension.csv
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
from v22_raw_scan import vpm_phase_shift, GEV2_TO_CM2, GEV_IN_G, C_KM_S
from numba import jit

cfg = load_config(__file__)
_bps = {bp['label']: bp for bp in
        GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1', 'MAP', 'MAP_relic']))}

OBSERVATIONS = GC.observations_as_tuples()


# ══════════════════════════════════════════════════════════════
#  CP-channel separated cross sections
# ══════════════════════════════════════════════════════════════
@jit(nopython=True, cache=True)
def sigma_channels(m_chi, m_phi, alpha, v_km_s):
    """
    Returns (sigma_even, sigma_odd) in cm²/g, BEFORE triplet weighting.

    σ_T^Maj = sigma_even + 3 * sigma_odd
    σ_T^Dir = sigma_even + sigma_odd
    """
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15:
        return 0.0, 0.0

    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000

    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)

    sum_even = 0.0
    sum_odd = 0.0
    total = 0.0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        contrib = (2*l + 1) * math.sin(delta)**2
        if l % 2 == 0:
            sum_even += contrib
        else:
            sum_odd += contrib
        total += contrib
        if l > int(kappa) + 1 and total > 0:
            if contrib / total < 1e-3:
                break

    prefactor = 4.0 * math.pi / (k * k) * GEV2_TO_CM2 / (m_chi * GEV_IN_G)
    return sum_even * prefactor, sum_odd * prefactor


# JIT warmup
sigma_channels(20.0, 10e-3, 1e-3, 100.0)


# ══════════════════════════════════════════════════════════════
#  Velocity grids
# ══════════════════════════════════════════════════════════════
N_PTS = 300
V_GRID = np.logspace(np.log10(3.0), np.log10(5000.0), N_PTS)


def main():
    t0 = time.time()
    hdr = "=" * 80
    print(hdr)
    print("  CP-Channel Decomposition at Tension Velocities")
    print("  σ_even (singlet) vs σ_odd (triplet) — is CP-odd peaked at 50-60 km/s?")
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

        # ── Dense velocity scan ──
        t1 = time.time()
        sig_even = np.empty(N_PTS)
        sig_odd = np.empty(N_PTS)
        for i, v in enumerate(V_GRID):
            se, so = sigma_channels(mc, mp_GeV, al, float(v))
            sig_even[i] = se
            sig_odd[i] = so

        sig_maj = sig_even + 3.0 * sig_odd
        sig_dir = sig_even + sig_odd

        # Triplet fraction R = 3·σ_odd / σ_Maj
        with np.errstate(divide='ignore', invalid='ignore'):
            R_triplet = np.where(sig_maj > 0, 3.0 * sig_odd / sig_maj, 0.0)
            ratio_maj_dir = np.where(sig_dir > 0, sig_maj / sig_dir, 1.0)

        dt_scan = time.time() - t1
        print(f"    Dense scan: {dt_scan:.1f}s  ({N_PTS} velocities)")

        # ── Observations table ──
        print(f"\n    {'System':<24s} {'v₀':>5s}  {'σ_even':>8s}  {'σ_odd':>8s}  "
              f"{'σ_Maj':>8s}  {'σ_Dir':>8s}  {'R_trip':>6s}  {'Maj/Dir':>7s}")
        print("    " + "-" * 85)

        obs_data = []
        for name, v, central, lo, hi, ref in OBSERVATIONS:
            se, so = sigma_channels(mc, mp_GeV, al, float(v))
            sm = se + 3.0 * so
            sd = se + so
            rt = 3.0 * so / sm if sm > 0 else 0.0
            md = sm / sd if sd > 0 else 1.0
            obs_data.append(dict(name=name, v=v, sig_even=se, sig_odd=so,
                                 sig_maj=sm, sig_dir=sd, R_triplet=rt,
                                 maj_dir_ratio=md))
            flag = "◄◄" if 40 <= v <= 80 else ""
            print(f"    {name:<24s} {v:>5.0f}  {se:>8.4f}  {so:>8.4f}  "
                  f"{sm:>8.4f}  {sd:>8.4f}  {rt:>6.1%}  {md:>7.3f}  {flag}")

        # Peak triplet fraction
        i_peak = np.argmax(R_triplet)
        v_peak = V_GRID[i_peak]
        R_peak = R_triplet[i_peak]
        # Average R at 40-80 km/s
        mask_tension = (V_GRID >= 40) & (V_GRID <= 80)
        R_mean_tension = np.mean(R_triplet[mask_tension]) if np.any(mask_tension) else 0
        R_mean_all = np.mean(R_triplet[sig_maj > 0])

        print(f"\n    Peak R_triplet: {R_peak:.1%} at v={v_peak:.1f} km/s")
        print(f"    Mean R_triplet (40-80 km/s): {R_mean_tension:.1%}")
        print(f"    Mean R_triplet (all v):      {R_mean_all:.1%}")
        print(f"    Is CP-odd enhanced at tension velocities? "
              f"{'YES — R(40-80) > R(all)' if R_mean_tension > R_mean_all + 0.02 else 'NO — CP-odd fraction is typical'}")

        all_results[label] = dict(
            sig_even=sig_even, sig_odd=sig_odd, sig_maj=sig_maj,
            sig_dir=sig_dir, R_triplet=R_triplet,
            ratio_maj_dir=ratio_maj_dir, obs_data=obs_data,
            v_peak=v_peak, R_peak=R_peak,
            R_tension=R_mean_tension, R_all=R_mean_all)

    # ──────────────────────────────────────────────────────────
    #  Figure: 3-panel per benchmark
    # ──────────────────────────────────────────────────────────
    n_bp = len(all_results)
    fig, axes = plt.subplots(n_bp, 3, figsize=(16, 4.5 * n_bp))
    if n_bp == 1:
        axes = axes.reshape(1, -1)

    for idx, (label, res) in enumerate(all_results.items()):
        # Panel 1: Channel cross sections
        ax1 = axes[idx, 0]
        ax1.loglog(V_GRID, res['sig_even'], 'b-', lw=1.5, label='σ_even (singlet)')
        ax1.loglog(V_GRID, 3.0 * res['sig_odd'], 'r-', lw=1.5, label='3·σ_odd (triplet)')
        ax1.loglog(V_GRID, res['sig_maj'], 'k-', lw=2, alpha=0.7, label='σ_T^Maj (total)')
        ax1.axvspan(40, 80, alpha=0.15, color='orange', label='tension region')
        # Mark observations
        for od in res['obs_data']:
            ax1.axvline(od['v'], color='gray', lw=0.5, alpha=0.3)
        ax1.set_xlabel('v [km/s]', fontsize=11)
        ax1.set_ylabel('σ/m [cm²/g]', fontsize=11)
        ax1.set_title(f'{label} — CP-channel decomposition', fontsize=12)
        ax1.legend(fontsize=8)
        ax1.set_xlim(3, 5000)

        # Panel 2: Triplet fraction
        ax2 = axes[idx, 1]
        ax2.semilogx(V_GRID, 100.0 * res['R_triplet'], 'r-', lw=1.5)
        ax2.axhline(75, color='gray', ls='--', lw=0.8, alpha=0.5,
                     label='75% (statistical)')
        ax2.axvspan(40, 80, alpha=0.15, color='orange', label='tension region')
        ax2.set_xlabel('v [km/s]', fontsize=11)
        ax2.set_ylabel('R_triplet = 3σ_odd/σ_Maj [%]', fontsize=11)
        ax2.set_title(f'{label} — triplet fraction', fontsize=12)
        ax2.legend(fontsize=9)
        ax2.set_xlim(3, 5000)
        ax2.set_ylim(0, 100)

        # Panel 3: Maj/Dir ratio
        ax3 = axes[idx, 2]
        ax3.semilogx(V_GRID, res['ratio_maj_dir'], 'k-', lw=1.5)
        ax3.axhline(1.0, color='gray', ls='--', lw=0.8)
        ax3.axvspan(40, 80, alpha=0.15, color='orange', label='tension region')
        ax3.set_xlabel('v [km/s]', fontsize=11)
        ax3.set_ylabel('σ_Maj/σ_Dir', fontsize=11)
        ax3.set_title(f'{label} — Majorana/Dirac ratio', fontsize=12)
        ax3.legend(fontsize=9)
        ax3.set_xlim(3, 5000)

    plt.tight_layout()
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'cp_channel_tension.{ext}'), dpi=200)
    plt.close(fig)
    print(f"\n  Figures saved to output/cp_channel_tension.{{png,pdf}}")

    # ──────────────────────────────────────────────────────────
    #  CSV
    # ──────────────────────────────────────────────────────────
    csv_path = os.path.join(out_dir, 'cp_channel_tension.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['benchmark', 'v_km_s', 'sigma_even', 'sigma_odd_raw',
                     'sigma_odd_x3', 'sigma_maj', 'sigma_dirac',
                     'R_triplet', 'maj_dir_ratio'])
        for label, res in all_results.items():
            for i in range(N_PTS):
                w.writerow([label, f"{V_GRID[i]:.4f}",
                            f"{res['sig_even'][i]:.6e}",
                            f"{res['sig_odd'][i]:.6e}",
                            f"{3.0 * res['sig_odd'][i]:.6e}",
                            f"{res['sig_maj'][i]:.6e}",
                            f"{res['sig_dir'][i]:.6e}",
                            f"{res['R_triplet'][i]:.6f}",
                            f"{res['ratio_maj_dir'][i]:.6f}"])
    print(f"  CSV saved to {csv_path}")

    dt = time.time() - t0
    print(f"\n  Total runtime: {dt:.1f}s")
    print(hdr)

    return all_results


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("✅ cp_channel_tension done!")
    except Exception:
        pass
