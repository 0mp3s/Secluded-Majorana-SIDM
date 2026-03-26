#!/usr/bin/env python3
"""
predictions/gravothermal/n_scatter_function.py
===============================================
N_scatter(v) as a continuous function — find "protected velocities" v*
where dN/dv = 0.

Physics:
  The scattering rate per Hubble time is:
      N(v) = (σ/m)(v) · ρ_s · v · t_age

  The derivative:
      dN/dv = ρ_s · t · [v · d(σ/m)/dv + σ/m]
            = ρ_s · t · (σ/m) · [η(v) + 1]

  where η(v) = d ln(σ/m) / d ln v is the logarithmic slope.

  Critical points:
      dN/dv = 0  ⟺  η(v*) = -1

  At v*, N_scatter has a LOCAL MINIMUM — halos with velocity dispersion
  near v* are "SIDM-protected": they experience the fewest scatterings
  and therefore retain their NFW cusp longest.

  Prediction: if a specific class of galaxies (with σ_v ≈ v*) shows cusps
  while neighbours with higher/lower σ_v show cores, this is a direct
  signature of velocity-dependent SIDM with resonant structure.

Produces:
  - Console: N(v) table, v* locations, regime classification per dSph
  - output/n_scatter_function.png  (2-panel)
  - output/n_scatter_function.pdf
  - output/n_scatter_function.csv
"""
import sys, os, math, time, csv
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

cfg = load_config(__file__)
_bps = {bp['label']: bp for bp in GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1', 'MAP']))}
BP1 = _bps['BP1']
MAP = _bps['MAP']

# ── Physical constants (CGS) ──
SEC_PER_GYR = 3.156e16
KM_S_TO_CM_S = 1e5
MSUN_G = 1.989e33
KPC_CM = 3.086e21

# ── Velocity grid ──
N_PTS = 500
V_MIN, V_MAX = 3.0, 5000.0
V_GRID = np.logspace(np.log10(V_MIN), np.log10(V_MAX), N_PTS)

# ── Representative dSph parameters (from dsphs_data.csv) ──
# Use Fornax as the fiducial halo — well-measured, clearly cored
FIDUCIAL = {
    'name': 'Fornax',
    'rho_s_Msun_kpc3': 5.0e7,
    't_age_Gyr': 12.0,
}


def rho_cgs(rho_Msun_kpc3):
    """Convert ρ [M_sun/kpc³] → [g/cm³]."""
    return rho_Msun_kpc3 * MSUN_G / KPC_CM ** 3


def compute_n_scatter(bp, rho_s, t_age_Gyr, v_grid):
    """Compute N_scatter(v) = (σ/m)(v) × ρ_s × v × t_age for each v in grid."""
    m_chi = bp['m_chi_GeV']
    m_phi = bp['m_phi_MeV'] / 1000.0
    alpha = bp['alpha']
    rho = rho_cgs(rho_s)
    t = t_age_Gyr * SEC_PER_GYR

    sigma = np.array([sigma_T_vpm(m_chi, m_phi, alpha, v) for v in v_grid])
    v_cm = v_grid * KM_S_TO_CM_S
    n_scatter = sigma * rho * v_cm * t
    return sigma, n_scatter


def find_protected_velocities(v_grid, n_scatter):
    """Find local minima of N(v) — "protected velocities" where η = -1."""
    log_n = np.log(np.maximum(n_scatter, 1e-30))
    log_v = np.log(v_grid)
    dn_dlogv = np.gradient(log_n, log_v)  # = η + 1

    # Local minima: dn_dlogv crosses zero from negative to positive
    minima = []
    for i in range(1, len(dn_dlogv) - 1):
        if dn_dlogv[i - 1] < 0 and dn_dlogv[i + 1] > 0:
            # Interpolate zero crossing
            f = -dn_dlogv[i - 1] / (dn_dlogv[i + 1] - dn_dlogv[i - 1])
            v_star = np.exp(log_v[i - 1] + f * (log_v[i + 1] - log_v[i - 1]))
            n_star = np.exp(np.interp(np.log(v_star), log_v, log_n))
            minima.append((v_star, n_star))

    # Also find local maxima (for context)
    maxima = []
    for i in range(1, len(dn_dlogv) - 1):
        if dn_dlogv[i - 1] > 0 and dn_dlogv[i + 1] < 0:
            f = dn_dlogv[i - 1] / (dn_dlogv[i - 1] - dn_dlogv[i + 1])
            v_peak = np.exp(log_v[i - 1] + f * (log_v[i + 1] - log_v[i - 1]))
            n_peak = np.exp(np.interp(np.log(v_peak), log_v, log_n))
            maxima.append((v_peak, n_peak))

    return minima, maxima


def load_dsphs():
    """Load dSph data."""
    csv_path = os.path.join(_DIR, cfg.get('dsphs_csv', 'dsphs_data.csv'))
    dsphs = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            dsphs.append({
                'name': row['name'],
                'sigma_v': float(row['sigma_v_km_s']),
                'rho_s': float(row['rho_s_Msun_kpc3']),
                't_age_Gyr': float(row['t_age_Gyr']),
                'core_observed': row['core_observed'],
            })
    return dsphs


def main():
    t0 = time.time()
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)
    dsphs = load_dsphs()

    benchmarks = {'BP1': BP1, 'MAP': MAP}
    results = {}

    print("=" * 80)
    print("  N_scatter(v) — Continuous Function Analysis")
    print("  Finding 'protected velocities' v* where dN/dv = 0 (η = -1)")
    print("=" * 80)

    for label, bp in benchmarks.items():
        lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000.0)
        print(f"\n  --- {label} (λ={lam:.1f}) with Fornax halo (ρ_s={FIDUCIAL['rho_s_Msun_kpc3']:.1e} M☉/kpc³) ---")

        sigma, n_scatter = compute_n_scatter(
            bp, FIDUCIAL['rho_s_Msun_kpc3'], FIDUCIAL['t_age_Gyr'], V_GRID)

        # η(v) = d ln(σ/m) / d ln v
        log_v = np.log(V_GRID)
        log_sigma = np.log(np.maximum(sigma, 1e-30))
        eta = np.gradient(log_sigma, log_v)

        results[label] = {
            'sigma': sigma, 'n_scatter': n_scatter,
            'eta': eta, 'lambda': lam,
        }

        # Find protected velocities
        minima, maxima = find_protected_velocities(V_GRID, n_scatter)
        results[label]['minima'] = minima
        results[label]['maxima'] = maxima

        if minima:
            print(f"\n    PROTECTED VELOCITIES (local minima of N_scatter):")
            for v_star, n_star in minima:
                regime = "CUSPY" if n_star < 1 else ("CORED" if n_star < 100 else "COLLAPSE")
                print(f"      v* = {v_star:.1f} km/s → N_scatter = {n_star:.2f} ({regime})")
                # Which dSphs are near this v*?
                for d in dsphs:
                    v_rel = d['sigma_v'] * math.sqrt(2)
                    if abs(np.log(v_rel / v_star)) < 0.3:  # within factor ~1.3
                        print(f"        → {d['name']} (v_rel={v_rel:.1f}) is near v*! "
                              f"Core observed: {d['core_observed']}")
        else:
            print(f"\n    No protected velocities found — N(v) monotonic in {V_MIN}–{V_MAX} km/s")

        if maxima:
            print(f"\n    PEAK SCATTERING VELOCITIES (local maxima of N_scatter):")
            for v_peak, n_peak in maxima:
                regime = "CUSPY" if n_peak < 1 else ("CORED" if n_peak < 100 else "COLLAPSE")
                print(f"      v_peak = {v_peak:.1f} km/s → N_scatter = {n_peak:.2f} ({regime})")

        # dSph predictions
        print(f"\n    dSph predictions at actual velocities:")
        print(f"    {'dSph':<12s} {'v_rel':>6s} {'σ/m':>8s} {'N_scat':>8s} {'η':>6s} {'Regime':>10s}")
        print("    " + "-" * 56)
        for d in dsphs:
            v_rel = d['sigma_v'] * math.sqrt(2)
            idx = np.argmin(np.abs(V_GRID - v_rel))
            s = sigma[idx]
            n = n_scatter[idx]
            e = eta[idx]
            regime = "CUSPY" if n < 1 else ("CORED" if n < 100 else "COLLAPSE")
            print(f"    {d['name']:<12s} {v_rel:>6.1f} {s:>8.3f} {n:>8.2f} {e:>6.2f} {regime:>10s}")

    # ── Plot: 2-panel ──
    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    colors = {'BP1': 'C0', 'MAP': 'C3'}

    # Panel (a): N_scatter(v) with regime bands
    ax = axes[0]
    ax.axhspan(0, 1, alpha=0.05, color='blue', label='CUSPY (N<1)')
    ax.axhspan(1, 100, alpha=0.08, color='green', label='CORED (1<N<100)')
    ax.axhspan(100, 1e5, alpha=0.05, color='red', label='COLLAPSE (N>100)')
    ax.axhline(1, ls='--', color='green', lw=0.7)
    ax.axhline(100, ls='--', color='red', lw=0.7)

    for label in benchmarks:
        lam = results[label]['lambda']
        ax.loglog(V_GRID, results[label]['n_scatter'], '-',
                  color=colors[label], lw=1.8,
                  label=f'{label} (λ={lam:.1f})')
        # Mark protected velocities
        for v_star, n_star in results[label]['minima']:
            ax.plot(v_star, n_star, 'v', color=colors[label], ms=12, zorder=5)
            ax.annotate(f'v*={v_star:.0f}', (v_star, n_star),
                        textcoords='offset points', xytext=(10, -15),
                        fontsize=8, color=colors[label])
        # Mark peak velocities
        for v_peak, n_peak in results[label]['maxima']:
            ax.plot(v_peak, n_peak, '^', color=colors[label], ms=10, zorder=5)

    # Mark actual dSphs
    for d in dsphs:
        v_rel = d['sigma_v'] * math.sqrt(2)
        ax.axvline(v_rel, ls=':', color='gray', lw=0.5, alpha=0.5)

    ax.set_ylabel(r'$N_{\rm scatter}(v) = (\sigma/m) \rho_s v \, t_{\rm age}$')
    ax.set_ylim(0.1, 1e4)
    ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.set_title(r'(a) $N_{\rm scatter}(v)$ — Fornax halo parameters')
    ax.grid(True, alpha=0.3, which='both')

    # Panel (b): η(v) with η = -1 line highlighted
    ax = axes[1]
    for label in benchmarks:
        ax.semilogx(V_GRID, results[label]['eta'], '-',
                     color=colors[label], lw=1.5, label=f'{label}')
        # Mark η = -1 crossings (= protected velocities)
        for v_star, _ in results[label]['minima']:
            ax.plot(v_star, -1, 'v', color=colors[label], ms=12, zorder=5)

    ax.axhline(-1, ls='-', color='orange', lw=1.5, alpha=0.7,
               label=r'$\eta = -1$ (protected velocity)')
    ax.axhline(0, ls='--', color='k', lw=0.5)
    ax.axhline(-2, ls=':', color='gray', lw=0.8, label=r'$\eta=-2$ (perturbative)')
    ax.set_ylabel(r'$\eta(v) = d\ln(\sigma/m) / d\ln v$')
    ax.set_xlabel(r'$v_{\rm rel}$ [km/s]')
    ax.set_ylim(-5, 3)
    ax.legend(fontsize=8, ncol=2, loc='lower left')
    ax.set_title(r'(b) Logarithmic slope — $\eta = -1$ marks $N_{\rm scatter}$ extrema')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'n_scatter_function.{ext}'), dpi=200)
    plt.close()
    print(f"\n  Saved: output/n_scatter_function.png")
    print(f"  Saved: output/n_scatter_function.pdf")

    # ── CSV ──
    csv_out = os.path.join(out_dir, 'n_scatter_function.csv')
    with open(csv_out, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['v_km_s']
        for label in benchmarks:
            header += [f'sigma_m_{label}', f'N_scatter_{label}', f'eta_{label}']
        writer.writerow(header)
        for i, v in enumerate(V_GRID):
            row = [f'{v:.4f}']
            for label in benchmarks:
                row.append(f'{results[label]["sigma"][i]:.6e}')
                row.append(f'{results[label]["n_scatter"][i]:.4f}')
                row.append(f'{results[label]["eta"][i]:.4f}')
            writer.writerow(row)
    print(f"  Saved: output/n_scatter_function.csv")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("✅ n_scatter_function done!")
    except Exception:
        pass
