#!/usr/bin/env python3
"""
cross_checks/averaging_error.py
================================
Velocity-averaging error  ε(v₀) = ⟨σ/m⟩_MB / σ(v₀) - 1

Physics:
  Observational constraints quote σ/m at a single representative velocity
  v₀ (e.g., 30 km/s for dSphs). But real halos have a Maxwell-Boltzmann
  velocity distribution. The averaging error ε quantifies the systematic
  bias from this point-estimate approximation:

      ε(v₀) = [ ⟨σ/m⟩_MB(v₀) - σ(v₀) ] / σ(v₀)

  where ⟨σ/m⟩_MB = ∫ σ(v) · v · f_MB(v; v₀/√2) dv  /  ∫ v · f_MB dv

  Key regimes:
    ε ≈ 0    point estimate is accurate — safe to compare theory vs data
    ε > 0    averaging enhances σ — point estimate is conservative
    ε < 0    averaging suppresses σ — point estimate is optimistic
    |ε| >> 1 near resonance nodes — σ(v₀) is completely misleading

  Relation to velocity slope η:
    For a power law σ ∝ v^η, the averaging correction is analytic:
      ⟨σ v⟩/⟨v⟩ = σ(v₀) · Γ((η+4)/2) · 2^{η/2} / Γ(2)
    So ε depends on η: flat regions (η≈0) have small ε; steep regions
    (|η| >> 1) have large ε. This script verifies and extends that
    relationship numerically.

  Practical consequence:
    For each observation in our dataset, we compute ε to quantify the
    systematic uncertainty from assuming σ(v_obs) = ⟨σ/m⟩. This is
    needed to determine if the chi-squared fit is valid at face value
    or requires a velocity-averaging correction.

Produces:
  - Console: ε table at observational velocities, correction factors
  - output/averaging_error.png  (3-panel: σ(v), ε(v), η-vs-ε scatter)
  - output/averaging_error.pdf
  - output/averaging_error.csv
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

# Warm up JIT (before threading)
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

cfg = load_config(__file__)
_bps = {bp['label']: bp for bp in GC.benchmarks_from_labels(cfg.get('benchmark_labels', ['BP1', 'MAP']))}
BP1 = _bps['BP1']
MAP = _bps['MAP']

# Load observations for reference
OBSERVATIONS = GC.observations_as_tuples()


# ══════════════════════════════════════════════════════════════
#  Pre-compute σ(v) on dense grid with multithreading,
#  then use interpolation for MB averaging — fast path.
# ══════════════════════════════════════════════════════════════
N_WORKERS = os.cpu_count() or 4

# Dense grid for interpolation (covers [1, 25000] to handle MB tail at 5×v_max)
N_DENSE = 500
V_DENSE = np.logspace(0, np.log10(25000.0), N_DENSE)


def _build_sigma_interp(m_chi, m_phi_GeV, alpha):
    """Compute σ(v) on dense grid using thread pool, return interpolator."""
    def _eval(v):
        return sigma_T_vpm(m_chi, m_phi_GeV, alpha, v)

    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        sig_dense = np.array(list(pool.map(_eval, V_DENSE)))

    # Log-log interpolation for smooth cross section
    valid = sig_dense > 0
    if np.sum(valid) < 10:
        return lambda v: 0.0, sig_dense
    log_interp = interp1d(np.log(V_DENSE[valid]), np.log(sig_dense[valid]),
                          kind='cubic', fill_value='extrapolate')
    def sigma_interp(v):
        return np.exp(log_interp(np.log(v)))
    return sigma_interp, sig_dense


def sigma_m_averaged_fast(sigma_func, v_char_km_s, n_gauss=30):
    """
    ⟨σ/m⟩ averaged over MB distribution using pre-computed interpolator.
    """
    v0 = v_char_km_s / np.sqrt(2.0)

    v_max = max(5.0 * v_char_km_s, 100.0)
    a, b = 1.0, v_max

    nodes, weights_gl = np.polynomial.legendre.leggauss(n_gauss)
    v_pts = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    jac = 0.5 * (b - a)

    sigma_vals = sigma_func(v_pts)  # vectorized interpolation
    w_vals = v_pts**3 * np.exp(-(v_pts / v0)**2 / 2.0)

    num = jac * np.sum(weights_gl * sigma_vals * w_vals)
    den = jac * np.sum(weights_gl * w_vals)

    if den < 1e-300:
        return 0.0
    return num / den


# ══════════════════════════════════════════════════════════════
#  Velocity grid and computation
# ══════════════════════════════════════════════════════════════
N_PTS = 150
V_MIN, V_MAX = 5.0, 5000.0
V_GRID = np.logspace(np.log10(V_MIN), np.log10(V_MAX), N_PTS)

# Observational reference velocities
V_OBS = sorted(set(v for _, v, *_ in OBSERVATIONS))


def compute_epsilon(bp):
    """Compute ε(v₀), σ(v₀), ⟨σ/m⟩_MB(v₀), and η(v₀) on the velocity grid."""
    m_chi = bp['m_chi_GeV']
    m_phi = bp['m_phi_MeV'] / 1000.0
    alpha = bp['alpha']

    # Step 1: build interpolator (dense grid + threading)
    sigma_func, _ = _build_sigma_interp(m_chi, m_phi, alpha)

    # Step 2: evaluate σ(v₀) at grid points (fast — from interpolator)
    sig_fixed = np.array([sigma_T_vpm(m_chi, m_phi, alpha, v) for v in V_GRID])

    # Step 3: MB averaging via interpolation (very fast)
    sig_avg = np.array([sigma_m_averaged_fast(sigma_func, v0) for v0 in V_GRID])

    epsilon = np.where(sig_fixed > 1e-30, sig_avg / sig_fixed - 1.0, np.nan)

    # Compute η via finite differences on log-log grid
    log_v = np.log(V_GRID)
    log_sig = np.where(sig_fixed > 0, np.log(sig_fixed), np.nan)
    eta = np.gradient(log_sig, log_v)

    return sig_fixed, sig_avg, epsilon, eta


def main():
    t0 = time.time()
    print("=" * 70)
    print("  Velocity-Averaging Error  ε(v₀) = ⟨σ/m⟩_MB / σ(v₀) - 1")
    print("=" * 70)

    results = {}
    for label, bp in [('BP1', BP1), ('MAP', MAP)]:
        lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000.0)
        print(f"\n  {label}:  m_χ={bp['m_chi_GeV']:.2f} GeV, "
              f"m_φ={bp['m_phi_MeV']:.2f} MeV, α={bp['alpha']:.4e}, λ={lam:.2f}")
        t1 = time.time()
        sig_fixed, sig_avg, epsilon, eta = compute_epsilon(bp)
        dt = time.time() - t1
        print(f"    Computed {N_PTS} velocities in {dt:.1f}s")
        results[label] = (sig_fixed, sig_avg, epsilon, eta)

        # Summary at observational velocities
        print(f"\n    {'System':<22s}  {'v₀':>6s}  {'σ(v₀)':>10s}  "
              f"{'⟨σ/m⟩_MB':>10s}  {'ε':>8s}  {'η':>6s}  {'Comment':s}")
        print("    " + "-" * 85)

        for name, v_obs, central, lo, hi, ref in OBSERVATIONS:
            # Interpolate from grid
            idx = np.argmin(np.abs(V_GRID - v_obs))
            sf = sig_fixed[idx]
            sa = sig_avg[idx]
            eps = epsilon[idx]
            et = eta[idx]

            # Flag significance
            if abs(eps) < 0.05:
                comment = "negligible"
            elif abs(eps) < 0.15:
                comment = "small correction"
            elif abs(eps) < 0.50:
                comment = "SYSTEMATIC"
            else:
                comment = "DANGEROUS"

            print(f"    {name:<22s}  {v_obs:>6d}  {sf:>10.4f}  "
                  f"{sa:>10.4f}  {eps:>+8.1%}  {et:>+6.2f}  {comment}")

    # ──────────────────────────────────────────────────────────
    #  3-panel figure
    # ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    colors = {'BP1': 'C0', 'MAP': 'C1'}
    obs_colors = {'dSph': 'green', 'TBTF': 'olive', 'LSB': 'teal',
                  'MW': 'orange', 'Cluster': 'red', 'Bullet': 'purple'}

    # Panel 1: σ(v₀) vs ⟨σ/m⟩_MB(v₀)
    ax1 = axes[0]
    for label in ['BP1', 'MAP']:
        sig_fixed, sig_avg, _, _ = results[label]
        ax1.loglog(V_GRID, sig_fixed, color=colors[label], ls='-', lw=2,
                   label=f'σ(v₀) {label}')
        ax1.loglog(V_GRID, sig_avg, color=colors[label], ls='--', lw=1.5,
                   label=f'⟨σ/m⟩_MB {label}')
    ax1.set_ylabel('σ_T/m [cm²/g]', fontsize=12)
    ax1.set_title('Point Estimate vs Maxwell-Boltzmann Average', fontsize=13)
    ax1.legend(fontsize=9)
    # Mark observational velocities
    for v_obs in V_OBS:
        ax1.axvline(v_obs, color='grey', alpha=0.2, ls=':', lw=0.8)

    # Panel 2: ε(v₀) — the averaging error as continuous function
    ax2 = axes[1]
    for label in ['BP1', 'MAP']:
        _, _, epsilon, _ = results[label]
        ax2.plot(V_GRID, epsilon * 100, color=colors[label], lw=2, label=label)
    ax2.axhline(0, color='grey', ls=':', lw=1)
    ax2.axhspan(-5, 5, color='green', alpha=0.1, label='|ε| < 5%')
    ax2.axhspan(-15, -5, color='yellow', alpha=0.08)
    ax2.axhspan(5, 15, color='yellow', alpha=0.08)
    ax2.set_xscale('log')
    ax2.set_ylabel('ε(v₀) [%]', fontsize=12)
    ax2.set_title('Velocity-Averaging Error ε(v₀) = ⟨σ/m⟩/σ(v₀) - 1', fontsize=13)
    ax2.legend(fontsize=10)
    # Mark each observation with a dot
    for label in ['BP1', 'MAP']:
        _, _, epsilon, _ = results[label]
        for v_obs in V_OBS:
            idx = np.argmin(np.abs(V_GRID - v_obs))
            ax2.plot(v_obs, epsilon[idx] * 100, 'o', color=colors[label],
                     ms=5, zorder=5)
    for v_obs in V_OBS:
        ax2.axvline(v_obs, color='grey', alpha=0.2, ls=':', lw=0.8)

    # Panel 3: η vs ε scatter — verify the theoretical relationship
    ax3 = axes[2]
    for label in ['BP1', 'MAP']:
        _, _, epsilon, eta = results[label]
        valid = ~np.isnan(epsilon) & ~np.isnan(eta)
        ax3.scatter(eta[valid], epsilon[valid] * 100, c=colors[label],
                    s=8, alpha=0.5, label=label)
    ax3.axhline(0, color='grey', ls=':', lw=1)
    ax3.axvline(0, color='grey', ls=':', lw=1)
    ax3.set_xlabel('η = d ln(σ/m) / d ln v', fontsize=12)
    ax3.set_ylabel('ε [%]', fontsize=12)
    ax3.set_title('Velocity Slope η vs Averaging Error ε (correlation check)', fontsize=13)
    ax3.legend(fontsize=10)

    # Overlay theoretical power-law prediction: for σ∝v^η, ε_theory
    eta_theory = np.linspace(-4, 2, 200)
    # ⟨v^η · v⟩/⟨v⟩ for MB distribution:
    # = 2^{η/2} · Γ((η+4)/2) / Γ(2)  [Γ(2) = 1]
    from scipy.special import gamma
    eps_theory = 2.0**(eta_theory / 2.0) * gamma((eta_theory + 4) / 2.0) / gamma(2.0) - 1.0
    ax3.plot(eta_theory, eps_theory * 100, 'k--', lw=1.5, alpha=0.7,
             label='Power-law prediction')
    ax3.legend(fontsize=9)

    plt.tight_layout()
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'averaging_error.{ext}'), dpi=200)
    plt.close(fig)
    print(f"\n  Figures saved to output/averaging_error.{{png,pdf}}")

    # ──────────────────────────────────────────────────────────
    #  CSV output
    # ──────────────────────────────────────────────────────────
    csv_path = os.path.join(out_dir, 'averaging_error.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['v0_km_s',
                     'BP1_sigma_fixed', 'BP1_sigma_avg', 'BP1_epsilon', 'BP1_eta',
                     'MAP_sigma_fixed', 'MAP_sigma_avg', 'MAP_epsilon', 'MAP_eta'])
        for i in range(N_PTS):
            row = [f'{V_GRID[i]:.4f}']
            for label in ['BP1', 'MAP']:
                sf, sa, eps, eta = results[label]
                row.extend([
                    f'{sf[i]:.8e}', f'{sa[i]:.8e}',
                    f'{eps[i]:.8e}', f'{eta[i]:.8e}'
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
        sig_fixed, sig_avg, epsilon, eta = results[label]
        valid = ~np.isnan(epsilon)

        # Max |ε| and where it occurs
        abs_eps = np.abs(epsilon[valid])
        max_idx = np.argmax(abs_eps)
        v_max_eps = V_GRID[valid][max_idx]
        max_eps = epsilon[valid][max_idx]

        # ε at each observational velocity
        obs_eps = []
        for v_obs in V_OBS:
            idx = np.argmin(np.abs(V_GRID - v_obs))
            obs_eps.append((v_obs, epsilon[idx]))

        # Count how many observations need correction (|ε| > 5%)
        n_dangerous = sum(1 for _, e in obs_eps if abs(e) > 0.05)

        print(f"\n  {label}:")
        print(f"    Max |ε|:  {max_eps:+.1%}  at v₀ = {v_max_eps:.1f} km/s")
        print(f"    Observations needing correction (|ε| > 5%):  {n_dangerous} / {len(V_OBS)}")

        # Per-observation ε
        for v_obs, eps in obs_eps:
            flag = " ⚠" if abs(eps) > 0.05 else ""
            print(f"      v₀ = {v_obs:>5d} km/s:  ε = {eps:+.1%}{flag}")

        # Correlation between η and ε
        valid_both = ~np.isnan(epsilon) & ~np.isnan(eta)
        if np.sum(valid_both) > 10:
            corr = np.corrcoef(eta[valid_both], epsilon[valid_both])[0, 1]
            print(f"    η–ε correlation:  r = {corr:.4f}")

    dt_total = time.time() - t0
    print(f"\n  Total runtime: {dt_total:.1f}s")
    print("=" * 70)


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("✅ averaging_error done!")
    except Exception:
        pass
