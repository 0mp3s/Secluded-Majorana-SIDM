#!/usr/bin/env python3
"""
vpm_scan/attractive_repulsive_ratio.py
======================================
Attractive vs Repulsive Yukawa cross-section ratio R(v) = sigma_att / sigma_rep
for BP1 and MAP benchmarks across SIDM-relevant velocities.

Physics:
  Attractive:  V(r) = -alpha * exp(-m_phi r) / r   [+lam in VPM ODE]
  Repulsive:   V(r) = +alpha * exp(-m_phi r) / r   [-lam in VPM ODE]

Attractive potentials support resonances (quasi-bound states) =>
  strong velocity-dependent enhancement at low v.
Repulsive potentials have no resonances =>
  monotonically decreasing sigma(v), Born expansion converges faster.

Produces:
  - Console table
  - output/attractive_repulsive_ratio.png
  - output/attractive_repulsive_ratio.pdf
"""
import sys, os, math, time
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)

from numba import jit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Constants ──
GEV2_TO_CM2 = 3.8938e-28
GEV_IN_G    = 1.78266e-24
C_KM_S      = 299792.458

# ── Load config ──
from config_loader import load_config
cfg = load_config(__file__)

_bps = {bp['label']: bp for bp in cfg.get('benchmark_points', [])}
BP1 = _bps['BP1']
MAP = _bps['MAP']

BENCHMARKS = [
    {"label": "BP1",
     "m_chi": BP1['m_chi_GeV'],
     "m_phi": BP1['m_phi_MeV'] * 1e-3,
     "alpha": BP1['alpha']},
    {"label": "MAP",
     "m_chi": MAP['m_chi_GeV'],
     "m_phi": MAP['m_phi_MeV'] * 1e-3,
     "alpha": MAP['alpha']},
]

VELOCITIES = [12, 30, 60, 100, 200, 500, 1000, 2000, 4700]

OUT_DIR = os.path.join(_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)


# ── Numba VPM solver ──
@jit(nopython=True, cache=True)
def sph_jn_numba(l, z):
    if z < 1e-30:
        return 1.0 if l == 0 else 0.0
    j0 = math.sin(z) / z
    if l == 0: return j0
    j1 = math.sin(z) / (z * z) - math.cos(z) / z
    if l == 1: return j1
    j_prev, j_curr = j0, j1
    for n in range(1, l):
        j_next = (2 * n + 1) / z * j_curr - j_prev
        j_prev, j_curr = j_curr, j_next
        if abs(j_curr) < 1e-300: return 0.0
    return j_curr

@jit(nopython=True, cache=True)
def sph_yn_numba(l, z):
    if z < 1e-30: return -1e300
    y0 = -math.cos(z) / z
    if l == 0: return y0
    y1 = -math.cos(z) / (z * z) - math.sin(z) / z
    if l == 1: return y1
    y_prev, y_curr = y0, y1
    for n in range(1, l):
        y_next = (2 * n + 1) / z * y_curr - y_prev
        y_prev, y_curr = y_curr, y_next
        if abs(y_curr) > 1e200: return y_curr
    return y_curr

@jit(nopython=True, cache=True)
def _vpm_rhs_signed(l, kappa, lam_signed, x, delta):
    """VPM RHS with signed lambda: +lam = attractive, -lam = repulsive."""
    if x < 1e-20: return 0.0
    z = kappa * x
    if z < 1e-20: return 0.0
    jl = sph_jn_numba(l, z)
    nl = sph_yn_numba(l, z)
    j_hat = z * jl
    n_hat = -z * nl
    cd = math.cos(delta)
    sd = math.sin(delta)
    bracket = j_hat * cd - n_hat * sd
    if not math.isfinite(bracket): return 0.0
    pot = lam_signed * math.exp(-x) / (kappa * x)
    val = pot * bracket * bracket
    return val if math.isfinite(val) else 0.0

@jit(nopython=True, cache=True)
def vpm_phase_shift_signed(l, kappa, lam_signed, x_max, N_steps):
    """Phase shift with signed coupling."""
    if abs(lam_signed) < 1e-30 or kappa < 1e-30: return 0.0
    x_min = max(1e-5, 0.05 / (kappa + 0.01))
    if l > 0:
        x_barrier = l / kappa
        if x_barrier > x_min:
            x_min = x_barrier
    h = (x_max - x_min) / N_steps
    delta = 0.0
    for i in range(N_steps):
        x = x_min + i * h
        k1 = _vpm_rhs_signed(l, kappa, lam_signed, x, delta)
        k2 = _vpm_rhs_signed(l, kappa, lam_signed, x + 0.5*h, delta + 0.5*h*k1)
        k3 = _vpm_rhs_signed(l, kappa, lam_signed, x + 0.5*h, delta + 0.5*h*k2)
        k4 = _vpm_rhs_signed(l, kappa, lam_signed, x + h, delta + h*k3)
        delta += h * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return delta

@jit(nopython=True, cache=True)
def sigma_T_signed(m_chi, m_phi, alpha, v_km_s, sign):
    """sigma_T/m [cm^2/g] with potential sign: +1 = attractive, -1 = repulsive."""
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    lam_signed = sign * lam
    if kappa < 1e-15: return 0.0

    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000

    l_max = min(max(3, int(kappa) + 3), 80)

    sigma_sum = 0.0
    for l in range(l_max + 1):
        delta = vpm_phase_shift_signed(l, kappa, lam_signed, x_max, N_steps)
        weight = 1.0 if l % 2 == 0 else 3.0
        contrib = weight * (2*l + 1) * math.sin(delta)**2
        sigma_sum += contrib
        if l > int(kappa) + 1 and sigma_sum > 0:
            if contrib / sigma_sum < 1e-3:
                break

    sigma_GeV2 = 2.0 * math.pi * sigma_sum / (k * k)
    sigma_cm2 = sigma_GeV2 * GEV2_TO_CM2
    return sigma_cm2 / (m_chi * GEV_IN_G)


def main():
    print("=" * 90)
    print("  Attractive vs Repulsive Yukawa -- Cross-Section Ratio R(v)")
    print("  R(v) = sigma_attractive / sigma_repulsive")
    print("=" * 90)

    # Warm up JIT
    sigma_T_signed(20.0, 10e-3, 1e-3, 100.0, +1.0)
    sigma_T_signed(20.0, 10e-3, 1e-3, 100.0, -1.0)

    results = {}  # results[label][v] = (sig_att, sig_rep, ratio)

    for bm in BENCHMARKS:
        label = bm['label']
        lam = bm['alpha'] * bm['m_chi'] / bm['m_phi']
        print(f"\n  --- {label}: m_chi={bm['m_chi']:.2f} GeV, "
              f"m_phi={bm['m_phi']*1e3:.2f} MeV, alpha={bm['alpha']:.3e}, "
              f"lambda={lam:.2f} ---")

        results[label] = {}
        print(f"  {'v [km/s]':>10}  {'sig_att/m':>12}  {'sig_rep/m':>12}  {'R=att/rep':>10}")
        print("  " + "-" * 52)

        for v in VELOCITIES:
            sig_att = sigma_T_signed(bm['m_chi'], bm['m_phi'], bm['alpha'],
                                     float(v), +1.0)
            sig_rep = sigma_T_signed(bm['m_chi'], bm['m_phi'], bm['alpha'],
                                     float(v), -1.0)
            ratio = sig_att / sig_rep if sig_rep > 0 else float('inf')
            results[label][v] = (sig_att, sig_rep, ratio)
            print(f"  {v:>10}  {sig_att:>12.6f}  {sig_rep:>12.6f}  {ratio:>10.3f}")

    # ── Publication figure ──
    print("\n  Generating attractive/repulsive ratio plot...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    colors = {'BP1': '#E91E63', 'MAP': '#2196F3'}
    markers = {'BP1': 'o', 'MAP': 's'}

    # Panel 1: sigma_att(v) and sigma_rep(v)
    ax = axes[0]
    for bm in BENCHMARKS:
        label = bm['label']
        vels = sorted(results[label].keys())
        sig_a = [results[label][v][0] for v in vels]
        sig_r = [results[label][v][1] for v in vels]
        ax.plot(vels, sig_a, '-', color=colors[label], lw=2,
                marker=markers[label], ms=5, label=f'{label} attractive')
        ax.plot(vels, sig_r, '--', color=colors[label], lw=2,
                marker=markers[label], ms=5, mfc='none', label=f'{label} repulsive')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$v$ [km/s]', fontsize=13)
    ax.set_ylabel('$\\sigma_T/m$ [cm$^2$/g]', fontsize=13)
    ax.set_title('Cross sections', fontsize=13)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')

    # Panel 2: Ratio R(v)
    ax = axes[1]
    for bm in BENCHMARKS:
        label = bm['label']
        vels = sorted(results[label].keys())
        ratios = [results[label][v][2] for v in vels]
        ax.plot(vels, ratios, '-', color=colors[label], lw=2.5,
                marker=markers[label], ms=6, label=label)
    ax.axhline(1.0, color='gray', ls=':', lw=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$v$ [km/s]', fontsize=13)
    ax.set_ylabel('$R(v) = \\sigma_{\\rm att}/\\sigma_{\\rm rep}$', fontsize=13)
    ax.set_title('Attractive / Repulsive ratio', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

    # Panel 3: Individual phase shifts l=0,1,2 for BP1 at sample velocity
    ax = axes[2]
    bm = BENCHMARKS[0]  # BP1
    v_sample = 60.0
    v = v_sample / C_KM_S
    mu = bm['m_chi'] / 2.0
    k = mu * v
    kappa = k / bm['m_phi']
    lam = bm['alpha'] * bm['m_chi'] / bm['m_phi']
    l_max = min(max(3, int(kappa) + 3), 30)
    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000

    ells = list(range(l_max + 1))
    d_att = [vpm_phase_shift_signed(l, kappa, +lam, x_max, N_steps) for l in ells]
    d_rep = [vpm_phase_shift_signed(l, kappa, -lam, x_max, N_steps) for l in ells]

    ax.bar(np.array(ells) - 0.15, [abs(d) for d in d_att], 0.3, color='#E91E63',
           alpha=0.7, label='Attractive')
    ax.bar(np.array(ells) + 0.15, [abs(d) for d in d_rep], 0.3, color='#2196F3',
           alpha=0.7, label='Repulsive')
    ax.set_xlabel('Partial wave $\\ell$', fontsize=13)
    ax.set_ylabel('$|\\delta_\\ell|$ [rad]', fontsize=13)
    ax.set_title(f'Phase shifts, BP1, $v={int(v_sample)}$ km/s', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle('Attractive vs Repulsive Yukawa: Secluded Majorana SIDM',
                 fontsize=14, y=1.02)
    fig.tight_layout()

    for ext in ['png', 'pdf']:
        path = os.path.join(OUT_DIR, f'attractive_repulsive_ratio.{ext}')
        fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: output/attractive_repulsive_ratio.png")
    print(f"  Saved: output/attractive_repulsive_ratio.pdf")

    # ── Summary ──
    print(f"\n{'='*90}")
    print("  SUMMARY")
    print(f"{'='*90}")
    for bm in BENCHMARKS:
        label = bm['label']
        r_low = results[label][30][2]
        r_high = results[label][1000][2]
        print(f"  {label}: R(30 km/s) = {r_low:.2f}, R(1000 km/s) = {r_high:.2f}")
    print(f"\n  Attractive Yukawa resonances drive R >> 1 at low velocities.")
    print(f"  At high v, both converge toward Born regime => R -> 1.")


if __name__ == '__main__':
    t0 = time.time()
    main()
    print(f"\n  Total runtime: {time.time()-t0:.1f}s")
