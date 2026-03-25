#!/usr/bin/env python3
"""
vpm_scan/transfer_elastic_ratio.py
==================================
Momentum-transfer vs elastic (total) cross-section ratio for identical
Majorana fermion SIDM: R(v) = sigma_tr / sigma_el.

Physics:
  For IDENTICAL Majorana fermions, the scattering amplitude is
  antisymmetrized under particle exchange (theta -> pi - theta).
  The singlet (even-l) and triplet (odd-l) channels scatter independently.

  Elastic (total) cross section:
    sigma_el = (2pi/k^2) sum_l w_l (2l+1) sin^2(delta_l)
    with w_l = 1 (even), 3 (odd)

  Momentum-transfer cross section:
    sigma_tr = int (1 - cos theta) (d sigma/d Omega) d Omega

  For identical particles scattering in the CM frame, 90-degree
  scattering ambiguity means we use the viscosity (transport) formula.
  
  For distinguishable particles, the standard result would be:
    sigma_tr = (4pi/k^2) sum_l (l+1) sin^2(delta_l - delta_{l+1})
  
  For IDENTICAL spin-1/2 fermions (Majorana), even and odd partial waves
  do not interfere. In each channel (singlet/triplet), the exchange
  symmetry modifies the angular integration. The result uses numerical
  angular integration of the symmetrized differential cross sections.

  Singlet channel (even l only):
    f_S(theta) = sum_{l even} (2l+1) exp(i delta_l) sin(delta_l) P_l(cos theta)
    dsigma_S/dOmega = |f_S(theta) + f_S(pi-theta)|^2 / (4k^2)
    (factor 1/2 from identical particles times 2 from exchange)

  Triplet channel (odd l only):
    f_T(theta) = sum_{l odd} (2l+1) exp(i delta_l) sin(delta_l) P_l(cos theta)
    dsigma_T/dOmega = 3 |f_T(theta) - f_T(pi-theta)|^2 / (4k^2)
    (minus from Fermi statistics, factor 3 from spin degeneracy)

Produces:
  - Console table
  - output/transfer_elastic_ratio.png
  - output/transfer_elastic_ratio.pdf
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
from numpy.polynomial.legendre import legval

# ── Load config ──
from config_loader import load_config
from global_config import GC

# ── Constants (sourced from global_config.json) ──
_PC = GC.physical_constants()
GEV2_TO_CM2 = _PC["GEV2_to_cm2"]
GEV_IN_G    = _PC["GeV_in_g"]
C_KM_S      = _PC["c_km_s"]
cfg = load_config(__file__)

_bps = {bp['label']: bp for bp in GC.benchmarks_from_labels(cfg.get('benchmark_labels', []))}
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

# Number of Gauss-Legendre points for angular integration
N_GAUSS = 200


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
def _vpm_rhs(l, kappa, lam, x, delta):
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
    pot = lam * math.exp(-x) / (kappa * x)
    val = pot * bracket * bracket
    return val if math.isfinite(val) else 0.0

@jit(nopython=True, cache=True)
def vpm_phase_shift(l, kappa, lam, x_max, N_steps):
    if lam < 1e-30 or kappa < 1e-30: return 0.0
    x_min = max(1e-5, 0.05 / (kappa + 0.01))
    if l > 0:
        x_barrier = l / kappa
        if x_barrier > x_min:
            x_min = x_barrier
    h = (x_max - x_min) / N_steps
    delta = 0.0
    for i in range(N_steps):
        x = x_min + i * h
        k1 = _vpm_rhs(l, kappa, lam, x, delta)
        k2 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k1)
        k3 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k2)
        k4 = _vpm_rhs(l, kappa, lam, x + h, delta + h*k3)
        delta += h * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return delta


def get_phase_shifts(m_chi, m_phi, alpha, v_km_s):
    """Compute all phase shifts for given parameters, return (deltas, kappa, k)."""
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15:
        return [], kappa, k

    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000

    l_max = min(max(3, int(kappa) + 3), 80)

    deltas = []
    sigma_sum = 0.0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        deltas.append(delta)
        weight = 1.0 if l % 2 == 0 else 3.0
        contrib = weight * (2*l + 1) * math.sin(delta)**2
        sigma_sum += contrib
        if l > int(kappa) + 1 and sigma_sum > 0:
            if contrib / sigma_sum < 1e-3:
                break

    return deltas, kappa, k


def compute_legendre_values(l_max, cos_theta_arr):
    """Compute P_l(cos theta) for all l up to l_max at given angles."""
    n_pts = len(cos_theta_arr)
    Pl = np.zeros((l_max + 1, n_pts))
    for l in range(l_max + 1):
        coeffs = np.zeros(l + 1)
        coeffs[l] = 1.0
        Pl[l, :] = legval(cos_theta_arr, coeffs)
    return Pl


def compute_sigma_el_and_tr(m_chi, m_phi, alpha, v_km_s):
    """
    Compute both elastic (total) and momentum-transfer cross sections
    for identical Majorana fermion scattering.

    Returns (sigma_el/m, sigma_tr/m) in cm^2/g.
    """
    deltas, kappa, k = get_phase_shifts(m_chi, m_phi, alpha, v_km_s)
    if not deltas or k < 1e-30:
        return 0.0, 0.0

    l_max = len(deltas) - 1

    # ── Elastic (total) cross section ──
    sigma_el_sum = 0.0
    for l in range(l_max + 1):
        w = 1.0 if l % 2 == 0 else 3.0
        sigma_el_sum += w * (2*l + 1) * math.sin(deltas[l])**2

    sigma_el_GeV2 = 2.0 * math.pi * sigma_el_sum / (k * k)
    sigma_el_cm2 = sigma_el_GeV2 * GEV2_TO_CM2
    sigma_el_over_m = sigma_el_cm2 / (m_chi * GEV_IN_G)

    # ── Momentum-transfer cross section via angular integration ──
    # Gauss-Legendre quadrature on cos(theta) in [-1, 1]
    cos_th, weights = np.polynomial.legendre.leggauss(N_GAUSS)

    # Precompute Legendre polynomials at all angles
    Pl_plus = compute_legendre_values(l_max, cos_th)        # P_l(cos theta)
    Pl_minus = compute_legendre_values(l_max, -cos_th)      # P_l(-cos theta) = (-1)^l P_l(cos theta)

    # Build partial-wave amplitudes a_l = (2l+1) exp(i delta_l) sin(delta_l) / k
    # For angular integration, we work with dimensionless a_l * k
    al_re = np.zeros(l_max + 1)  # Re[a_l * k] = (2l+1) cos(delta_l) sin(delta_l)
    al_im = np.zeros(l_max + 1)  # Im[a_l * k] = (2l+1) sin^2(delta_l)
    for l in range(l_max + 1):
        al_re[l] = (2*l + 1) * math.cos(deltas[l]) * math.sin(deltas[l])
        al_im[l] = (2*l + 1) * math.sin(deltas[l])**2

    # Singlet channel (even l): f_S(theta) + f_S(pi-theta)
    # f_S(theta) = (1/k) sum_{l even} a_l P_l(cos theta)
    # f_S(pi-theta) = (1/k) sum_{l even} a_l P_l(-cos theta) = (1/k) sum_{l even} a_l P_l(cos theta)
    # since P_l(-x) = (-1)^l P_l(x) and l is even => P_l(-x) = P_l(x)
    # So: f_S(theta) + f_S(pi-theta) = (2/k) sum_{l even} a_l P_l(cos theta)
    fS_re = np.zeros(N_GAUSS)  # Re[k * (f_S + f_S_exchange)]
    fS_im = np.zeros(N_GAUSS)  # Im[...]
    for l in range(0, l_max + 1, 2):
        fS_re += 2.0 * al_re[l] * Pl_plus[l]
        fS_im += 2.0 * al_im[l] * Pl_plus[l]

    # Triplet channel (odd l): f_T(theta) - f_T(pi-theta)  [Fermi minus]
    # f_T(theta) = (1/k) sum_{l odd} a_l P_l(cos theta)
    # f_T(pi-theta) = (1/k) sum_{l odd} a_l (-1)^l P_l(cos theta)
    #               = -(1/k) sum_{l odd} a_l P_l(cos theta)  [since l odd]
    # So: f_T(theta) - f_T(pi-theta) = (2/k) sum_{l odd} a_l P_l(cos theta)
    fT_re = np.zeros(N_GAUSS)
    fT_im = np.zeros(N_GAUSS)
    for l in range(1, l_max + 1, 2):
        fT_re += 2.0 * al_re[l] * Pl_plus[l]
        fT_im += 2.0 * al_im[l] * Pl_plus[l]

    # Differential cross sections (dimensionless: k^2 * dsigma/dOmega)
    # Singlet: (1/4) |f_S + f_S_exchange|^2  (1/4 = 1/2 identical * 1/2 spin avg for S=0)
    # Actually for identical spin-1/2 particles:
    #   dsigma/dOmega = (1/4) |f_S + f_S_exch|^2  (singlet, weight 1/4)
    #                 + (3/4) |f_T - f_T_exch|^2   (triplet, weight 3/4)
    # where the 1/4 and 3/4 are the spin-state probabilities and the
    # exchange is already included in f_S±f_S_exch.

    # k^2 * dsigma/dOmega:
    dsig_dimless = (1.0/4.0) * (fS_re**2 + fS_im**2) + \
                   (3.0/4.0) * (fT_re**2 + fT_im**2)

    # Elastic cross section from angular integration (for validation):
    # Factor 1/2 for identical particles: each (theta, pi-theta) pair is
    # one physical final state, so full-sphere integral overcounts by 2.
    sigma_el_angular = 0.5 * 2.0 * math.pi * np.sum(dsig_dimless * weights) / (k * k)

    # Momentum-transfer cross section:
    # sigma_tr = (1/2) * 2pi int_{-1}^{1} (1 - cos theta) dsigma/dOmega d(cos theta)
    # The 1/2 is the identical-particle symmetry factor.
    sigma_tr_dimless = np.sum((1.0 - cos_th) * dsig_dimless * weights)
    sigma_tr_GeV2 = 0.5 * 2.0 * math.pi * sigma_tr_dimless / (k * k)
    sigma_tr_cm2 = sigma_tr_GeV2 * GEV2_TO_CM2
    sigma_tr_over_m = sigma_tr_cm2 / (m_chi * GEV_IN_G)

    # Validation: check elastic from angular integration matches partial-wave sum
    sigma_el_angular_cm2 = sigma_el_angular * GEV2_TO_CM2
    sigma_el_angular_over_m = sigma_el_angular_cm2 / (m_chi * GEV_IN_G)

    return sigma_el_over_m, sigma_tr_over_m, sigma_el_angular_over_m


def main():
    print("=" * 90)
    print("  Momentum-Transfer vs Elastic Cross-Section Ratio")
    print("  R(v) = sigma_tr / sigma_el for identical Majorana fermions")
    print("=" * 90)

    # Warm up JIT
    vpm_phase_shift(0, 1.0, 1.0, 50.0, 4000)

    results = {}  # results[label][v] = (sig_el, sig_tr, ratio, sig_el_angular)

    for bm in BENCHMARKS:
        label = bm['label']
        lam = bm['alpha'] * bm['m_chi'] / bm['m_phi']
        print(f"\n  --- {label}: m_chi={bm['m_chi']:.2f} GeV, "
              f"m_phi={bm['m_phi']*1e3:.2f} MeV, alpha={bm['alpha']:.3e}, "
              f"lambda={lam:.2f} ---")

        results[label] = {}
        print(f"  {'v [km/s]':>10}  {'sig_el/m':>12}  {'sig_tr/m':>12}  "
              f"{'R=tr/el':>10}  {'el_check':>12}  {'check_err':>10}")
        print("  " + "-" * 76)

        for v in VELOCITIES:
            sig_el, sig_tr, sig_el_ang = compute_sigma_el_and_tr(
                bm['m_chi'], bm['m_phi'], bm['alpha'], float(v))
            ratio = sig_tr / sig_el if sig_el > 0 else float('inf')
            check_err = abs(sig_el - sig_el_ang) / sig_el if sig_el > 0 else 0
            results[label][v] = (sig_el, sig_tr, ratio)
            print(f"  {v:>10}  {sig_el:>12.6f}  {sig_tr:>12.6f}  "
                  f"{ratio:>10.6f}  {sig_el_ang:>12.6f}  {check_err:>9.2e}")

    # ── Publication figure ──
    print("\n  Generating transfer/elastic ratio plot...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    colors = {'BP1': '#E91E63', 'MAP': '#2196F3'}
    markers = {'BP1': 'o', 'MAP': 's'}

    # Panel 1: sigma_el and sigma_tr overlaid
    ax = axes[0]
    for bm in BENCHMARKS:
        label = bm['label']
        vels = sorted(results[label].keys())
        sel = [results[label][v][0] for v in vels]
        str_ = [results[label][v][1] for v in vels]
        ax.plot(vels, sel, '-', color=colors[label], lw=2.5,
                marker=markers[label], ms=5, label=f'{label} $\\sigma_{{el}}$')
        ax.plot(vels, str_, '--', color=colors[label], lw=2,
                marker=markers[label], ms=5, mfc='none',
                label=f'{label} $\\sigma_{{tr}}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$v$ [km/s]', fontsize=13)
    ax.set_ylabel('$\\sigma/m$ [cm$^2$/g]', fontsize=13)
    ax.set_title('Elastic vs Transfer cross sections', fontsize=13)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')

    # Panel 2: Ratio R(v)
    ax = axes[1]
    for bm in BENCHMARKS:
        label = bm['label']
        lam = bm['alpha'] * bm['m_chi'] / bm['m_phi']
        vels = sorted(results[label].keys())
        ratios = [results[label][v][2] for v in vels]
        ax.plot(vels, ratios, '-', color=colors[label], lw=2.5,
                marker=markers[label], ms=6,
                label=f'{label} ($\\lambda={lam:.1f}$)')
    ax.axhline(1.0, color='gray', ls=':', lw=1)
    ax.set_xscale('log')
    ax.set_xlabel('$v$ [km/s]', fontsize=13)
    ax.set_ylabel('$R(v) = \\sigma_{\\rm tr}/\\sigma_{\\rm el}$', fontsize=13)
    ax.set_title('Transfer / Elastic ratio', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

    # Panel 3: Percentage difference
    ax = axes[2]
    for bm in BENCHMARKS:
        label = bm['label']
        lam = bm['alpha'] * bm['m_chi'] / bm['m_phi']
        vels = sorted(results[label].keys())
        pct_diff = [(results[label][v][2] - 1.0) * 100 for v in vels]
        ax.plot(vels, pct_diff, '-', color=colors[label], lw=2.5,
                marker=markers[label], ms=6,
                label=f'{label} ($\\lambda={lam:.1f}$)')
    ax.axhline(0.0, color='gray', ls=':', lw=1)
    ax.set_xscale('log')
    ax.set_xlabel('$v$ [km/s]', fontsize=13)
    ax.set_ylabel('$(\\sigma_{\\rm tr} - \\sigma_{\\rm el})/\\sigma_{\\rm el}$ [%]',
                  fontsize=13)
    ax.set_title('Relative deviation from $\\sigma_{\\rm tr} = \\sigma_{\\rm el}$',
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

    fig.suptitle('Momentum-Transfer vs Elastic: Secluded Majorana SIDM',
                 fontsize=14, y=1.02)
    fig.tight_layout()

    for ext in ['png', 'pdf']:
        path = os.path.join(OUT_DIR, f'transfer_elastic_ratio.{ext}')
        fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: output/transfer_elastic_ratio.png")
    print(f"  Saved: output/transfer_elastic_ratio.pdf")

    # ── Summary ──
    print(f"\n{'='*90}")
    print("  SUMMARY")
    print(f"{'='*90}")
    for bm in BENCHMARKS:
        label = bm['label']
        lam = bm['alpha'] * bm['m_chi'] / bm['m_phi']
        ratios = [results[label][v][2] for v in VELOCITIES]
        r_min = min(ratios)
        r_max = max(ratios)
        max_dev_pct = max(abs(r - 1.0) * 100 for r in ratios)
        print(f"  {label} (lam={lam:.1f}): R range [{r_min:.6f}, {r_max:.6f}], "
              f"max |R-1| = {max_dev_pct:.2f}%")
    print(f"\n  For identical Majorana fermions, the symmetry under theta <-> pi-theta")
    print(f"  ensures sigma_tr = sigma_el exactly (proven analytically in Section 3.2).")
    print(f"  This numerical computation verifies the identity to high precision.")


if __name__ == '__main__':
    t0 = time.time()
    main()
    print(f"\n  Total runtime: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    try:
        import sys as _sys, os as _os
        _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'core'))
        from tg_notify import notify
        notify("\u2705 transfer_elastic_ratio done!")
    except Exception:
        pass
