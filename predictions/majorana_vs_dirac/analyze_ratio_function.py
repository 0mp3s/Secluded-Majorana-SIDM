#!/usr/bin/env python3
"""
predictions/majorana_vs_dirac/analyze_ratio_function.py
=======================================================
Full calculus-style analysis of R(v) = sigma_T^Maj / sigma_T^Dir.

Treats R(v) as a proper function of velocity and extracts:
  1. Domain & range
  2. Monotonicity regions (dR/dv > 0 vs < 0)
  3. All local extrema (maxima, minima) with locations
  4. Inflection points (d2R/dv2 = 0)
  5. Asymptotic behavior: v->0 (Born) and v->inf (semi-classical)
  6. Partial-wave decomposition: which l contributes at each v
  7. Comparison across lambda values (BP1, MAP, hypothetical old-MAP)

Output
------
- r_v_analysis.png  : multi-panel calculus analysis figure
- r_v_analysis.csv  : R(v), R'(v), R''(v) tabulated
"""
import sys, os, math
import numpy as np

_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from v22_raw_scan import vpm_phase_shift, sigma_T_vpm, GEV2_TO_CM2, GEV_IN_G, C_KM_S
from numba import jit

# JIT warmup
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)


@jit(nopython=True, cache=True)
def sigma_T_dirac(m_chi, m_phi, alpha, v_km_s):
    """Dirac sigma_T/m [cm2/g] -- all weights = 1, factor 4pi/k2."""
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
    return sigma_GeV2 * GEV2_TO_CM2 / (m_chi * GEV_IN_G)


@jit(nopython=True, cache=True)
def partial_wave_decomposition(m_chi, m_phi, alpha, v_km_s):
    """Return fractions of (2l+1)sin2(delta_l) for even/odd l."""
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15:
        return 1.0, 0.0, 0
    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000
    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)

    sum_even = 0.0
    sum_odd  = 0.0
    total    = 0.0
    n_ell    = 0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        contrib = (2*l + 1) * math.sin(delta)**2
        total += contrib
        if l % 2 == 0:
            sum_even += contrib
        else:
            sum_odd += contrib
        n_ell = l
        if l > int(kappa) + 1 and total > 0:
            if contrib / total < 1e-3:
                break
    if total > 0:
        return sum_even / total, sum_odd / total, n_ell
    return 1.0, 0.0, n_ell


# JIT warmup
sigma_T_dirac(20.0, 10e-3, 1e-3, 100.0)
partial_wave_decomposition(20.0, 10e-3, 1e-3, 100.0)


# ---- Benchmarks — loaded from global_config.json ----
from global_config import GC
BENCHMARKS = []
for _lbl, _clr in [("BP1", "#2196F3"), ("MAP", "#4CAF50")]:
    _b = GC.benchmark(_lbl)
    BENCHMARKS.append({"label": _lbl, "m_chi": _b["m_chi_GeV"], "m_phi_MeV": _b["m_phi_MeV"], "alpha": _b["alpha"], "color": _clr})
# old-MAP is a hypothetical comparison point, not a global benchmark
BENCHMARKS.append({"label": "old-MAP", "m_chi": 90.64, "m_phi_MeV": 13.85, "alpha": 2.546e-2, "color": "#E91E63"})

V_GRID = np.logspace(np.log10(2.0), np.log10(3000.0), 500)

V_SCALES = {'dSph\n30': 30.0, 'MW\n220': 220.0, 'Cluster\n1200': 1200.0}


def compute_ratio(bp):
    m_chi = bp['m_chi']
    m_phi = bp['m_phi_MeV'] / 1000.0
    alpha = bp['alpha']
    lam = alpha * m_chi / m_phi
    sig_maj = np.zeros(len(V_GRID))
    sig_dir = np.zeros(len(V_GRID))
    for i, v in enumerate(V_GRID):
        sig_maj[i] = sigma_T_vpm(m_chi, m_phi, alpha, v)
        sig_dir[i] = sigma_T_dirac(m_chi, m_phi, alpha, v)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(sig_dir > 1e-20, sig_maj / sig_dir, np.nan)
    return sig_maj, sig_dir, ratio, lam


def compute_pw_fractions(bp):
    m_chi = bp['m_chi']
    m_phi = bp['m_phi_MeV'] / 1000.0
    alpha = bp['alpha']
    f_even = np.zeros(len(V_GRID))
    f_odd  = np.zeros(len(V_GRID))
    n_ell  = np.zeros(len(V_GRID), dtype=int)
    for i, v in enumerate(V_GRID):
        fe, fo, nl = partial_wave_decomposition(m_chi, m_phi, alpha, v)
        f_even[i] = fe
        f_odd[i]  = fo
        n_ell[i]  = nl
    return f_even, f_odd, n_ell


def numerical_derivatives(v, R):
    log_v = np.log10(v)
    dR_dlogv = np.gradient(R, log_v)
    d2R_dlogv2 = np.gradient(dR_dlogv, log_v)
    return dR_dlogv, d2R_dlogv2


def find_extrema_and_inflections(v, R, dR, d2R):
    extrema_max, extrema_min, inflections = [], [], []
    for i in range(1, len(v) - 1):
        if np.isnan(dR[i-1]) or np.isnan(dR[i+1]):
            continue
        if dR[i-1] > 0 and dR[i+1] < 0:
            extrema_max.append((v[i], R[i]))
        elif dR[i-1] < 0 and dR[i+1] > 0:
            extrema_min.append((v[i], R[i]))
        if np.isnan(d2R[i-1]) or np.isnan(d2R[i+1]):
            continue
        if d2R[i-1] * d2R[i+1] < 0:
            inflections.append((v[i], R[i]))
    return extrema_max, extrema_min, inflections


def main():
    out_dir = os.path.join(_DIR, 'output')
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("  R(v) = sigma_T^Maj / sigma_T^Dir  --  Full Calculus Analysis")
    print("=" * 80)

    all_results = {}

    for bp in BENCHMARKS:
        label = bp['label']
        sig_maj, sig_dir, ratio, lam = compute_ratio(bp)
        dR, d2R = numerical_derivatives(V_GRID, ratio)
        emax, emin, infl = find_extrema_and_inflections(V_GRID, ratio, dR, d2R)
        f_even, f_odd, n_ell = compute_pw_fractions(bp)

        all_results[label] = {
            'ratio': ratio, 'lam': lam, 'color': bp['color'],
            'dR': dR, 'd2R': d2R,
            'emax': emax, 'emin': emin, 'infl': infl,
            'f_even': f_even, 'f_odd': f_odd, 'n_ell': n_ell,
            'sig_maj': sig_maj, 'sig_dir': sig_dir,
        }

        print(f"\n  -- {label} (lambda = {lam:.2f}) --")
        valid = ratio[~np.isnan(ratio)]
        print(f"    R(v_min) = {valid[0]:.4f}   (Born limit: 0.5)")
        print(f"    R(v_max) = {valid[-1]:.4f}   (semi-classical: ~1)")
        print(f"    Range:     [{np.nanmin(ratio):.4f}, {np.nanmax(ratio):.4f}]")

        if emax:
            for v_ex, r_ex in emax:
                exceeds = " ** R > 1 **" if r_ex > 1.0 else ""
                print(f"    LOCAL MAX:  R = {r_ex:.4f} at v = {v_ex:.1f} km/s{exceeds}")
        if emin:
            for v_ex, r_ex in emin:
                print(f"    LOCAL MIN:  R = {r_ex:.4f} at v = {v_ex:.1f} km/s")
        if not emax and not emin:
            print(f"    MONOTONIC (no local extrema)")

        print(f"    Inflection points: {len(infl)}")
        for v_inf, r_inf in infl[:5]:
            print(f"      v = {v_inf:.1f} km/s, R = {r_inf:.4f}")

        # Partial wave at key velocities
        idx_30  = np.argmin(np.abs(V_GRID - 30))
        idx_220 = np.argmin(np.abs(V_GRID - 220))
        idx_1200 = np.argmin(np.abs(V_GRID - 1200))
        print(f"    Partial wave fractions (even / odd):")
        for name, idx in [('30', idx_30), ('220', idx_220), ('1200', idx_1200)]:
            fe, fo = f_even[idx], f_odd[idx]
            R_check = (fe + 3*fo) / 2.0
            print(f"      v={name:>4s}: even={fe:.3f}, odd={fo:.3f}, "
                  f"n_ell={n_ell[idx]:2d}, R_formula={R_check:.4f}")

        # Where does f_odd peak?
        peak_idx = np.argmax(f_odd)
        print(f"    Max f_odd = {f_odd[peak_idx]:.4f} at v = {V_GRID[peak_idx]:.1f} km/s")
        if f_odd[peak_idx] > 0.5:
            print(f"    ** f_odd > 1/2 => R > 1 is possible! **")
        else:
            print(f"    f_odd < 1/2 everywhere => R < 1 always")

    # ================================================================
    #  6-panel figure
    # ================================================================
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.25)

    # Panel 1: R(v) for all lambda
    ax1 = fig.add_subplot(gs[0, 0])
    for label, r in all_results.items():
        lbl = f"{label} ($\\lambda$={r['lam']:.1f})"
        ax1.semilogx(V_GRID, r['ratio'], color=r['color'], lw=2.5, label=lbl)
        for v_ex, r_ex in r['emax']:
            ax1.plot(v_ex, r_ex, 'v', color=r['color'], ms=10, zorder=5)
        for v_ex, r_ex in r['emin']:
            ax1.plot(v_ex, r_ex, '^', color=r['color'], ms=10, zorder=5)
    ax1.axhline(0.5, color='blue', ls='--', lw=1, alpha=0.4, label='Born limit (1/2)')
    ax1.axhline(1.0, color='gray', ls=':', lw=1, alpha=0.5, label='R = 1')
    for vn, vv in V_SCALES.items():
        ax1.axvline(vv, color='gray', ls=':', lw=0.7, alpha=0.3)
        ax1.text(vv, 0.37, vn, fontsize=7, ha='center', color='gray')
    ax1.set_xlabel(r'$v_{\rm rel}$ [km/s]')
    ax1.set_ylabel(r'$R(v) = \sigma_T^{\rm Maj} / \sigma_T^{\rm Dir}$')
    ax1.set_title(r'$R(v)$ across $\lambda$ regimes', fontsize=12)
    ax1.legend(fontsize=8, loc='lower right')
    ax1.set_xlim(2, 3000)
    ax1.set_ylim(0.35, max(1.35, max(np.nanmax(r['ratio']) for r in all_results.values()) * 1.1))

    # Panel 2: dR/d(log v)
    ax2 = fig.add_subplot(gs[0, 1])
    for label, r in all_results.items():
        ax2.semilogx(V_GRID, r['dR'], color=r['color'], lw=2.0, label=label)
    ax2.axhline(0, color='gray', ls='-', lw=0.5)
    for vn, vv in V_SCALES.items():
        ax2.axvline(vv, color='gray', ls=':', lw=0.7, alpha=0.3)
    ax2.set_xlabel(r'$v_{\rm rel}$ [km/s]')
    ax2.set_ylabel(r'$dR/d(\log v)$')
    ax2.set_title(r"$R'(v)$ -- velocity sensitivity", fontsize=12)
    ax2.legend(fontsize=8)
    ax2.set_xlim(2, 3000)

    # Panel 3: d2R/d(logv)2
    ax3 = fig.add_subplot(gs[1, 0])
    for label, r in all_results.items():
        ax3.semilogx(V_GRID, r['d2R'], color=r['color'], lw=2.0, label=label)
        for v_inf, r_inf in r['infl'][:6]:
            ax3.plot(v_inf, 0, 'o', color=r['color'], ms=6, zorder=5)
    ax3.axhline(0, color='gray', ls='-', lw=0.5)
    for vn, vv in V_SCALES.items():
        ax3.axvline(vv, color='gray', ls=':', lw=0.7, alpha=0.3)
    ax3.set_xlabel(r'$v_{\rm rel}$ [km/s]')
    ax3.set_ylabel(r"$d^2R/d(\log v)^2$")
    ax3.set_title(r"$R''(v)$ -- curvature (inflection points marked)", fontsize=12)
    ax3.legend(fontsize=8)
    ax3.set_xlim(2, 3000)

    # Panel 4: Partial wave fractions - MAP
    ax4 = fig.add_subplot(gs[1, 1])
    r_map = all_results['MAP']
    ax4.semilogx(V_GRID, r_map['f_even'], color='#2196F3', lw=2.5,
                 label=r'Even $\ell$ (weight 1)')
    ax4.semilogx(V_GRID, r_map['f_odd'], color='#E91E63', lw=2.5,
                 label=r'Odd $\ell$ (weight 3)')
    ax4.fill_between(V_GRID, r_map['f_even'], alpha=0.15, color='#2196F3')
    ax4.fill_between(V_GRID, r_map['f_odd'], alpha=0.15, color='#E91E63')
    ax4.axhline(0.5, color='gray', ls='--', lw=1, alpha=0.5,
                label=r'$f_{\rm odd}=1/2$ (R=1 threshold)')
    for vn, vv in V_SCALES.items():
        ax4.axvline(vv, color='gray', ls=':', lw=0.7, alpha=0.3)
    ax4.set_xlabel(r'$v_{\rm rel}$ [km/s]')
    ax4.set_ylabel('Fraction of total cross section')
    ax4.set_title(r'Partial-wave decomposition -- MAP ($\lambda=48.6$)', fontsize=12)
    ax4.legend(fontsize=8)
    ax4.set_xlim(2, 3000); ax4.set_ylim(0, 1.05)

    # Panel 5: R(v) anatomy - decomposed
    ax5 = fig.add_subplot(gs[2, 0])
    for label in ['MAP', 'old-MAP']:
        r = all_results[label]
        term_even = r['f_even'] / 2.0
        term_odd  = 3 * r['f_odd'] / 2.0
        ls = '-' if label == 'MAP' else '--'
        ax5.semilogx(V_GRID, term_even, color='#2196F3', lw=2, ls=ls,
                     label=f'$f_{{\\rm even}}/2$ ({label})')
        ax5.semilogx(V_GRID, term_odd, color='#E91E63', lw=2, ls=ls,
                     label=f'$3f_{{\\rm odd}}/2$ ({label})')
    ax5.axhline(0.5, color='blue', ls=':', lw=1, alpha=0.3)
    ax5.axhline(1.0, color='gray', ls=':', lw=1, alpha=0.3)
    for vn, vv in V_SCALES.items():
        ax5.axvline(vv, color='gray', ls=':', lw=0.7, alpha=0.3)
    ax5.set_xlabel(r'$v_{\rm rel}$ [km/s]')
    ax5.set_ylabel('Contribution to R')
    ax5.set_title(r'Anatomy: $R = f_{\rm even}/2 + 3f_{\rm odd}/2$', fontsize=12)
    ax5.legend(fontsize=7, loc='center right', ncol=2)
    ax5.set_xlim(2, 3000)

    # Panel 6: n_ell vs velocity
    ax6 = fig.add_subplot(gs[2, 1])
    for label, r in all_results.items():
        ax6.semilogx(V_GRID, r['n_ell'], color=r['color'], lw=2.0, label=label)
    for vn, vv in V_SCALES.items():
        ax6.axvline(vv, color='gray', ls=':', lw=0.7, alpha=0.3)
    ax6.set_xlabel(r'$v_{\rm rel}$ [km/s]')
    ax6.set_ylabel(r'$\ell_{\rm max}$ converged')
    ax6.set_title('Number of contributing partial waves', fontsize=12)
    ax6.legend(fontsize=8)
    ax6.set_xlim(2, 3000)

    fig_path = os.path.join(out_dir, 'r_v_analysis.png')
    fig.savefig(fig_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved: {fig_path}")

    # ---- Analytical insight summary ----
    print("\n" + "=" * 80)
    print("  ANALYTICAL INSIGHTS")
    print("=" * 80)
    print("""
    R(v) = sigma_T^Maj / sigma_T^Dir = (f_even + 3*f_odd) / 2

    where f_even + f_odd = 1 are the fractions of the total
    sum((2l+1)*sin2(delta_l)) from even vs odd partial waves.

    LIMITING CASES:
    -------------------------------------------------------
    v -> 0 (Born):  Only l=0 (even) => f_even=1, f_odd=0
                    R = (1 + 0)/2 = 1/2 exactly

    v -> inf (semi-classical): Many l, even/odd roughly equal
                    f_even ~ f_odd ~ 1/2
                    R = (1/2 + 3/2)/2 = 1 exactly

    CRITICAL CONDITION for R > 1:
    -------------------------------------------------------
    R > 1  <=>  f_odd > 1/2

    Odd-l must contribute MORE than half the total cross section.
    This requires specific odd-l resonances to dominate, 
    which happens in the deep resonant regime (large lambda).

    MAP (lambda=48.6): f_odd never exceeds 1/2 => R < 1 everywhere
    old-MAP (lambda=167): f_odd CAN exceed 1/2 => R > 1 possible
    """)

    # CSV
    csv_path = os.path.join(out_dir, 'r_v_analysis.csv')
    cols = [V_GRID]
    header = 'v_km_s'
    for label, r in all_results.items():
        s = label.replace(' ', '_').replace('-', '_')
        header += f',R_{s},dRdlogv_{s},d2Rdlogv2_{s},f_even_{s},f_odd_{s}'
        cols.extend([r['ratio'], r['dR'], r['d2R'], r['f_even'], r['f_odd']])
    np.savetxt(csv_path, np.column_stack(cols), delimiter=',',
               header=header, comments='')
    print(f"\n  Saved: {csv_path}")


if __name__ == '__main__':
    main()
    try:
        from tg_notify import notify
        notify("\u2705 analyze_ratio_function done!")
    except Exception:
        pass
