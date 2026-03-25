#!/usr/bin/env python3
"""
predictions/majorana_vs_dirac/predict_maj_vs_dir.py
====================================================
Majorana vs Dirac σ_T(v) comparison — the unique fingerprint
of the Secluded Mixed-Majorana SIDM model.

Identical Majorana fermions:
    partial-wave weights  w_even = 1 (singlet),  w_odd = 3 (triplet)

Dirac fermion + same scalar mediator:
    partial-wave weights  w_all  = 1

Same VPM phase shifts δ_l(k), same (m_χ, m_φ, α) — only the statistical
weights differ.  The ratio σ_T^Maj / σ_T^Dir oscillates with velocity,
providing a distinctive discriminant at ≥ 3 velocity scales.

Output
------
- majorana_vs_dirac.png        : 4-panel figure (BP1 abs / MAP abs / BP1 ratio / MAP ratio)
- maj_vs_dir_data.csv          : tabulated σ_T(v) for both statistics, both BPs
"""
import sys, os, math
import numpy as np

_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Import the JIT-compiled VPM primitives from the core solver
from v22_raw_scan import (
    vpm_phase_shift, sigma_T_vpm,
    GEV2_TO_CM2, GEV_IN_G, C_KM_S,
)
from numba import jit

# ── JIT warmup ──
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)


# ══════════════════════════════════════════════════════════════
#  Dirac σ_T  (weight = 1 for ALL partial waves)
# ══════════════════════════════════════════════════════════════
@jit(nopython=True, cache=True)
def sigma_T_dirac(m_chi, m_phi, alpha, v_km_s):
    """Transfer cross section σ_T/m [cm²/g] for distinguishable Dirac fermion.

    All partial-wave weights = 1  (no identical-particle exchange).
    Same VPM phase shifts as Majorana — only weights differ.
    """
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

    l_max = min(max(3, int(kappa) + 3), 80)

    sigma_sum = 0.0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        weight = 1.0                          # <-- Dirac: all ℓ equal
        contrib = weight * (2*l + 1) * math.sin(delta)**2
        sigma_sum += contrib
        if l > int(kappa) + 1 and sigma_sum > 0:
            if contrib / sigma_sum < 1e-3:
                break

    # Dirac: 4π/k² (distinguishable particles — no identical-particle 1/2)
    sigma_GeV2 = 4.0 * math.pi * sigma_sum / (k * k)
    sigma_cm2 = sigma_GeV2 * GEV2_TO_CM2
    return sigma_cm2 / (m_chi * GEV_IN_G)


# ── JIT warmup for Dirac ──
sigma_T_dirac(20.0, 10e-3, 1e-3, 100.0)


# ══════════════════════════════════════════════════════════════
#  Benchmark points — loaded from global_config.json
# ══════════════════════════════════════════════════════════
from global_config import GC
BPS = []
for _lbl in ["BP1", "MAP"]:
    _b = GC.benchmark(_lbl)
    BPS.append({"label": _lbl, "m_chi": _b["m_chi_GeV"], "m_phi_MeV": _b["m_phi_MeV"], "alpha": _b["alpha"]})

# Velocity grid: 5 → 2000 km/s  (dwarfs → clusters)
V_GRID = np.logspace(np.log10(5.0), np.log10(2000.0), 120)

# Astrophysical velocity scales (for annotation)
V_SCALES = {
    r'dSph ($\sim$30)':     30.0,
    r'MW ($\sim$220)':     220.0,
    r'Cluster ($\sim$1200)': 1200.0,
}


def main():
    out_dir = os.path.join(_DIR, 'output')
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("  Majorana vs Dirac σ_T(v) — unique fingerprint")
    print("=" * 80)
    print()

    results = {}

    for bp in BPS:
        label = bp['label']
        m_chi = bp['m_chi']
        m_phi = bp['m_phi_MeV'] / 1000.0      # MeV → GeV
        alpha = bp['alpha']
        lam   = alpha * m_chi / m_phi

        print(f"  {label}: m_χ = {m_chi} GeV, m_φ = {bp['m_phi_MeV']} MeV, "
              f"α = {alpha:.4e}, λ = {lam:.2f}")

        sig_maj = np.zeros(len(V_GRID))
        sig_dir = np.zeros(len(V_GRID))

        for i, v in enumerate(V_GRID):
            sig_maj[i] = sigma_T_vpm(m_chi, m_phi, alpha, v)
            sig_dir[i] = sigma_T_dirac(m_chi, m_phi, alpha, v)

        ratio = np.where(sig_dir > 0, sig_maj / sig_dir, np.nan)

        results[label] = {
            'sig_maj': sig_maj,
            'sig_dir': sig_dir,
            'ratio':   ratio,
            'm_chi':   m_chi,
            'm_phi':   m_phi,
            'alpha':   alpha,
            'lam':     lam,
        }

        # Print key velocity points
        for vname, vval in V_SCALES.items():
            idx = np.argmin(np.abs(V_GRID - vval))
            sm = sig_maj[idx]
            sd = sig_dir[idx]
            r  = ratio[idx]
            print(f"    v = {vval:6.0f} km/s:  σ_Maj = {sm:8.3f}  "
                  f"σ_Dir = {sd:8.3f}  ratio = {r:.4f}")
        print()

    # ── Save CSV ──
    csv_path = os.path.join(out_dir, 'maj_vs_dir_data.csv')
    header = 'v_km_s'
    data_cols = [V_GRID]
    for label in ['BP1', 'MAP']:
        r = results[label]
        header += f',sigma_maj_{label},sigma_dir_{label},ratio_{label}'
        data_cols.extend([r['sig_maj'], r['sig_dir'], r['ratio']])
    np.savetxt(csv_path, np.column_stack(data_cols), delimiter=',',
               header=header, comments='')
    print(f"  Saved: {csv_path}")

    # ══════════════════════════════════════════════════════════
    #  4-panel figure
    # ══════════════════════════════════════════════════════════
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    bp_colors = {'BP1': ('#2196F3', '#E91E63'), 'MAP': ('#4CAF50', '#FF9800')}
    titles = {}
    for _lbl in ['BP1', 'MAP']:
        _b = GC.benchmark(_lbl)
        if _lbl == 'BP1':
            titles[_lbl] = r'BP1: $m_\chi$={mc:.2f} GeV, $m_\phi$={mp:.2f} MeV, $\lambda$={{lam:.2f}}'.format(mc=_b['m_chi_GeV'], mp=_b['m_phi_MeV'])
        else:
            titles[_lbl] = r'MAP: $m_\chi$={mc:.2f} GeV, $m_\phi$={mp:.2f} MeV, $\lambda$={{lam:.1f}}'.format(mc=_b['m_chi_GeV'], mp=_b['m_phi_MeV'])

    for col, label in enumerate(['BP1', 'MAP']):
        r = results[label]
        c_maj, c_dir = bp_colors[label]
        ax_abs = axes[0, col]
        ax_rat = axes[1, col]

        # ── Top: absolute σ_T(v) ──
        ax_abs.loglog(V_GRID, r['sig_maj'], color=c_maj, lw=2.5,
                      label=r'Majorana $w_\ell = (1,3)$')
        ax_abs.loglog(V_GRID, r['sig_dir'], color=c_dir, lw=2.5, ls='--',
                      label=r'Dirac $w_\ell = (1,1)$')

        # Astrophysical bands
        for vname, vval in V_SCALES.items():
            ax_abs.axvline(vval, color='gray', ls=':', lw=0.8, alpha=0.5)
            ax_abs.text(vval*1.05, ax_abs.get_ylim()[0] if ax_abs.get_ylim()[0] > 0 else 0.01,
                        vname, fontsize=7, rotation=90, va='bottom', color='gray')

        ax_abs.set_xlabel(r'$v_{\rm rel}$ [km/s]', fontsize=11)
        ax_abs.set_ylabel(r'$\sigma_T/m_\chi$ [cm$^2$/g]', fontsize=11)
        ax_abs.set_title(titles[label].format(lam=r['lam']), fontsize=11)
        ax_abs.legend(fontsize=10, loc='upper right')
        ax_abs.set_xlim(5, 2000)

        # SIDM constraint band
        ax_abs.axhspan(1, 10, color='green', alpha=0.07)
        ax_abs.axhline(1, color='green', ls='--', lw=0.7, alpha=0.4)
        ax_abs.axhline(10, color='green', ls='--', lw=0.7, alpha=0.4)

        # ── Bottom: ratio ──
        ax_rat.semilogx(V_GRID, r['ratio'], color='k', lw=2.5)

        # Reference lines
        ax_rat.axhline(0.5, color='blue', ls='--', lw=1, alpha=0.5,
                        label=r'Born limit: $\frac{1}{2}$ (even-$\ell$ only)')
        ax_rat.axhline(1.0, color='gray', ls=':', lw=1, alpha=0.5,
                        label='Equal')

        # Color-fill above/below 1
        ax_rat.fill_between(V_GRID, r['ratio'], 1.0,
                            where=r['ratio'] > 1.0,
                            color='red', alpha=0.15, label='Majorana > Dirac')
        ax_rat.fill_between(V_GRID, r['ratio'], 1.0,
                            where=r['ratio'] < 1.0,
                            color='blue', alpha=0.10)

        for vname, vval in V_SCALES.items():
            ax_rat.axvline(vval, color='gray', ls=':', lw=0.8, alpha=0.5)

        ax_rat.set_xlabel(r'$v_{\rm rel}$ [km/s]', fontsize=11)
        ax_rat.set_ylabel(r'$\sigma_T^{\rm Maj} / \sigma_T^{\rm Dir}$', fontsize=11)
        ax_rat.set_title(f'{label}: Majorana / Dirac ratio', fontsize=11)
        ax_rat.legend(fontsize=8, loc='upper right')
        ax_rat.set_xlim(5, 2000)
        ax_rat.set_ylim(0, max(1.5, np.nanmax(r['ratio']) * 1.15))

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'majorana_vs_dirac.png')
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # ══════════════════════════════════════════════════════════
    #  Physics summary
    # ══════════════════════════════════════════════════════════
    print()
    print("  ── Physics summary ──")
    print("  Born regime (low v, s-wave dominant):")
    print("    Only even-ℓ contribute → Majorana = ½ × Dirac (exactly)")
    print()
    print("  Resonant regime (moderate λ):")
    print("    Odd-ℓ resonances get weight 3 → Majorana can EXCEED Dirac")
    print("    Ratio oscillates between ~0.5 and >1 depending on which")
    print("    partial waves are on resonance")
    print()
    print("  Discriminant: measure σ_T at ≥ 3 velocity scales")
    print("  (dSph ~30, MW ~220, cluster ~1200 km/s)")
    print("  Non-monotonic ratio → Majorana signature")
    print()
    print("  Done.")


if __name__ == '__main__':
    main()


if __name__ == '__main__':
    try:
        import sys as _sys, os as _os
        _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'core'))
        from tg_notify import notify
        notify("\u2705 predict_maj_vs_dir done!")
    except Exception:
        pass
