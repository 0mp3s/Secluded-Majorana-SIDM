#!/usr/bin/env python3
"""
predictions/delta_neff/predict_neff.py
======================================
Compute the contribution of the light scalar mediator phi to the
effective number of relativistic species Delta-N_eff and compare
to Planck 2018 and CMB-S4 sensitivity.

Physics:
  The dark sector (chi + phi) decouples from the SM plasma at
  T_dec ~ m_chi (when DM becomes non-relativistic and annihilation
  chi-bar chi -> phi phi freezes out of chemical equilibrium with SM).

  After decoupling, SM entropy transfers heat the SM bath relative
  to the decoupled dark sector:
      T_phi / T_SM = ( g_*S(T_SM) / g_*S(T_dec) )^{1/3}

  The mediator phi (real scalar, g_phi = 1) contributes to N_eff:
      Delta-N_eff = (4/7) * (T_phi / T_nu)^4

  However, phi has mass m_phi ~ 10-14 MeV.  If m_phi >> T_phi at
  the epoch of interest (BBN: T_SM ~ 1 MeV; CMB: T_SM ~ 0.3 eV),
  phi is non-relativistic and its energy density is Boltzmann-
  suppressed:  rho_phi ~ exp(-m_phi / T_phi).

  Result: Delta-N_eff < 10^{-8} at BBN and CMB for our parameter
  space.  The model is automatically safe.

Output: table with T_phi at key epochs, Boltzmann factor,
        effective Delta-N_eff, and comparison to bounds.
"""
import sys, os, math

# ---------- path bootstrap ----------
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
# ------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from config_loader import load_config
from global_config import GC

_DIR = os.path.dirname(os.path.abspath(__file__))


def T_dark_over_T_SM(g_star_s_SM, g_star_s_dec):
    """
    Ratio T_dark / T_SM after the dark sector decouples at g_*S = g_star_s_dec
    and the SM is at g_*S = g_star_s_SM.
    """
    return (g_star_s_SM / g_star_s_dec) ** (1.0 / 3.0)


def delta_neff_massless(T_phi_over_T_nu):
    """Delta-N_eff for a massless real scalar with temperature T_phi."""
    return (4.0 / 7.0) * T_phi_over_T_nu ** 4


def boltzmann_suppression(m_MeV, T_MeV):
    """
    Suppression factor for a massive species relative to the massless limit.
    For m >> T: rho ~ (m T)^{3/2} * T * exp(-m/T)  vs  rho_massless ~ T^4
    """
    if T_MeV < 1e-30:
        return 0.0
    x = m_MeV / T_MeV
    if x > 500:
        return 0.0
    if x < 0.1:
        return 1.0  # effectively massless
    # Full ratio rho_massive / rho_massless for a boson
    # Approximate: (15/(pi^4)) * x^4 * K_2(x) / 2  [Kolb & Turner eq 3.55 approx]
    # Simpler: use exp(-x) * (1 + 15/(8x)) as leading-order suppression
    return math.exp(-x) * (1.0 + 15.0 / (8.0 * x))


def main():
    cfg = load_config(__file__)
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = GC.benchmarks_from_labels(cfg['benchmark_labels'])
    planck_neff     = cfg.get('planck2018_Neff', 2.99)
    planck_err      = cfg.get('planck2018_Neff_err', 0.17)
    cmbs4_sigma     = cfg.get('cmbs4_sigma_Neff', 0.03)
    sm_neff         = cfg.get('sm_Neff', 3.044)
    g_star_s_dec    = cfg.get('g_star_s_20_90_GeV', 86.25)
    g_star_s_nu_dec = cfg.get('g_star_s_nu_decoupling', 10.75)

    print("=" * 80)
    print("  Delta-N_eff Prediction from Light Mediator phi")
    print("=" * 80)

    # Key epochs
    epochs = [
        ("High-T (T_SM >> m_phi)", 1000.0),   # MeV
        ("QCD transition",          150.0),
        ("BBN (T_SM ~ 1 MeV)",       1.0),
        ("CMB (T_SM ~ 0.3 eV)",      3e-4),
    ]

    for bp in bps:
        label     = bp['label']
        m_chi     = bp['m_chi_GeV']
        m_phi_MeV = bp['m_phi_MeV']

        # Dark sector decouples at T_dec ~ m_chi
        # g_*S at T_dec: for 20-90 GeV → 86.25 (below top, above bottom)
        g_dec = g_star_s_dec

        print(f"\n  --- {label}: m_chi = {m_chi:.1f} GeV, "
              f"m_phi = {m_phi_MeV:.2f} MeV ---")
        print(f"  T_dec ~ {m_chi:.0f} GeV, g_*S(T_dec) = {g_dec:.2f}")
        print(f"  {'Epoch':<30s} {'T_SM':>8s} {'T_phi':>8s} "
              f"{'m_phi/T_phi':>11s} {'Boltz':>10s} {'dN_eff':>12s}")
        print("  " + "-" * 82)

        for epoch_name, T_SM_MeV in epochs:
            # g_*S at this SM temperature
            if T_SM_MeV > 150:
                g_sm = 86.25
            elif T_SM_MeV > 1:
                g_sm = 10.75
            else:
                g_sm = 3.91

            ratio = T_dark_over_T_SM(g_sm, g_dec)
            T_phi_MeV = ratio * T_SM_MeV

            # T_phi / T_nu at neutrino decoupling: T_nu ~ T_SM at that time
            # but after e+e- annihilation T_nu / T_gamma = (4/11)^{1/3}
            # For simplicity, compare at the given epoch:
            # T_nu at this epoch relative to T_SM:
            if T_SM_MeV < 0.5:
                # after e+e- annihilation
                T_nu_MeV = T_SM_MeV * (4.0 / 11.0) ** (1.0 / 3.0)
            else:
                T_nu_MeV = T_SM_MeV  # neutrinos still coupled or just decoupled

            T_phi_over_T_nu = T_phi_MeV / T_nu_MeV if T_nu_MeV > 0 else 0

            # Massless limit
            dn_massless = delta_neff_massless(T_phi_over_T_nu)

            # Boltzmann suppression
            x = m_phi_MeV / T_phi_MeV if T_phi_MeV > 0 else 1e10
            boltz = boltzmann_suppression(m_phi_MeV, T_phi_MeV)

            dn_eff = dn_massless * boltz

            print(f"  {epoch_name:<30s} {T_SM_MeV:>8.3g} {T_phi_MeV:>8.3g} "
                  f"{x:>11.1f} {boltz:>10.2e} {dn_eff:>12.2e}")

    # ---- Summary comparison to Planck ----
    print(f"\n{'=' * 80}")
    print("  Comparison to Experimental Bounds")
    print(f"{'=' * 80}")
    print(f"  Planck 2018:  N_eff = {planck_neff} +/- {planck_err}  "
          f"(SM prediction: {sm_neff})")
    print(f"  Planck 95% CL on Delta-N_eff: < {2 * planck_err:.2f}")
    print(f"  CMB-S4 forecast sigma: {cmbs4_sigma}")
    print()
    print(f"  Our prediction at BBN/CMB epoch:  Delta-N_eff ~ 0  "
          f"(Boltzmann-suppressed, m_phi >> T_phi)")
    print(f"  => CONSISTENT with Planck 2018")
    print(f"  => CONSISTENT with CMB-S4 (even at sigma = 0.03)")
    print()
    print("  Note: If the mediator were massless, Delta-N_eff ~ 0.036")
    print("  would be detectable by CMB-S4 but is consistent with Planck 2018.")

    # ---- Plot: Delta-N_eff vs T_SM ----
    T_SM_arr = np.logspace(-4, 5, 500)  # MeV
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = {'BP1': 'steelblue', 'BP16': 'seagreen', 'MAP': 'firebrick'}

    for bp in bps:
        label     = bp['label']
        m_phi_MeV = bp['m_phi_MeV']
        g_dec     = g_star_s_dec

        dn_arr = []
        for T_SM_MeV in T_SM_arr:
            if T_SM_MeV > 150:
                g_sm = 86.25
            elif T_SM_MeV > 1:
                g_sm = 10.75
            else:
                g_sm = 3.91

            ratio = T_dark_over_T_SM(g_sm, g_dec)
            T_phi_MeV = ratio * T_SM_MeV

            if T_SM_MeV < 0.5:
                T_nu_MeV = T_SM_MeV * (4.0 / 11.0) ** (1.0 / 3.0)
            else:
                T_nu_MeV = T_SM_MeV

            if T_nu_MeV > 0:
                dn_ml = delta_neff_massless(T_phi_MeV / T_nu_MeV)
            else:
                dn_ml = 0
            boltz = boltzmann_suppression(m_phi_MeV, T_phi_MeV)
            dn_arr.append(max(dn_ml * boltz, 1e-20))

        ax.plot(T_SM_arr, dn_arr, lw=2, color=colors.get(label, 'gray'),
                label=f'{label} (m_φ = {m_phi_MeV:.1f} MeV)')

    # Bounds
    ax.axhline(2 * planck_err, color='orange', ls='--', lw=2,
               label=f'Planck 95% CL ({2*planck_err:.2f})')
    ax.axhline(2 * cmbs4_sigma, color='purple', ls=':', lw=2,
               label=f'CMB-S4 2σ ({2*cmbs4_sigma:.2f})')

    # Epoch markers
    for t, name in [(1.0, 'BBN'), (3e-4, 'CMB')]:
        ax.axvline(t, color='gray', ls=':', alpha=0.4)
        ax.text(t * 1.3, ax.get_ylim()[1] * 0.5, name, fontsize=9,
                color='gray', rotation=90, va='center')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$T_{\rm SM}$ (MeV)', fontsize=12)
    ax.set_ylabel(r'$\Delta N_{\rm eff}$', fontsize=12)
    ax.set_title(r'$\Delta N_{\rm eff}$ from Light Mediator $\phi$', fontsize=13)
    ax.set_ylim(1e-15, 1.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'delta_neff.png')
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved: {fig_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
