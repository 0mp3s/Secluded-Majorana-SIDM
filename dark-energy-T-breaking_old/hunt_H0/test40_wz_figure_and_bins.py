"""
Test 40 — w(z) Figure for Paper 2  +  DESI DR2 Binned Comparison
==================================================================

Two outputs:
  (A) Publication-quality w(z) figure showing:
      - Model w(z) at best-fit ε=-0.12, θ_i=3.026  (solid)
      - Baseline ε=0 (dashed)
      - ΛCDM w=-1 reference
      - DESI DR2 CPL bands (PP, U3, DY5)
      - DESI BAO effective redshifts marked

  (B) Binned comparison table:
      w_model(z_eff) vs w_DESI(z_eff) at each DESI BAO bin
"""

import numpy as np
import sys, os

sys.path.insert(0, os.path.dirname(__file__))

import layer8_cosmic_ode as L8
import desi_comparison as DC

from layer8_cosmic_ode import (
    M_PL, find_Lambda_d_for_H0, H0_PLANCK_KMS, OMEGA_R_H2, OMEGA_B_H2, _RHO_UNIT,
)
from desi_comparison import extract_w_of_a, fit_cpl

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ── DESI DR2 targets (arXiv:2503.14738v3) ──────────────────────────────
DESI = {
    'PP':  {'w0': -0.838, 'w0_err': 0.055, 'wa': -0.62,  'wa_err': 0.205,
            'label': 'DESI+CMB+PP',     'color': '#1f77b4', 'alpha': 0.18},
    'U3':  {'w0': -0.667, 'w0_err': 0.088, 'wa': -1.09,  'wa_err': 0.290,
            'label': 'DESI+CMB+Union3', 'color': '#ff7f0e', 'alpha': 0.15},
    'DY5': {'w0': -0.752, 'w0_err': 0.057, 'wa': -0.86,  'wa_err': 0.215,
            'label': 'DESI+CMB+DESY5',  'color': '#2ca02c', 'alpha': 0.15},
}

# DESI BAO effective redshifts (DR2, Table I)
Z_EFF = {
    'BGS':       0.295,
    'LRG1':      0.510,
    'LRG2':      0.706,
    'LRG3+ELG1': 0.934,
    'ELG2':      1.317,
    'QSO':       1.491,
    'Ly-α':      2.330,
}

# MAP point
M_CHI   = 98.19
ALPHA_D = 3.274e-3
F_REF   = 0.27 * M_PL
OMEGA   = 0.120


def install_harmonic_potential(eps):
    """Monkey-patch V and dV with higher harmonics."""
    def V_harm(sigma, f, Lambda_d):
        theta = sigma / f
        return Lambda_d**4 * ((1.0 - np.cos(theta)) + eps * (1.0 - np.cos(2.0 * theta)))

    def dV_harm(sigma, f, Lambda_d):
        theta = sigma / f
        return Lambda_d**4 / f * (np.sin(theta) + 2.0 * eps * np.sin(2.0 * theta))

    L8.V_sigma  = V_harm
    L8.dV_sigma = dV_harm
    DC.V_sigma  = V_harm


def restore_potential(V_orig, dV_orig, V_dc_orig):
    L8.V_sigma  = V_orig
    L8.dV_sigma = dV_orig
    DC.V_sigma  = V_dc_orig


def w_cpl(z, w0, wa):
    """CPL parametrization: w(z) = w0 + wa * z/(1+z)."""
    return w0 + wa * z / (1.0 + z)


def sigma_w_cpl(z, w0_err, wa_err):
    """
    Error on w(z) from CPL, assuming w0 and wa uncorrelated (conservative).
    σ_w(z) = sqrt(σ_w0² + (z/(1+z))² σ_wa²)
    """
    f = z / (1.0 + z)
    return np.sqrt(w0_err**2 + f**2 * wa_err**2)


def extract_wz_at_redshifts(res, f, Lambda_d, z_targets):
    """
    Extract w at specific redshift values from ODE solution.
    Returns dict {z: w_value}.
    """
    a_arr, w_arr = extract_w_of_a(res, f, Lambda_d)
    result = {}
    for z in z_targets:
        a_target = 1.0 / (1.0 + z)
        idx = np.argmin(np.abs(a_arr - a_target))
        result[z] = w_arr[idx]
    return result


def run_model(eps, theta_i, label=""):
    """Run Layer 8 ODE and return (a_arr, w_arr, w0_cpl, wa_cpl, Lambda_d)."""
    if eps != 0.0:
        install_harmonic_potential(eps)

    Ld, res = find_Lambda_d_for_H0(
        H0_PLANCK_KMS, M_CHI, ALPHA_D, F_REF, theta_i,
        omega_chi_h2=OMEGA, tol=0.05)

    if res is None or res.H0_kms is None:
        raise RuntimeError(f"Failed for ε={eps}, θ_i={theta_i}")

    a_arr, w_arr = extract_w_of_a(res, F_REF, Ld)
    w0, wa, _, _ = fit_cpl(a_arr, w_arr, z_max_fit=2.0)

    if label:
        print(f"  {label}:")
        print(f"    ε={eps:.3f}, θ_i={theta_i:.4f}, Λ_d={Ld:.4e} GeV, H₀={res.H0_kms:.2f}")
        print(f"    w₀={w0:.5f}, w_a={wa:.5f}, w(z=0)={res.w_sigma:.6f}")

    return a_arr, w_arr, w0, wa, Ld, res


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    V_orig   = L8.V_sigma
    dV_orig  = L8.dV_sigma
    V_dc_orig = DC.V_sigma

    print("=" * 90)
    print("  Test 40 — w(z) Figure + DESI DR2 Binned Comparison")
    print("=" * 90)

    # ── (A) Compute model w(z) ──────────────────────────────────────────
    print("\n  Running models...")

    # Best-fit harmonic (from Test 39)
    a_h, w_h, w0_h, wa_h, Ld_h, res_h = run_model(
        eps=-0.12, theta_i=3.026, label="Best-fit harmonic (ε=-0.12)")

    # Restore and run baseline
    restore_potential(V_orig, dV_orig, V_dc_orig)

    a_b, w_b, w0_b, wa_b, Ld_b, res_b = run_model(
        eps=0.0, theta_i=2.887, label="Baseline (ε=0)")

    # Restore again
    restore_potential(V_orig, dV_orig, V_dc_orig)

    # Convert to z arrays
    z_h = 1.0 / a_h - 1.0
    z_b = 1.0 / a_b - 1.0

    # ── (B) Binned comparison ───────────────────────────────────────────
    print(f"\n{'='*90}")
    print(f"  PART B: Binned w(z) Comparison at DESI BAO Effective Redshifts")
    print(f"{'='*90}")

    z_list = sorted(Z_EFF.values())

    # Get w at those redshifts from harmonic model
    install_harmonic_potential(-0.12)
    Ld_h2, res_h2 = find_Lambda_d_for_H0(
        H0_PLANCK_KMS, M_CHI, ALPHA_D, F_REF, 3.026,
        omega_chi_h2=OMEGA, tol=0.05)
    wz_model = extract_wz_at_redshifts(res_h2, F_REF, Ld_h2, z_list)
    restore_potential(V_orig, dV_orig, V_dc_orig)

    print(f"\n  {'Tracer':>12}  {'z_eff':>5}  {'w_model':>8}  ", end="")
    for ds_name in DESI:
        print(f"{'w_'+ds_name:>9}  {'±':>5}  {'Δσ':>5}  ", end="")
    print()
    print("  " + "-" * 110)

    for tracer, z in sorted(Z_EFF.items(), key=lambda x: x[1]):
        w_mod = wz_model[z]
        print(f"  {tracer:>12}  {z:5.3f}  {w_mod:8.5f}  ", end="")
        for ds_name, ds in DESI.items():
            w_desi = w_cpl(z, ds['w0'], ds['wa'])
            sig_desi = sigma_w_cpl(z, ds['w0_err'], ds['wa_err'])
            tension = abs(w_mod - w_desi) / sig_desi if sig_desi > 0 else 0
            print(f"{w_desi:9.4f}  {sig_desi:5.3f}  {tension:5.2f}  ", end="")
        print()

    # ── Summary per DESI dataset ────────────────────────────────────────
    print(f"\n  Summary (average |Δσ| across z bins):")
    for ds_name, ds in DESI.items():
        tensions = []
        for z in z_list:
            w_mod = wz_model[z]
            w_d = w_cpl(z, ds['w0'], ds['wa'])
            sig_d = sigma_w_cpl(z, ds['w0_err'], ds['wa_err'])
            tensions.append(abs(w_mod - w_d) / sig_d if sig_d > 0 else 0)
        avg_t = np.mean(tensions)
        max_t = np.max(tensions)
        print(f"    {ds['label']:25s}: ⟨|Δσ|⟩ = {avg_t:.2f},  max = {max_t:.2f}")

    # ── CPL comparison table ────────────────────────────────────────────
    print(f"\n  CPL Parameters Comparison:")
    print(f"  {'':20s}  {'w₀':>8}  {'w_a':>8}")
    print(f"  {'-'*40}")
    print(f"  {'Our model (ε=-0.12)':20s}  {w0_h:8.4f}  {wa_h:8.4f}")
    print(f"  {'Baseline (ε=0)':20s}  {w0_b:8.4f}  {wa_b:8.4f}")
    for ds_name, ds in DESI.items():
        print(f"  {'DESI '+ds_name:20s}  {ds['w0']:8.3f}  {ds['wa']:8.2f}")
    print(f"  {'ΛCDM':20s}  {-1.000:8.3f}  {0.00:8.2f}")

    # ── (C) Plot ────────────────────────────────────────────────────────
    print(f"\n  Generating figure...")

    fig, ax = plt.subplots(figsize=(10, 6.5))

    z_plot = np.linspace(0.01, 3.0, 500)

    # DESI DR2 CPL bands
    for ds_name, ds in DESI.items():
        w_central = w_cpl(z_plot, ds['w0'], ds['wa'])
        sig_1 = sigma_w_cpl(z_plot, ds['w0_err'], ds['wa_err'])
        sig_2 = 2.0 * sig_1

        ax.fill_between(z_plot, w_central - sig_2, w_central + sig_2,
                         alpha=ds['alpha'] * 0.5, color=ds['color'], linewidth=0)
        ax.fill_between(z_plot, w_central - sig_1, w_central + sig_1,
                         alpha=ds['alpha'], color=ds['color'], linewidth=0,
                         label=f'{ds["label"]} (1σ/2σ)')
        ax.plot(z_plot, w_central, '--', color=ds['color'], alpha=0.5, linewidth=0.8)

    # ΛCDM
    ax.axhline(y=-1.0, color='gray', linestyle=':', linewidth=1.0, alpha=0.7, label='ΛCDM (w = −1)')

    # Our model: baseline ε=0
    mask_b = (z_b > 0) & (z_b < 3.0)
    sort_b = np.argsort(z_b[mask_b])
    ax.plot(z_b[mask_b][sort_b], w_b[mask_b][sort_b],
            '--', color='gray', linewidth=1.5, alpha=0.7,
            label=f'Baseline ε=0, θ_i=2.887 (w₀={w0_b:.3f}, w_a={wa_b:.3f})')

    # Our model: best-fit harmonic
    mask_h = (z_h > 0) & (z_h < 3.0)
    sort_h = np.argsort(z_h[mask_h])
    ax.plot(z_h[mask_h][sort_h], w_h[mask_h][sort_h],
            '-', color='#d62728', linewidth=2.5,
            label=f'This work: ε=−0.12, θ_i=3.026 (w₀={w0_h:.3f}, w_a={wa_h:.3f})')

    # Mark DESI BAO effective redshifts
    for tracer, z_eff in Z_EFF.items():
        w_at_z = wz_model.get(z_eff, None)
        if w_at_z is not None:
            ax.plot(z_eff, w_at_z, 'o', color='#d62728', markersize=6,
                    markeredgecolor='black', markeredgewidth=0.5, zorder=5)
            ax.annotate(tracer, (z_eff, w_at_z),
                        textcoords="offset points", xytext=(5, 8),
                        fontsize=7, color='#333333', style='italic')

    ax.set_xlabel('Redshift  z', fontsize=13)
    ax.set_ylabel('w(z)', fontsize=13)
    ax.set_title(r'Dark QCD quintessence:  $V(\theta) = \Lambda_d^4[(1-\cos\theta) + \varepsilon(1-\cos 2\theta)]$'
                 '\nvs DESI DR2 (arXiv:2503.14738)', fontsize=12)
    ax.set_xlim(0.0, 2.8)
    ax.set_ylim(-1.8, -0.2)
    ax.legend(loc='lower left', fontsize=8.5, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=11)

    # Inset text box
    txt = (f'MAP point: m_χ=98.2 GeV, α_D=3.27×10⁻³\n'
           f'f = 0.27 M_Pl, Ω_χh²=0.120\n'
           f'Λ_d = {Ld_h:.2e} GeV ≈ 2 meV')
    ax.text(0.97, 0.97, txt, transform=ax.transAxes,
            fontsize=8, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.7))

    plt.tight_layout()
    fig_path = os.path.join(os.path.dirname(__file__), 'test40_wz_desi_dr2.png')
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"  Figure saved: {fig_path}")

    # Also save PDF for paper
    pdf_path = fig_path.replace('.png', '.pdf')
    fig.savefig(pdf_path, bbox_inches='tight')
    print(f"  PDF saved:    {pdf_path}")

    plt.close()
    print(f"\n{'='*90}")
    print(f"  Test 40 COMPLETE")
    print(f"{'='*90}")
