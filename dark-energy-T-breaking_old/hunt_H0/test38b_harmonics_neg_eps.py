"""
Test 38b — Higher Harmonics with NEGATIVE ε  (corrected physics)
=================================================================

Test 38 revealed: ε > 0 makes hilltop FLATTER (smaller |V''(π)|),
which REDUCES |w_a|. Wrong direction!

Physics at θ=π:
  V''(π) = −Λ_d⁴/f² (1 − 4ε)
  
  ε < 0: |V''(π)| INCREASES → steeper hilltop → faster rolling → LARGER |w_a|
  ε > −1/4: stability bound (potential stays positive everywhere)

V(θ) = Λ_d⁴ [(1−cosθ) + ε(1−cos 2θ)]  with ε ∈ (−0.25, 0)

At θ=π: V(π) = 2Λ⁴ (independent of ε, since cos2π=1)
Positivity: V > 0 for all θ iff ε > −1/4.

Physical motivation: dilute instanton contributions with alternating signs,
or multi-fermion condensates in the dark sector.
"""

import numpy as np
import sys, os, time

sys.path.insert(0, os.path.dirname(__file__))

import layer8_cosmic_ode as L8
import desi_comparison as DC

from layer8_cosmic_ode import (
    M_PL, H_100_GEV, OMEGA_B_H2, OMEGA_R_H2,
    find_Lambda_d_for_H0, H0_PLANCK_KMS,
)
from desi_comparison import extract_w_of_a, fit_cpl

# ── DESI DR2 targets ────────────────────────────────────────────────────
DESI_TARGETS = {
    'PP':   {'w0': -0.838, 'w0_err': 0.055, 'wa': -0.62,  'wa_err': 0.205},
    'U3':   {'w0': -0.667, 'w0_err': 0.088, 'wa': -1.09,  'wa_err': 0.290},
    'DY5':  {'w0': -0.752, 'w0_err': 0.057, 'wa': -0.86,  'wa_err': 0.215},
}

# ── MAP parameters ──────────────────────────────────────────────────────
M_CHI   = 98.19
ALPHA_D = 3.274e-3
F_REF   = 0.27 * M_PL
OMEGA   = 0.120


def chi2(w0, wa, ds):
    return ((w0 - ds['w0'])/ds['w0_err'])**2 + ((wa - ds['wa'])/ds['wa_err'])**2


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


if __name__ == '__main__':
    V_orig  = L8.V_sigma
    dV_orig = L8.dV_sigma
    V_dc_orig = DC.V_sigma

    print("=" * 78)
    print("  Test 38b — Higher Harmonics: V = Λ⁴[(1-cosθ) + ε(1-cos2θ)]")
    print("  NEGATIVE ε: steeper hilltop → larger |w_a|")
    print("=" * 78)

    # Hilltop analysis
    print("\n  Hilltop analysis at θ=π:  V''(π) = −Λ⁴/f² × (1 − 4ε)")
    print(f"  {'ε':>7}  {'|V″|/Λ⁴f⁻²':>12}  {'vs baseline':>12}")
    print("  " + "-" * 38)
    for eps_show in [0, -0.02, -0.05, -0.08, -0.10, -0.15, -0.20, -0.24]:
        m2_abs = abs(1.0 - 4*eps_show)
        ratio = m2_abs / 1.0
        print(f"  {eps_show:7.2f}  {m2_abs:12.4f}  {ratio:10.1f}×")

    # Scan: negative ε, wider θ_i range
    eps_values = [0.0, -0.02, -0.04, -0.06, -0.08, -0.10, -0.12, -0.15, -0.18, -0.20, -0.22, -0.24]
    theta_values = [2.950, 2.960, 2.970, 2.980, 2.990, 3.000, 3.010, 3.020]

    results = []
    t0 = time.time()

    print(f"\n{'='*95}")
    print(f"  {'ε':>6}  {'θ_i':>5}  {'Λ_d':>11}  {'H₀':>6}  "
          f"{'w₀_CPL':>8}  {'w_a_CPL':>8}  {'w(0)':>8}  "
          f"{'χ²_PP':>6}  {'χ²_U3':>6}  {'χ²_DY5':>6}")
    print("  " + "-" * 90)

    for eps in eps_values:
        install_harmonic_potential(eps)
        for theta_i in theta_values:
            try:
                Ld, res = find_Lambda_d_for_H0(
                    H0_PLANCK_KMS, M_CHI, ALPHA_D, F_REF, theta_i,
                    omega_chi_h2=OMEGA, tol=0.05)
            except Exception:
                continue

            if res is None or res.H0_kms is None or abs(res.H0_kms - 67.4) > 1.0:
                continue

            a_arr, w_arr = extract_w_of_a(res, F_REF, Ld)
            w0_cpl, wa_cpl, _, _ = fit_cpl(a_arr, w_arr, z_max_fit=2.0)

            c2_pp  = chi2(w0_cpl, wa_cpl, DESI_TARGETS['PP'])
            c2_u3  = chi2(w0_cpl, wa_cpl, DESI_TARGETS['U3'])
            c2_dy5 = chi2(w0_cpl, wa_cpl, DESI_TARGETS['DY5'])

            marker = ""
            if c2_pp < 4.0 or c2_dy5 < 4.0:
                marker = " ★"
            if c2_pp < 1.0 or c2_dy5 < 1.0:
                marker = " ★★"

            print(f"  {eps:6.2f}  {theta_i:5.3f}  {Ld:11.4e}  {res.H0_kms:6.2f}  "
                  f"{w0_cpl:8.4f}  {wa_cpl:8.4f}  {res.w_sigma:8.5f}  "
                  f"{c2_pp:6.2f}  {c2_u3:6.2f}  {c2_dy5:6.2f}{marker}")

            results.append({
                'eps': eps, 'theta_i': theta_i, 'Lambda_d': Ld,
                'H0': res.H0_kms, 'w0': w0_cpl, 'wa': wa_cpl,
                'w_today': res.w_sigma, 'Omega_DE': res.Omega_DE,
                'chi2_PP': c2_pp, 'chi2_U3': c2_u3, 'chi2_DY5': c2_dy5,
            })

    elapsed = time.time() - t0

    # Restore
    L8.V_sigma  = V_orig
    L8.dV_sigma = dV_orig
    DC.V_sigma  = V_dc_orig

    print(f"\n{'='*95}")
    print(f"  SCAN COMPLETE ({elapsed:.0f}s, {len(results)} valid points)")
    print(f"{'='*95}")

    if not results:
        print("  No valid points!")
        sys.exit(1)

    # Best per dataset
    for ds_name, ds in DESI_TARGETS.items():
        key = f'chi2_{ds_name}'
        best = min(results, key=lambda r: r[key])
        c2 = best[key]
        sigma = np.sqrt(c2)
        print(f"\n  Best for DESI+CMB+{ds_name}:")
        print(f"    ε = {best['eps']:.2f},  θ_i = {best['theta_i']:.3f}")
        print(f"    w₀ = {best['w0']:.4f}  (DESI: {ds['w0']:.3f} ± {ds['w0_err']:.3f})")
        print(f"    w_a = {best['wa']:.4f}  (DESI: {ds['wa']:.2f} ± {ds['wa_err']:.3f})")
        print(f"    χ² = {c2:.2f}  →  {sigma:.1f}σ combined")
        print(f"    Λ_d = {best['Lambda_d']:.4e} GeV,  H₀ = {best['H0']:.2f}")

    # Global best
    def total_chi2(r):
        return r['chi2_PP'] + r['chi2_U3'] + r['chi2_DY5']

    gb = min(results, key=total_chi2)
    tc = total_chi2(gb)
    print(f"\n  GLOBAL BEST (Σχ²={tc:.1f}):")
    print(f"    ε = {gb['eps']:.2f},  θ_i = {gb['theta_i']:.3f}")
    print(f"    w₀ = {gb['w0']:.4f},  w_a = {gb['wa']:.4f}")
    print(f"    Λ_d = {gb['Lambda_d']:.4e} GeV,  H₀ = {gb['H0']:.2f}")

    # Baseline comparison
    base = [r for r in results if r['eps'] == 0.0]
    if base:
        b0 = min(base, key=total_chi2)
        print(f"\n  BASELINE (ε=0):  w₀ = {b0['w0']:.4f},  w_a = {b0['wa']:.4f},  Σχ² = {total_chi2(b0):.1f}")

    # ── Summary table: best χ²_PP per ε ─────────────────────────────────
    print(f"\n  {'ε':>6}  {'best θ_i':>8}  {'w₀':>8}  {'w_a':>8}  {'χ²_PP':>6}  {'χ²_DY5':>7}")
    print("  " + "-" * 52)
    for eps in eps_values:
        pts = [r for r in results if r['eps'] == eps]
        if pts:
            bp = min(pts, key=lambda r: r['chi2_PP'])
            print(f"  {eps:6.2f}  {bp['theta_i']:8.3f}  {bp['w0']:8.4f}  {bp['wa']:8.4f}  {bp['chi2_PP']:6.2f}  {bp['chi2_DY5']:7.2f}")

    print(f"\n{'='*95}")
    print(f"  Test 38b COMPLETE")
    print(f"{'='*95}")
