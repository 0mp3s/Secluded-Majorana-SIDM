"""
Test 38 — Higher Harmonics: Can ε(1−cos 2θ) fix w_a?
=====================================================

Test 37 found w_a ≈ −0.18 vs DESI DR2 wants −0.6 to −1.1.
The cosine potential is too flat near the hilltop (θ≈π).

Dark QCD instanton corrections naturally generate higher harmonics:
  V(θ) = Λ_d⁴ [ (1−cosθ) + ε(1−cos 2θ) ]

At the hilltop (θ=π):  V''(π) = −Λ_d⁴/f² (1 − 4ε)
  → ε > 0: hilltop STEEPER → faster rolling → larger |w_a|
  → ε > 1/4: hilltop becomes local minimum (excluded)

Strategy: Scan ε from 0 to 0.24 at θ_i=2.97, 2.98, 2.99, 3.00.
For each (θ_i, ε), find Λ_d → H₀=67.4, extract CPL, compare DESI DR2.
"""

import numpy as np
import sys, os, time

sys.path.insert(0, os.path.dirname(__file__))

# Import the module (not just functions) so we can monkey-patch
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
    """
    Monkey-patch V_sigma and dV_sigma in both layer8_cosmic_ode and
    desi_comparison to use V = Λ⁴[(1-cosθ) + ε(1-cos2θ)].
    """
    def V_harm(sigma, f, Lambda_d):
        theta = sigma / f
        return Lambda_d**4 * ((1.0 - np.cos(theta)) + eps * (1.0 - np.cos(2.0 * theta)))

    def dV_harm(sigma, f, Lambda_d):
        theta = sigma / f
        return Lambda_d**4 / f * (np.sin(theta) + 2.0 * eps * np.sin(2.0 * theta))

    # Patch both modules
    L8.V_sigma  = V_harm
    L8.dV_sigma = dV_harm
    DC.V_sigma  = V_harm


# ══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    # Save originals
    V_orig  = L8.V_sigma
    dV_orig = L8.dV_sigma
    V_dc_orig = DC.V_sigma

    print("=" * 78)
    print("  Test 38 — Higher Harmonics: V = Λ⁴[(1-cosθ) + ε(1-cos2θ)]")
    print("  Scanning ε to find w_a match with DESI DR2")
    print("=" * 78)

    # Physics at hilltop:
    print("\n  Hilltop analysis at θ=π:")
    print(f"  {'ε':>6}  {'V_max/Λ⁴':>10}  {'m²_eff/Λ⁴f⁻²':>14}  {'steeper?':>10}")
    print("  " + "-" * 48)
    for eps_show in [0, 0.05, 0.10, 0.15, 0.20, 0.24]:
        V_max = 2.0 + 0  # (1-cos π) + ε(1-cos 2π) = 2 + 0
        m2_eff = -(1.0 - 4*eps_show)  # V''(π)/Λ⁴f⁻²
        steep = "STEEPER" if eps_show > 0 else "baseline"
        if eps_show >= 0.25:
            steep = "FLAT (min!)"
        print(f"  {eps_show:6.2f}  {V_max:10.4f}  {m2_eff:14.4f}  {steep:>10}")

    # Scan
    eps_values = [0.0, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.24]
    theta_values = [2.970, 2.980, 2.990, 3.000]

    results = []
    t0 = time.time()

    print(f"\n{'='*78}")
    print(f"  {'ε':>5}  {'θ_i':>5}  {'Λ_d':>11}  {'H₀':>6}  "
          f"{'w₀_CPL':>8}  {'w_a_CPL':>8}  {'w(0)':>8}  "
          f"{'χ²_PP':>6}  {'χ²_U3':>6}  {'χ²_DY5':>6}")
    print("  " + "-" * 95)

    for eps in eps_values:
        install_harmonic_potential(eps)
        for theta_i in theta_values:
            try:
                Ld, res = find_Lambda_d_for_H0(
                    H0_PLANCK_KMS, M_CHI, ALPHA_D, F_REF, theta_i,
                    omega_chi_h2=OMEGA, tol=0.05)
            except Exception:
                print(f"  {eps:5.2f}  {theta_i:5.3f}  FAILED (solver)")
                continue

            if res is None or res.H0_kms is None or abs(res.H0_kms - 67.4) > 1.0:
                H0_show = res.H0_kms if res and res.H0_kms else 0
                print(f"  {eps:5.2f}  {theta_i:5.3f}  FAILED (H₀={H0_show:.1f})")
                continue

            a_arr, w_arr = extract_w_of_a(res, F_REF, Ld)
            w0_cpl, wa_cpl, _, _ = fit_cpl(a_arr, w_arr, z_max_fit=2.0)

            c2_pp  = chi2(w0_cpl, wa_cpl, DESI_TARGETS['PP'])
            c2_u3  = chi2(w0_cpl, wa_cpl, DESI_TARGETS['U3'])
            c2_dy5 = chi2(w0_cpl, wa_cpl, DESI_TARGETS['DY5'])

            marker = ""
            if c2_dy5 < 4.0:  # within 2σ combined
                marker = " ★"
            elif c2_pp < 4.0:
                marker = " ★"

            print(f"  {eps:5.2f}  {theta_i:5.3f}  {Ld:11.4e}  {res.H0_kms:6.2f}  "
                  f"{w0_cpl:8.4f}  {wa_cpl:8.4f}  {res.w_sigma:8.5f}  "
                  f"{c2_pp:6.2f}  {c2_u3:6.2f}  {c2_dy5:6.2f}{marker}")

            results.append({
                'eps': eps, 'theta_i': theta_i, 'Lambda_d': Ld,
                'H0': res.H0_kms, 'w0': w0_cpl, 'wa': wa_cpl,
                'w_today': res.w_sigma, 'Omega_DE': res.Omega_DE,
                'chi2_PP': c2_pp, 'chi2_U3': c2_u3, 'chi2_DY5': c2_dy5,
            })

    elapsed = time.time() - t0

    # Restore originals
    L8.V_sigma  = V_orig
    L8.dV_sigma = dV_orig
    DC.V_sigma  = V_dc_orig

    print(f"\n{'='*78}")
    print(f"  SCAN COMPLETE ({elapsed:.0f}s, {len(results)} valid points)")
    print(f"{'='*78}")

    if not results:
        print("  No valid points!")
        sys.exit(1)

    # Best fit per dataset
    for ds_name, ds in DESI_TARGETS.items():
        key = f'chi2_{ds_name}'
        best = min(results, key=lambda r: r[key])
        c2 = best[key]
        print(f"\n  Best for DESI+CMB+{ds_name}:")
        print(f"    ε = {best['eps']:.2f},  θ_i = {best['theta_i']:.3f}")
        print(f"    w₀ = {best['w0']:.4f}  (DESI: {ds['w0']:.3f} ± {ds['w0_err']:.3f})")
        print(f"    w_a = {best['wa']:.4f}  (DESI: {ds['wa']:.2f} ± {ds['wa_err']:.3f})")
        print(f"    χ² = {c2:.2f}  →  {np.sqrt(c2):.1f}σ combined")
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

    # ε=0 baseline for comparison
    base = [r for r in results if r['eps'] == 0.0]
    if base:
        b0 = min(base, key=total_chi2)
        print(f"\n  BASELINE (ε=0):  w₀ = {b0['w0']:.4f},  w_a = {b0['wa']:.4f},  Σχ² = {total_chi2(b0):.1f}")

    print(f"\n{'='*78}")
    print(f"  Test 38 COMPLETE")
    print(f"{'='*78}")
