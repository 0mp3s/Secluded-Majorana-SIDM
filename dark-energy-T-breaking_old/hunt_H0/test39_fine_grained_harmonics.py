"""
Test 39 — Fine-Grained Harmonic Scan
=====================================

Test 38b identified the breakthrough region:
  ε ∈ [-0.12, -0.04],  θ_i ∈ [2.94, 3.03]

This scan zooms in with:
  ε  : -0.03 to -0.13 in steps of 0.01
  θ_i: 2.940 to 3.030 in steps of 0.002

Goal: find the global minimum of the combined χ² = χ²_PP + χ²_U3 + χ²_DY5
and map the 2σ-allowed region precisely.
"""

import numpy as np
import sys, os, time

sys.path.insert(0, os.path.dirname(__file__))

import layer8_cosmic_ode as L8
import desi_comparison as DC

from layer8_cosmic_ode import (
    M_PL, find_Lambda_d_for_H0, H0_PLANCK_KMS,
)
from desi_comparison import extract_w_of_a, fit_cpl

# ── DESI DR2 targets ────────────────────────────────────────────────────
DESI_TARGETS = {
    'PP':   {'w0': -0.838, 'w0_err': 0.055, 'wa': -0.62,  'wa_err': 0.205},
    'U3':   {'w0': -0.667, 'w0_err': 0.088, 'wa': -1.09,  'wa_err': 0.290},
    'DY5':  {'w0': -0.752, 'w0_err': 0.057, 'wa': -0.86,  'wa_err': 0.215},
}

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

    print("=" * 100)
    print("  Test 39 — Fine-Grained Harmonic Scan")
    print("  ε: -0.03 → -0.13 (step 0.01),  θ_i: 2.940 → 3.030 (step 0.002)")
    print("=" * 100)

    eps_values   = np.round(np.arange(-0.03, -0.131, -0.01), 3)
    theta_values = np.round(np.arange(2.940, 3.031, 0.002), 4)

    results = []
    t0 = time.time()

    print(f"\n{'='*105}")
    print(f"  {'ε':>6}  {'θ_i':>6}  {'Λ_d':>11}  {'H₀':>6}  "
          f"{'w₀_CPL':>8}  {'w_a_CPL':>8}  "
          f"{'χ²_PP':>7}  {'χ²_U3':>7}  {'χ²_DY5':>7}  {'Σχ²':>7}")
    print("  " + "-" * 100)

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
            c2_tot = c2_pp + c2_u3 + c2_dy5

            marker = ""
            if c2_tot < 15:
                marker = " ★"
            if c2_tot < 12:
                marker = " ★★"
            if c2_tot < 10:
                marker = " ★★★"

            print(f"  {eps:6.2f}  {theta_i:6.3f}  {Ld:11.4e}  {res.H0_kms:6.2f}  "
                  f"{w0_cpl:8.4f}  {wa_cpl:8.4f}  "
                  f"{c2_pp:7.2f}  {c2_u3:7.2f}  {c2_dy5:7.2f}  {c2_tot:7.2f}{marker}")

            results.append({
                'eps': eps, 'theta_i': theta_i, 'Lambda_d': Ld,
                'H0': res.H0_kms, 'w0': w0_cpl, 'wa': wa_cpl,
                'w_today': res.w_sigma,
                'chi2_PP': c2_pp, 'chi2_U3': c2_u3, 'chi2_DY5': c2_dy5,
                'chi2_tot': c2_tot,
            })

    elapsed = time.time() - t0

    # Restore
    L8.V_sigma  = V_orig
    L8.dV_sigma = dV_orig
    DC.V_sigma  = V_dc_orig

    print(f"\n{'='*105}")
    print(f"  SCAN COMPLETE ({elapsed:.0f}s, {len(results)} valid points)")
    print(f"{'='*105}")

    if not results:
        print("  No valid points!")
        sys.exit(1)

    # ── Best per dataset ────────────────────────────────────────────────
    for ds_name, ds in DESI_TARGETS.items():
        key = f'chi2_{ds_name}'
        best = min(results, key=lambda r: r[key])
        c2 = best[key]
        sigma = np.sqrt(c2)
        print(f"\n  Best for DESI+CMB+{ds_name}:")
        print(f"    ε = {best['eps']:.3f},  θ_i = {best['theta_i']:.4f}")
        print(f"    w₀ = {best['w0']:.4f},  w_a = {best['wa']:.4f}")
        print(f"    χ² = {c2:.2f}  →  {sigma:.1f}σ combined")

    # ── Global best (minimum total χ²) ──────────────────────────────────
    best_tot = min(results, key=lambda r: r['chi2_tot'])
    print(f"\n  {'='*60}")
    print(f"  GLOBAL BEST (min Σχ²):")
    print(f"    ε = {best_tot['eps']:.3f},  θ_i = {best_tot['theta_i']:.4f}")
    print(f"    w₀ = {best_tot['w0']:.4f},  w_a = {best_tot['wa']:.4f}")
    print(f"    χ²_PP = {best_tot['chi2_PP']:.2f},  χ²_U3 = {best_tot['chi2_U3']:.2f},  χ²_DY5 = {best_tot['chi2_DY5']:.2f}")
    print(f"    Σχ² = {best_tot['chi2_tot']:.2f}")
    print(f"    Λ_d = {best_tot['Lambda_d']:.4e} GeV,  H₀ = {best_tot['H0']:.2f}")
    print(f"  {'='*60}")

    # ── Top 10 points ───────────────────────────────────────────────────
    top10 = sorted(results, key=lambda r: r['chi2_tot'])[:10]
    print(f"\n  TOP 10 POINTS (by Σχ²):")
    print(f"  {'#':>3}  {'ε':>6}  {'θ_i':>6}  {'w₀':>8}  {'w_a':>8}  {'χ²_PP':>7}  {'χ²_U3':>7}  {'χ²_DY5':>7}  {'Σχ²':>7}")
    print("  " + "-" * 75)
    for i, r in enumerate(top10, 1):
        print(f"  {i:3d}  {r['eps']:6.3f}  {r['theta_i']:6.3f}  {r['w0']:8.4f}  {r['wa']:8.4f}  "
              f"{r['chi2_PP']:7.2f}  {r['chi2_U3']:7.2f}  {r['chi2_DY5']:7.2f}  {r['chi2_tot']:7.2f}")

    # ── Save to CSV ─────────────────────────────────────────────────────
    csv_path = os.path.join(os.path.dirname(__file__), 'test39_fine_results.csv')
    with open(csv_path, 'w') as f:
        f.write('eps,theta_i,Lambda_d,H0,w0,wa,chi2_PP,chi2_U3,chi2_DY5,chi2_tot\n')
        for r in sorted(results, key=lambda x: x['chi2_tot']):
            f.write(f"{r['eps']:.3f},{r['theta_i']:.4f},{r['Lambda_d']:.6e},"
                    f"{r['H0']:.2f},{r['w0']:.5f},{r['wa']:.5f},"
                    f"{r['chi2_PP']:.3f},{r['chi2_U3']:.3f},{r['chi2_DY5']:.3f},{r['chi2_tot']:.3f}\n")
    print(f"\n  Results saved to {csv_path}")
