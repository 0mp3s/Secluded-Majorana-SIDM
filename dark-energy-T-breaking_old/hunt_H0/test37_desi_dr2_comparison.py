"""
Test 37 — DESI DR2 Comparison: Scan θ_i to find best match for w₀, w_a
=========================================================================

DESI DR2 (arXiv:2503.14738, Phys. Rev. D 112, 083515) best-fit values:

  DESI+CMB+Pantheon+:  w₀ = -0.838 ± 0.055,  w_a = -0.62 (+0.22/-0.19)
  DESI+CMB+Union3:     w₀ = -0.667 ± 0.088,  w_a = -1.09 (+0.31/-0.27)
  DESI+CMB+DESY5:      w₀ = -0.752 ± 0.057,  w_a = -0.86 (+0.23/-0.20)
  DESI+CMB (no SNe):   w₀ = -0.42  ± 0.21,   w_a = -1.75 ± 0.58

Strategy:
---------
For each θ_i, use find_Lambda_d_for_H0 (analytic guess + broad scan +
binary refinement) to get Λ_d → H₀ = 67.4, then extract w(a) and fit CPL.
Reuses existing infrastructure from layer8_cosmic_ode and desi_comparison.
"""

import numpy as np
import sys, os, time

sys.path.insert(0, os.path.dirname(__file__))
from layer8_cosmic_ode import (
    M_PL, H_100_GEV, OMEGA_B_H2, OMEGA_R_H2,
    find_Lambda_d_for_H0, H0_PLANCK_KMS,
)
from desi_comparison import extract_w_of_a, fit_cpl

# ──────────────────────────────────────────────────────────────────────────
#  DESI DR2 targets (Section VII.1, arXiv:2503.14738, Eq 26-28)
# ──────────────────────────────────────────────────────────────────────────
DESI_TARGETS = {
    'DESI+CMB+PP':   {'w0': -0.838, 'w0_err': 0.055,
                      'wa': -0.62,  'wa_err': 0.205, 'sig': 2.8},
    'DESI+CMB+U3':   {'w0': -0.667, 'w0_err': 0.088,
                      'wa': -1.09,  'wa_err': 0.290, 'sig': 3.8},
    'DESI+CMB+DESY5':{'w0': -0.752, 'w0_err': 0.057,
                      'wa': -0.86,  'wa_err': 0.215, 'sig': 4.2},
}

# ── MAP benchmark parameters (from global_config.json) ──────────────────
M_CHI   = 98.19       # GeV
ALPHA_D = 3.274e-3
F_REF   = 0.27 * M_PL # GeV
OMEGA   = 0.120


def chi2_for_dataset(w0, wa, dataset):
    dw0 = (w0 - dataset['w0']) / dataset['w0_err']
    dwa = (wa - dataset['wa']) / dataset['wa_err']
    return dw0**2 + dwa**2


# ══════════════════════════════════════════════════════════════════════════
#  Main scan
# ══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("="*72)
    print("  Test 37 — DESI DR2 w₀-w_a comparison")
    print("  Scanning θ_i → find Λ_d(H₀=67.4) → CPL fit → vs DESI DR2")
    print("="*72)

    theta_values = np.arange(2.5, 3.15, 0.05)

    results = []
    t0 = time.time()

    print(f"\n  {'θ_i':>6s}  {'Λ_d [GeV]':>12s}  {'H₀':>6s}  "
          f"{'w₀(CPL)':>8s}  {'wₐ(CPL)':>8s}  {'Ω_DE':>5s}  "
          f"{'χ²_PP':>6s}  {'χ²_U3':>6s}  {'χ²_DY5':>6s}")
    print("  " + "-" * 90)

    for theta_i in theta_values:
        try:
            Ld, res = find_Lambda_d_for_H0(
                H0_PLANCK_KMS, M_CHI, ALPHA_D, F_REF, theta_i,
                omega_chi_h2=OMEGA, tol=0.05)
        except Exception:
            print(f"  {theta_i:6.3f}  FAILED (solver error)")
            continue

        if res is None or res.H0_kms is None:
            print(f"  {theta_i:6.3f}  FAILED (no Λ_d found)")
            continue

        if abs(res.H0_kms - H0_PLANCK_KMS) > 1.0:
            print(f"  {theta_i:6.3f}  FAILED (H₀={res.H0_kms:.1f}, too far)")
            continue

        # Extract w(a) and fit CPL
        a_arr, w_arr = extract_w_of_a(res, F_REF, Ld)
        w0_cpl, wa_cpl, _, _ = fit_cpl(a_arr, w_arr, z_max_fit=2.0)

        chi2_vals = {n: chi2_for_dataset(w0_cpl, wa_cpl, d)
                     for n, d in DESI_TARGETS.items()}

        print(f"  {theta_i:6.3f}  {Ld:12.4e}  {res.H0_kms:6.2f}  "
              f"{w0_cpl:8.5f}  {wa_cpl:8.4f}  {res.Omega_DE:5.3f}  "
              f"{chi2_vals['DESI+CMB+PP']:6.2f}  "
              f"{chi2_vals['DESI+CMB+U3']:6.2f}  "
              f"{chi2_vals['DESI+CMB+DESY5']:6.2f}")

        results.append({
            'theta_i': theta_i,
            'Lambda_d': Ld,
            'H0': res.H0_kms,
            'w0_cpl': w0_cpl,
            'wa_cpl': wa_cpl,
            'w_today': res.w_sigma,
            'Omega_DE': res.Omega_DE,
            'chi2': chi2_vals,
        })

    elapsed = time.time() - t0

    # ── Summary ──
    print(f"\n{'='*72}")
    print(f"  SCAN COMPLETE ({elapsed:.1f}s)")
    print(f"{'='*72}\n")

    if not results:
        print("  No valid points found!")
        sys.exit(1)

    # Find best θ_i for each dataset
    for ds_name in DESI_TARGETS:
        best = min(results, key=lambda r: r['chi2'][ds_name])
        chi2_val = best['chi2'][ds_name]
        dsig = DESI_TARGETS[ds_name]['sig']
        print(f"  Best for {ds_name} ({dsig}σ vs ΛCDM):")
        print(f"    θ_i = {best['theta_i']:.3f} rad")
        print(f"    Λ_d = {best['Lambda_d']:.4e} GeV")
        print(f"    w₀  = {best['w0_cpl']:.4f}  (DESI: {DESI_TARGETS[ds_name]['w0']:.3f} ± {DESI_TARGETS[ds_name]['w0_err']:.3f})")
        print(f"    w_a = {best['wa_cpl']:.4f}  (DESI: {DESI_TARGETS[ds_name]['wa']:.2f} ± {DESI_TARGETS[ds_name]['wa_err']:.3f})")
        print(f"    χ²  = {chi2_val:.2f}  →  Δσ = {np.sqrt(chi2_val):.1f}σ combined")
        print(f"    H₀  = {best['H0']:.2f} km/s/Mpc,  Ω_DE = {best['Omega_DE']:.4f}")
        print()

    # Global best (min total χ² across all 3 datasets)
    def total_chi2(r):
        return sum(r['chi2'].values())

    global_best = min(results, key=total_chi2)
    tc = total_chi2(global_best)
    print(f"  GLOBAL BEST FIT (sum of all χ²):")
    print(f"    θ_i    = {global_best['theta_i']:.3f} rad")
    print(f"    Λ_d    = {global_best['Lambda_d']:.4e} GeV")
    print(f"    w₀     = {global_best['w0_cpl']:.4f}")
    print(f"    w_a    = {global_best['wa_cpl']:.4f}")
    print(f"    H₀     = {global_best['H0']:.2f} km/s/Mpc")
    print(f"    Ω_DE   = {global_best['Omega_DE']:.4f}")
    print(f"    Σχ²    = {tc:.2f}")

    # Full scan table
    print(f"\n  Full θ_i scan table:")
    print(f"  {'θ_i':>6}  {'Λ_d':>11}  {'H₀':>6}  {'w₀_CPL':>8}  {'w_a_CPL':>8}  {'w(z=0)':>8}  {'Ω_DE':>6}")
    print(f"  {'─'*6}  {'─'*11}  {'─'*6}  {'─'*8}  {'─'*8}  {'─'*8}  {'─'*6}")
    for r in results:
        print(f"  {r['theta_i']:6.3f}  {r['Lambda_d']:11.4e}  {r['H0']:6.2f}  "
              f"{r['w0_cpl']:8.4f}  {r['wa_cpl']:8.4f}  {r['w_today']:8.6f}  "
              f"{r['Omega_DE']:6.4f}")

    print(f"\n  DESI DR2 targets for reference:")
    for name, ds in DESI_TARGETS.items():
        print(f"    {name}: w₀ = {ds['w0']:.3f} ± {ds['w0_err']:.3f}, "
              f"w_a = {ds['wa']:.2f} ± {ds['wa_err']:.3f}  ({ds['sig']}σ)")

    print(f"\n{'='*72}")
    print(f"  Test 37 COMPLETE")
    print(f"{'='*72}")
