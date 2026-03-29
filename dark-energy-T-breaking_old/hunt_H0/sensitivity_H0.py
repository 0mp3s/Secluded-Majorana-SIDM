"""
Test 33 — Sensitivity analysis: ∂H₀/∂(parameters)
====================================================

QUESTION: How sensitive is the H₀ prediction to each free parameter?
What error bars should we quote on H₀?

Our model has 2 truly free dark-QCD parameters: {α_d, θ_i}
(Λ_d is derived from α_d via transmutation; f is fixed by the ODE fit.)

But observationally, the direct inputs to the ODE are {f, Λ_d, θ_i},
so we compute sensitivities to all three.

PLAN:
    Part 1: ∂H₀/∂θ_i  (central finite differences)
    Part 2: ∂H₀/∂Λ_d  (central finite differences)
    Part 3: ∂H₀/∂f    (central finite differences)
    Part 4: Error propagation → σ(H₀) from parameter uncertainties
    Part 5: Jacobian summary — which parameter matters most?
"""

import numpy as np
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, '..', '..', 'core'))

from layer8_cosmic_ode import (
    solve_layer8, find_Lambda_d_for_H0, compute_omega_chi_h2,
    M_PL, H_100_GEV, H0_PLANCK_KMS, H0_SHOES_KMS
)

# ═══════════════════════════════════════════════════════════════════════════
#  Parameters (from best-fit benchmarks)
# ═══════════════════════════════════════════════════════════════════════════
M_CHI      = 98.19
M_PHI      = 9.66e-3
ALPHA_D    = 3.274e-3
LAMBDA_D   = 2.0e-12       # 2 meV in GeV
F_DEFAULT  = 0.27 * M_PL   # from PI-7 / Test 28
THETA_I    = 3.0            # near-hilltop

# Benchmark points
BENCHMARKS = {
    'Hilltop (θ=3.0)': {
        'theta_i': 3.0, 'Lambda_d': 2.0e-12, 'f': 0.27 * M_PL,
    },
    'θ=3.1 (SH0ES)': {
        'theta_i': 3.1, 'Lambda_d': 2.0e-12, 'f': 0.27 * M_PL,
    },
}


def print_header(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


def compute_H0(theta_i, Lambda_d, f, omega_chi_h2):
    """Wrapper around solve_layer8 returning H₀ in km/s/Mpc."""
    r = solve_layer8(M_CHI, M_PHI, ALPHA_D, f, Lambda_d, theta_i,
                     omega_chi_h2=omega_chi_h2, verbose=False)
    return r.H0_kms if r.H0_kms is not None else np.nan


def numerical_derivative(func, x0, dx, omega_chi_h2, other_params):
    """Central finite difference: df/dx at x0 with step dx."""
    f_plus  = func(x0 + dx, omega_chi_h2=omega_chi_h2, **other_params)
    f_minus = func(x0 - dx, omega_chi_h2=omega_chi_h2, **other_params)
    if np.isnan(f_plus) or np.isnan(f_minus):
        return np.nan
    return (f_plus - f_minus) / (2.0 * dx)


# ═══════════════════════════════════════════════════════════════════════════
#  Part 1: ∂H₀/∂θ_i
# ═══════════════════════════════════════════════════════════════════════════
def part1_dH0_dtheta(omega_chi_h2):
    print_header("Part 1: ∂H₀/∂θ_i — sensitivity to misalignment angle")

    for name, params in BENCHMARKS.items():
        theta_i = params['theta_i']
        Ld      = params['Lambda_d']
        f       = params['f']

        print(f"  ── {name} ──")
        print(f"  θ_i = {theta_i},  f/M_Pl = {f/M_PL:.2f},  Λ_d = {Ld*1e12:.1f} meV")

        H0_center = compute_H0(theta_i, Ld, f, omega_chi_h2)
        print(f"  H₀(center) = {H0_center:.3f} km/s/Mpc")

        # Multiple step sizes for convergence check
        print(f"\n  {'δθ':>10}  {'H₀(θ+δ)':>12}  {'H₀(θ-δ)':>12}  "
              f"{'∂H₀/∂θ':>12}  {'∂ln H₀/∂θ':>12}")
        print(f"  {'─'*10}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*12}")

        best_deriv = None
        for dth in [0.1, 0.05, 0.01, 0.005, 0.001]:
            H_plus  = compute_H0(theta_i + dth, Ld, f, omega_chi_h2)
            H_minus = compute_H0(theta_i - dth, Ld, f, omega_chi_h2)
            dHdth = (H_plus - H_minus) / (2.0 * dth)
            dln = dHdth / H0_center if H0_center > 0 else np.nan
            print(f"  {dth:>10.4f}  {H_plus:>12.3f}  {H_minus:>12.3f}  "
                  f"{dHdth:>12.3f}  {dln:>12.4f}")
            if dth == 0.01:
                best_deriv = dHdth

        if best_deriv is not None:
            # What θ_i shift gives ΔH₀ = 1 km/s/Mpc?
            dtheta_for_1kms = 1.0 / abs(best_deriv) if best_deriv != 0 else np.inf
            print(f"\n  ∂H₀/∂θ_i = {best_deriv:.3f} km/s/Mpc per radian  (at δθ=0.01)")
            print(f"  To shift H₀ by 1 km/s/Mpc: Δθ_i = {dtheta_for_1kms:.4f} rad")
            print(f"  To bridge Hubble tension (5.6 km/s/Mpc): "
                  f"Δθ_i = {5.64/abs(best_deriv):.4f} rad")
        print()


# ═══════════════════════════════════════════════════════════════════════════
#  Part 2: ∂H₀/∂Λ_d
# ═══════════════════════════════════════════════════════════════════════════
def part2_dH0_dLambda(omega_chi_h2):
    print_header("Part 2: ∂H₀/∂Λ_d — sensitivity to confinement scale")

    for name, params in BENCHMARKS.items():
        theta_i = params['theta_i']
        Ld      = params['Lambda_d']
        f       = params['f']

        print(f"  ── {name} ──")
        print(f"  θ_i = {theta_i},  f/M_Pl = {f/M_PL:.2f},  Λ_d = {Ld*1e12:.1f} meV")

        H0_center = compute_H0(theta_i, Ld, f, omega_chi_h2)
        print(f"  H₀(center) = {H0_center:.3f} km/s/Mpc")

        print(f"\n  {'δΛ_d/Λ_d':>10}  {'H₀(Λ+δ)':>12}  {'H₀(Λ-δ)':>12}  "
              f"{'∂H₀/∂ln Λ_d':>14}  {'∂ln H₀/∂ln Λ_d':>16}")
        print(f"  {'─'*10}  {'─'*12}  {'─'*12}  {'─'*14}  {'─'*16}")

        best_deriv = None
        for frac in [0.1, 0.05, 0.01, 0.005, 0.001]:
            dLd = frac * Ld
            H_plus  = compute_H0(theta_i, Ld + dLd, f, omega_chi_h2)
            H_minus = compute_H0(theta_i, Ld - dLd, f, omega_chi_h2)
            dH_dLd = (H_plus - H_minus) / (2.0 * dLd)
            # Convert to ∂H₀/∂(ln Λ_d) = Λ_d · ∂H₀/∂Λ_d
            dH_dlnLd = Ld * dH_dLd
            dln = dH_dlnLd / H0_center if H0_center > 0 else np.nan
            print(f"  {frac:>10.4f}  {H_plus:>12.3f}  {H_minus:>12.3f}  "
                  f"{dH_dlnLd:>14.3f}  {dln:>16.4f}")
            if frac == 0.01:
                best_deriv = dH_dlnLd

        if best_deriv is not None:
            pct_for_1kms = 1.0 / abs(best_deriv) * 100 if best_deriv != 0 else np.inf
            print(f"\n  ∂H₀/∂(ln Λ_d) = {best_deriv:.3f} km/s/Mpc per e-fold")
            print(f"  To shift H₀ by 1 km/s/Mpc: δΛ_d/Λ_d = {pct_for_1kms:.2f}%")
        print()


# ═══════════════════════════════════════════════════════════════════════════
#  Part 3: ∂H₀/∂f
# ═══════════════════════════════════════════════════════════════════════════
def part3_dH0_df(omega_chi_h2):
    print_header("Part 3: ∂H₀/∂f — sensitivity to decay constant")

    for name, params in BENCHMARKS.items():
        theta_i = params['theta_i']
        Ld      = params['Lambda_d']
        f       = params['f']

        print(f"  ── {name} ──")
        print(f"  θ_i = {theta_i},  f/M_Pl = {f/M_PL:.2f},  Λ_d = {Ld*1e12:.1f} meV")

        H0_center = compute_H0(theta_i, Ld, f, omega_chi_h2)
        print(f"  H₀(center) = {H0_center:.3f} km/s/Mpc")

        print(f"\n  {'δf/f':>10}  {'H₀(f+δ)':>12}  {'H₀(f-δ)':>12}  "
              f"{'∂H₀/∂ln f':>14}  {'∂ln H₀/∂ln f':>16}")
        print(f"  {'─'*10}  {'─'*12}  {'─'*12}  {'─'*14}  {'─'*16}")

        best_deriv = None
        for frac in [0.1, 0.05, 0.01, 0.005, 0.001]:
            df = frac * f
            H_plus  = compute_H0(theta_i, Ld, f + df, omega_chi_h2)
            H_minus = compute_H0(theta_i, Ld, f - df, omega_chi_h2)
            dH_df   = (H_plus - H_minus) / (2.0 * df)
            dH_dlnf = f * dH_df
            dln = dH_dlnf / H0_center if H0_center > 0 else np.nan
            print(f"  {frac:>10.4f}  {H_plus:>12.3f}  {H_minus:>12.3f}  "
                  f"{dH_dlnf:>14.3f}  {dln:>16.4f}")
            if frac == 0.01:
                best_deriv = dH_dlnf

        if best_deriv is not None:
            pct_for_1kms = 1.0 / abs(best_deriv) * 100 if best_deriv != 0 else np.inf
            print(f"\n  ∂H₀/∂(ln f) = {best_deriv:.3f} km/s/Mpc per e-fold")
            print(f"  To shift H₀ by 1 km/s/Mpc: δf/f = {pct_for_1kms:.2f}%")
        print()


# ═══════════════════════════════════════════════════════════════════════════
#  Part 4: Error propagation
# ═══════════════════════════════════════════════════════════════════════════
def part4_error_propagation(omega_chi_h2):
    print_header("Part 4: Error propagation → σ(H₀)")

    print("  We ask: what is the theoretical uncertainty on H₀?")
    print("  This requires estimating uncertainties on input parameters.\n")

    # Use the Hilltop benchmark
    theta_i = 3.0
    Ld      = 2.0e-12
    f       = 0.27 * M_PL

    H0_center = compute_H0(theta_i, Ld, f, omega_chi_h2)
    print(f"  Central value: H₀ = {H0_center:.3f} km/s/Mpc\n")

    # Compute partial derivatives (δ = 1%)
    dth = 0.01
    H_p = compute_H0(theta_i + dth, Ld, f, omega_chi_h2)
    H_m = compute_H0(theta_i - dth, Ld, f, omega_chi_h2)
    dH_dtheta = (H_p - H_m) / (2.0 * dth)

    dLd = 0.01 * Ld
    H_p = compute_H0(theta_i, Ld + dLd, f, omega_chi_h2)
    H_m = compute_H0(theta_i, Ld - dLd, f, omega_chi_h2)
    dH_dlnLd = Ld * (H_p - H_m) / (2.0 * dLd)

    df = 0.01 * f
    H_p = compute_H0(theta_i, Ld, f + df, omega_chi_h2)
    H_m = compute_H0(theta_i, Ld, f - df, omega_chi_h2)
    dH_dlnf = f * (H_p - H_m) / (2.0 * df)

    print(f"  Jacobian at (θ_i=3.0, Λ_d=2 meV, f=0.27 M_Pl):")
    print(f"    ∂H₀/∂θ_i      = {dH_dtheta:+.3f}  km/s/Mpc / rad")
    print(f"    ∂H₀/∂(ln Λ_d) = {dH_dlnLd:+.3f}  km/s/Mpc / e-fold")
    print(f"    ∂H₀/∂(ln f)   = {dH_dlnf:+.3f}  km/s/Mpc / e-fold")

    print(f"\n  Assumed parameter uncertainties (theory-motivated):")

    # θ_i: varies over [0, π], our prediction is θ_i ~ 3.0
    # Uncertainty: θ_i is a cosmological initial condition set by inflation
    # Reasonable range: ±0.1 rad (3% of π)
    sigma_theta = 0.1  # rad
    # Λ_d: set by transmutation, α_d has ~1% 2-loop uncertainty
    # From Test 31: amplification factor = 1002 → 1% α_d → O(1) Λ_d
    # We take δΛ_d/Λ_d ~ 50% (generous, covers 2-loop ambiguity)
    sigma_lnLd = 0.5   # e-folds (50%)
    # f: determined by V(σ)=ρ_Λ constraint, uncertainty ~10%
    sigma_lnf  = 0.1   # e-folds (10%)

    print(f"    σ(θ_i)       = {sigma_theta:.2f} rad")
    print(f"    σ(ln Λ_d)    = {sigma_lnLd:.2f}  (δΛ_d/Λ_d ~ {sigma_lnLd*100:.0f}%)")
    print(f"    σ(ln f)      = {sigma_lnf:.2f}  (δf/f ~ {sigma_lnf*100:.0f}%)")

    # Error propagation: σ²(H₀) = Σ (∂H₀/∂x_i)² σ²(x_i)
    sigma_H0_theta = abs(dH_dtheta) * sigma_theta
    sigma_H0_Ld    = abs(dH_dlnLd) * sigma_lnLd
    sigma_H0_f     = abs(dH_dlnf) * sigma_lnf
    sigma_H0_total = np.sqrt(sigma_H0_theta**2 + sigma_H0_Ld**2 + sigma_H0_f**2)

    print(f"\n  Error contributions:")
    print(f"    From θ_i:    σ(H₀)|_θ  = {sigma_H0_theta:.2f} km/s/Mpc")
    print(f"    From Λ_d:    σ(H₀)|_Λ  = {sigma_H0_Ld:.2f} km/s/Mpc")
    print(f"    From f:      σ(H₀)|_f  = {sigma_H0_f:.2f} km/s/Mpc")
    print(f"    ──────────────────────────────────────────")
    print(f"    Total:       σ(H₀)     = {sigma_H0_total:.2f} km/s/Mpc")
    print(f"\n  → H₀ = {H0_center:.1f} ± {sigma_H0_total:.1f} km/s/Mpc")

    # Relative errors
    print(f"\n  Fractional contributions to σ²(H₀):")
    total_var = sigma_H0_total**2
    if total_var > 0:
        print(f"    θ_i:  {sigma_H0_theta**2/total_var*100:.1f}%")
        print(f"    Λ_d:  {sigma_H0_Ld**2/total_var*100:.1f}%")
        print(f"    f:    {sigma_H0_f**2/total_var*100:.1f}%")

    return dH_dtheta, dH_dlnLd, dH_dlnf


# ═══════════════════════════════════════════════════════════════════════════
#  Part 5: Jacobian summary
# ═══════════════════════════════════════════════════════════════════════════
def part5_jacobian_summary(omega_chi_h2, dH_dtheta, dH_dlnLd, dH_dlnf):
    print_header("Part 5: Jacobian summary — which parameter controls H₀?")

    print("  ┌────────────────────────────────────────────────────────────┐")
    print("  │  Parameter   │ ∂H₀/∂x           │ Nature           │ Role │")
    print("  ├──────────────┼───────────────────┼──────────────────┼──────┤")
    print(f"  │  θ_i         │ {dH_dtheta:>+8.3f} /rad     │ Initial condition│ Dial │")
    print(f"  │  ln Λ_d      │ {dH_dlnLd:>+8.3f} /e-fold  │ Derived (α_d)    │ Amp  │")
    print(f"  │  ln f        │ {dH_dlnf:>+8.3f} /e-fold  │ Fixed by V=ρ_Λ   │ Mass │")
    print("  └────────────────────────────────────────────────────────────┘")

    print(f"\n  Interpretation:")
    print(f"  • θ_i is the 'dial': small changes in θ_i tune H₀ precisely")
    print(f"  • Λ_d sets the amplitude: V ∝ Λ_d⁴")
    print(f"  • f sets the dynamics: m_σ = Λ_d²/f controls whether σ oscillates")

    # Compare: Planck (θ=3.0) vs SH0ES (θ=3.1)
    print(f"\n  ── Hubble tension in θ_i units ──")
    H0_30 = compute_H0(3.0, 2.0e-12, 0.27 * M_PL, omega_chi_h2)
    H0_31 = compute_H0(3.1, 2.0e-12, 0.27 * M_PL, omega_chi_h2)
    print(f"  H₀(θ_i=3.0) = {H0_30:.2f} km/s/Mpc")
    print(f"  H₀(θ_i=3.1) = {H0_31:.2f} km/s/Mpc")
    print(f"  δθ = 0.1 rad → δH₀ = {H0_31 - H0_30:.2f} km/s/Mpc")
    if abs(dH_dtheta) > 0:
        dtheta_tension = (H0_SHOES_KMS - H0_PLANCK_KMS) / dH_dtheta
        print(f"  To go from Planck→SH0ES: Δθ_i = {dtheta_tension:.3f} rad")
        print(f"  (i.e. θ_i goes from {3.0:.2f} to {3.0+dtheta_tension:.2f})")

    print(f"\n  ┌──────────────────────────────────────────────────────────┐")
    print(f"  │  CONCLUSION:                                            │")
    print(f"  │  H₀ is a calculable output of θ_i (set by inflation)  │")
    print(f"  │  and Λ_d (set by transmutation).                       │")
    print(f"  │                                                         │")
    print(f"  │  The model has 2 genuinely free parameters: α_d, θ_i  │")
    print(f"  │  α_d → Λ_d (via transmutation, known b₀)             │")
    print(f"  │  f → fixed by condition V(σ_today) = ρ_Λ             │")
    print(f"  │  θ_i → directly maps to H₀ once Λ_d, f are set       │")
    print(f"  │                                                         │")
    print(f"  │  This is NOT fine-tuning: θ_i is a continuous initial  │")
    print(f"  │  condition, and Λ_d is completely determined by α_d.   │")
    print(f"  └──────────────────────────────────────────────────────────┘")


# ═══════════════════════════════════════════════════════════════════════════
#  main
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("╔════════════════════════════════════════════════════════════════════╗")
    print("║  Test 33: Sensitivity analysis — ∂H₀/∂(parameters)             ║")
    print("╚════════════════════════════════════════════════════════════════════╝")

    # Compute relic density once
    omega_chi = compute_omega_chi_h2(M_CHI, ALPHA_D)
    print(f"\n  Ω_χ h² = {omega_chi:.4f}")

    part1_dH0_dtheta(omega_chi)
    part2_dH0_dLambda(omega_chi)
    part3_dH0_df(omega_chi)
    dH_dtheta, dH_dlnLd, dH_dlnf = part4_error_propagation(omega_chi)
    part5_jacobian_summary(omega_chi, dH_dtheta, dH_dlnLd, dH_dlnf)

    print("\n  Done.  Test 33 complete.")
