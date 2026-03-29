"""
Test 31 — 2-loop b₀ correction to dimensional transmutation
=============================================================

QUESTION: The transmutation formula Λ_d = μ·exp(−2π/(b₀·α_d(μ))) has
exponential sensitivity to b₀. How much does the 2-loop correction shift Λ_d?

1-loop:  b₀ = (11/3)C₂(G) − (2/3)T(R)N_f = 22/3 − N_f/3

2-loop:  β(α) = −b₀α²/(2π) − b₁α³/(4π²)

where b₁ = (34/3)C₂(G)² − [(20/3)C₂(G) + 4C₂(R)]T(R)N_f

For SU(2): C₂(G) = 2, C₂(R=fund) = 3/4, T(R=fund) = 1/2

PLAN:
    Part 1: Compute b₁ for SU(2)_d with N_f = 3 Majorana
    Part 2: Solve 2-loop RG equation numerically
    Part 3: Compare Λ_d at 1-loop vs 2-loop
    Part 4: Uncertainty band on α_d needed for Λ_d = 2 meV
"""

import numpy as np
import math
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from rg_transmutation import (
    beta_coeff_SU2, transmutation_scale, alpha_needed_for_Lambda,
    M_PL, M_CHI, M_PHI, ALPHA_D, LAMBDA_D_TARGET
)

# ═══════════════════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════════════════
N_F = 3  # Majorana dark quarks


def print_header(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


# ═══════════════════════════════════════════════════════════════════════════
#  Part 1: 2-loop coefficient b₁
# ═══════════════════════════════════════════════════════════════════════════
def part1_two_loop_coefficient():
    print_header("Part 1: 2-loop β-function coefficient b₁")

    # SU(N) group theory:
    C2_G = 2.0       # SU(2): N
    C2_R = 3.0/4.0   # SU(2) fundamental: (N²-1)/(2N) = 3/4
    T_R  = 0.5       # fundamental: 1/2

    b0 = beta_coeff_SU2(N_F)

    # 2-loop coefficient for Dirac fermions:
    # b₁ = (34/3)C₂(G)² − [(20/3)C₂(G) + 4C₂(R)]T(R)·n_D
    #
    # For Majorana fermions: each counts as 1/2 Dirac → n_D = N_f/2
    # OR equivalently: replace T(R)·n_D → T(R)·N_f/2
    # (Majorana contributes half the fermion loop of Dirac)
    #
    # More precisely for the 2-loop coefficient with Weyl/Majorana:
    # b₁ = (34/3)C₂(G)² − [(10/3)C₂(G) + 2C₂(R)]T(R)·N_f_Majorana

    n_D_equiv = N_F / 2.0  # Majorana → half-Dirac

    b1_dirac_formula = (34.0/3.0) * C2_G**2 - \
                       ((20.0/3.0) * C2_G + 4.0 * C2_R) * T_R * n_D_equiv

    print(f"  SU(2)_d group theory:")
    print(f"    C₂(G) = {C2_G}")
    print(f"    C₂(R=fund) = {C2_R}")
    print(f"    T(R=fund) = {T_R}")
    print(f"    N_f = {N_F} Majorana = {n_D_equiv} Dirac equivalent")
    print()
    print(f"  1-loop: b₀ = {b0:.4f}  ({int(22-N_F)}/3)")
    print()

    # Detailed calculation
    term1 = (34.0/3.0) * C2_G**2
    term2 = (20.0/3.0) * C2_G * T_R * n_D_equiv
    term3 = 4.0 * C2_R * T_R * n_D_equiv
    b1 = term1 - term2 - term3

    print(f"  2-loop: b₁ = (34/3)C₂(G)² − [(20/3)C₂(G) + 4C₂(R)]T(R)·n_D")
    print(f"    Term 1: (34/3)·{C2_G}² = {term1:.4f}")
    print(f"    Term 2: (20/3)·{C2_G}·{T_R}·{n_D_equiv} = {term2:.4f}")
    print(f"    Term 3: 4·{C2_R}·{T_R}·{n_D_equiv} = {term3:.4f}")
    print(f"    b₁ = {term1:.4f} − {term2:.4f} − {term3:.4f} = {b1:.4f}")
    print()
    print(f"  Ratio b₁/b₀² = {b1/b0**2:.4f}")
    print(f"  (This controls the relative size of 2-loop correction)")

    return b0, b1


# ═══════════════════════════════════════════════════════════════════════════
#  Part 2: Solve 2-loop RG numerically
# ═══════════════════════════════════════════════════════════════════════════
def part2_two_loop_rg(b0, b1):
    print_header("Part 2: 2-loop RG running")

    # β(g) = −b₀g³/(16π²) − b₁g⁵/(16π²)²
    # In terms of α = g²/(4π):
    # dα/d(ln μ) = β(α) = −b₀α²/(2π) − b₁α³/(4π²)
    #
    # Or equivalently:
    # d(1/α)/d(ln μ) = b₀/(2π) + b₁α/(4π²)

    # Numerical integration from μ = m_χ downward to find Λ_d
    # (where α → ∞, i.e., 1/α → 0)

    # Method: integrate d(1/α)/dt where t = ln(μ/m_χ), from t=0 downward

    from scipy.integrate import solve_ivp

    def rhs_1loop(t, inv_alpha):
        return [b0 / (2 * math.pi)]

    def rhs_2loop(t, inv_alpha):
        if inv_alpha[0] <= 0:
            return [0.0]
        alpha = 1.0 / inv_alpha[0]
        return [b0 / (2 * math.pi) + b1 * alpha / (4 * math.pi**2)]

    # Start at μ = m_χ with α_d needed for Λ_d = 2 meV
    alpha_start = alpha_needed_for_Lambda(LAMBDA_D_TARGET, M_CHI, b0)
    inv_alpha_start = 1.0 / alpha_start

    print(f"  Starting point: μ = m_χ = {M_CHI} GeV")
    print(f"  α_d(m_χ) for Λ_d = 2 meV (1-loop): {alpha_start:.6f}")
    print(f"  1/α_d = {inv_alpha_start:.2f}")
    print()

    # Integrate from t = 0 (μ = m_χ) to large negative t (μ → 0)
    # Λ_d is where 1/α → 0 (confinement)
    # ln(Λ_d/m_χ) = t at which 1/α = 0

    # 1-loop analytic: t* = −2π/(b₀·α) = −inv_alpha·2π/b₀
    t_star_1loop = -inv_alpha_start * 2 * math.pi / b0
    Lambda_1loop = M_CHI * math.exp(t_star_1loop)

    print(f"  1-loop analytic:")
    print(f"    t* = ln(Λ_d/m_χ) = {t_star_1loop:.2f}")
    print(f"    Λ_d = {Lambda_1loop:.4e} GeV = {Lambda_1loop*1e12:.4f} meV")
    print()

    # 2-loop numerical: integrate until 1/α = 0
    t_span = [0, t_star_1loop * 2]  # go further than 1-loop estimate

    # Event: 1/α crosses zero
    def inv_alpha_zero(t, y):
        return y[0]
    inv_alpha_zero.terminal = True
    inv_alpha_zero.direction = -1

    sol_1loop = solve_ivp(rhs_1loop, t_span, [inv_alpha_start],
                          events=inv_alpha_zero, rtol=1e-12, max_step=0.1)
    sol_2loop = solve_ivp(rhs_2loop, t_span, [inv_alpha_start],
                          events=inv_alpha_zero, rtol=1e-12, max_step=0.1)

    if sol_1loop.t_events[0].size > 0:
        t_conf_1loop = sol_1loop.t_events[0][0]
        Lambda_1loop_num = M_CHI * math.exp(t_conf_1loop)
    else:
        t_conf_1loop = t_star_1loop
        Lambda_1loop_num = Lambda_1loop

    if sol_2loop.t_events[0].size > 0:
        t_conf_2loop = sol_2loop.t_events[0][0]
        Lambda_2loop = M_CHI * math.exp(t_conf_2loop)
    else:
        # If no crossing, extrapolate
        t_last = sol_2loop.t[-1]
        inv_last = sol_2loop.y[0, -1]
        # Linear extrapolation
        if len(sol_2loop.t) > 1:
            dt = sol_2loop.t[-1] - sol_2loop.t[-2]
            dinv = sol_2loop.y[0, -1] - sol_2loop.y[0, -2]
            if dinv < 0:
                t_conf_2loop = t_last - inv_last * dt / dinv
                Lambda_2loop = M_CHI * math.exp(t_conf_2loop)
            else:
                Lambda_2loop = 0
                t_conf_2loop = -1e10
        else:
            Lambda_2loop = 0
            t_conf_2loop = -1e10

    print(f"  Numerical integration results:")
    print(f"  {'':>4} {'1-loop':>15} {'2-loop':>15} {'Ratio':>10}")
    print(f"  {'':>4} {'─'*15} {'─'*15} {'─'*10}")
    print(f"  {'t*':>4} {t_conf_1loop:>15.4f} {t_conf_2loop:>15.4f}"
          f" {t_conf_2loop/t_conf_1loop:>10.6f}")
    print(f"  {'Λ_d':>4} {Lambda_1loop_num:>15.4e} {Lambda_2loop:>15.4e}"
          f" {Lambda_2loop/Lambda_1loop_num:>10.4f}")
    print(f"  {'meV':>4} {Lambda_1loop_num*1e12:>15.4f}"
          f" {Lambda_2loop*1e12:>15.4f}"
          f" {Lambda_2loop/Lambda_1loop_num:>10.4f}")

    shift = (Lambda_2loop - Lambda_1loop_num) / Lambda_1loop_num * 100
    print(f"\n  2-loop shift: {shift:+.2f}%")
    print(f"  Λ_d(2-loop) / Λ_d(1-loop) = {Lambda_2loop/Lambda_1loop_num:.4f}")

    return Lambda_1loop_num, Lambda_2loop


# ═══════════════════════════════════════════════════════════════════════════
#  Part 3: Impact on α_d needed for Λ_d = 2 meV
# ═══════════════════════════════════════════════════════════════════════════
def part3_alpha_correction(b0, b1):
    print_header("Part 3: α_d correction for target Λ_d = 2 meV")

    # At 1-loop: α_d(m_χ) = 2π / (b₀ · ln(m_χ/Λ_d))
    alpha_1loop = alpha_needed_for_Lambda(LAMBDA_D_TARGET, M_CHI, b0)

    # At 2-loop: need to solve iteratively
    # The 2-loop Λ parameter is:
    # Λ_2loop = μ · exp(−2π/(b₀α)) · (b₀α/(2π))^(−b₁/(2b₀²))
    #
    # Setting Λ_2loop = Λ_d_target and solving for α:
    # This is transcendental → iterate

    from scipy.optimize import brentq
    from scipy.integrate import solve_ivp

    def rhs_2loop(t, inv_alpha):
        if inv_alpha[0] <= 0:
            return [0.0]
        alpha = 1.0 / inv_alpha[0]
        return [b0 / (2 * math.pi) + b1 * alpha / (4 * math.pi**2)]

    def Lambda_from_alpha_2loop(alpha_start):
        """Find Λ_d from 2-loop RG starting at α(m_χ) = alpha_start."""
        inv_a = 1.0 / alpha_start
        t_span = [0, -inv_a * 2 * math.pi / b0 * 2]

        def inv_alpha_zero(t, y):
            return y[0]
        inv_alpha_zero.terminal = True
        inv_alpha_zero.direction = -1

        sol = solve_ivp(rhs_2loop, t_span, [inv_a],
                        events=inv_alpha_zero, rtol=1e-12, max_step=0.1)
        if sol.t_events[0].size > 0:
            return M_CHI * math.exp(sol.t_events[0][0])
        return 0.0

    # Find α that gives Λ_d = 2 meV at 2-loop
    def residual(alpha):
        L = Lambda_from_alpha_2loop(alpha)
        if L <= 0:
            return -1
        return math.log(L) - math.log(LAMBDA_D_TARGET)

    # Bracket: try around 1-loop value
    alpha_2loop = brentq(residual, alpha_1loop * 0.5, alpha_1loop * 2.0,
                          xtol=1e-10)

    shift_pct = (alpha_2loop - alpha_1loop) / alpha_1loop * 100

    print(f"  Target: Λ_d = {LAMBDA_D_TARGET*1e12:.1f} meV")
    print(f"  μ = m_χ = {M_CHI:.2f} GeV")
    print()
    print(f"  {'':>25} {'1-loop':>15} {'2-loop':>15} {'Shift':>10}")
    print(f"  {'':>25} {'─'*15} {'─'*15} {'─'*10}")
    print(f"  {'α_d(m_χ)':>25} {alpha_1loop:>15.6f} {alpha_2loop:>15.6f}"
          f" {shift_pct:>+9.2f}%")
    print(f"  {'1/α_d':>25} {1/alpha_1loop:>15.2f} {1/alpha_2loop:>15.2f}")
    print()
    print(f"  The 2-loop correction shifts α_d by {shift_pct:+.2f}%")
    print(f"  This is {'negligible' if abs(shift_pct) < 5 else 'significant'}"
          f" given the exponential sensitivity.")

    return alpha_1loop, alpha_2loop


# ═══════════════════════════════════════════════════════════════════════════
#  Part 4: Exponential sensitivity analysis
# ═══════════════════════════════════════════════════════════════════════════
def part4_sensitivity(b0, b1):
    print_header("Part 4: Exponential sensitivity — δΛ_d/Λ_d vs δα_d/α_d")

    # From Λ_d = μ·exp(−2π/(b₀·α)):
    # δ(ln Λ_d) = (2π/(b₀·α²)) · δα = (1/(b₀·α)) · (2π/α) · δα
    # δΛ_d/Λ_d = (2π/(b₀·α²)) · δα

    alpha_needed = alpha_needed_for_Lambda(LAMBDA_D_TARGET, M_CHI, b0)

    amplification = 2 * math.pi / (b0 * alpha_needed**2)

    print(f"  α_d(m_χ) = {alpha_needed:.6f}")
    print(f"  Amplification factor: ∂ln(Λ_d)/∂α = 2π/(b₀α²) = {amplification:.1f}")
    print(f"  → A 1% change in α_d shifts Λ_d by {amplification:.0f}%")
    print()

    # Table: δα_d → δΛ_d
    print(f"  {'δα_d/α_d':>12} {'δΛ_d/Λ_d':>15} {'Λ_d [meV]':>12}")
    print(f"  {'─'*12} {'─'*15} {'─'*12}")

    for delta_pct in [-5, -2, -1, -0.5, 0, 0.5, 1, 2, 5]:
        alpha_new = alpha_needed * (1 + delta_pct/100)
        Lambda_new = transmutation_scale(alpha_new, M_CHI, b0)
        delta_Lambda_pct = (Lambda_new - LAMBDA_D_TARGET) / LAMBDA_D_TARGET * 100
        print(f"  {delta_pct:>+11.1f}% {delta_Lambda_pct:>+14.1f}%"
              f" {Lambda_new*1e12:>12.4f}")

    print()
    print(f"  ┌──────────────────────────────────────────────────────────┐")
    print(f"  │  CONCLUSION:                                            │")
    print(f"  │  2-loop corrections modify α_d by a few %              │")
    print(f"  │  Due to exponential sensitivity, Λ_d shifts by O(1)    │")
    print(f"  │  This is a theoretical uncertainty, not a crisis:       │")
    print(f"  │  α_d is a FREE PARAMETER — 2-loop just shifts the      │")
    print(f"  │  numerical value needed, not the physics.              │")
    print(f"  │                                                         │")
    print(f"  │  Key: Λ_d ∝ exp(−2π/(b₀α))                           │")
    print(f"  │  Amplification factor = {amplification:.0f} for α = {alpha_needed:.4f}        │")
    print(f"  └──────────────────────────────────────────────────────────┘")


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("╔" + "═"*68 + "╗")
    print("║  Test 31: 2-loop β-function correction to Λ_d                    ║")
    print("╚" + "═"*68 + "╝")

    b0, b1 = part1_two_loop_coefficient()
    Lambda_1, Lambda_2 = part2_two_loop_rg(b0, b1)
    alpha_1, alpha_2 = part3_alpha_correction(b0, b1)
    part4_sensitivity(b0, b1)

    print("\n  Done.  Test 31 complete.")
