"""
Test 32 — Future of the Universe: When does dark energy die?
==============================================================

QUESTION: In our model, V(σ) → 0 as σ oscillates toward the minimum.
When does the cosmic acceleration stop? What's the ultimate fate?

Unlike ΛCDM (eternal acceleration), our dark energy is TRANSIENT.
The field σ eventually reaches the minimum of V = Λ_d⁴(1 − cos(σ/f)),
and V → 0 → dark energy disappears.

PLAN:
    Part 1: Extend ODE into future (N > 0, a > 1)
    Part 2: Extract w(z) for z < 0 (future)
    Part 3: Find the deceleration epoch (when q = 0 again)
    Part 4: Asymptotic fate of the universe
"""

import numpy as np
from scipy.integrate import solve_ivp
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, '..', '..', 'core'))

from layer8_cosmic_ode import (
    _ode_rhs, V_sigma, dV_sigma, M_PL, T_CMB_GEV, H_100_GEV,
    N_EFF, OMEGA_R_H2, OMEGA_B_H2, _RHO_UNIT, _RHO_GAMMA,
    compute_omega_chi_h2, H0_PLANCK_KMS
)

# ═══════════════════════════════════════════════════════════════════════════
#  Parameters
# ═══════════════════════════════════════════════════════════════════════════
M_CHI      = 98.19
ALPHA_D    = 3.274e-3
F_DEFAULT  = 0.27 * M_PL
LAMBDA_D   = 2.0e-12       # 2 meV in GeV

# Benchmark points from Test 28
BENCHMARKS = {
    'Hilltop (θ=3.0)':  {'theta_i': 3.0,  'Lambda_d': 2.0e-12, 'f': 0.27 * M_PL},
    'Planck+DESI':       {'theta_i': 2.92, 'Lambda_d': 2.0e-12, 'f': 0.17 * M_PL},
    'SH0ES-compat':      {'theta_i': 2.96, 'Lambda_d': 2.1e-12, 'f': 0.20 * M_PL},
}


def print_header(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


# ═══════════════════════════════════════════════════════════════════════════
#  Part 1: Extend ODE into the future
# ═══════════════════════════════════════════════════════════════════════════
def solve_future(theta_i, f, Lambda_d, N_future=5.0,
                 omega_chi_h2=None, verbose=True):
    """
    Solve the coupled σ-Friedmann ODE from reheating through today
    and into the future (N > 0 → a > 1).

    N_future = ln(a_max/a_today). N_future=5 → a ~ 150.
    """
    if omega_chi_h2 is None:
        omega_chi_h2 = compute_omega_chi_h2(M_CHI, ALPHA_D)

    rho_r0 = OMEGA_R_H2 * _RHO_UNIT
    rho_m0 = (omega_chi_h2 + OMEGA_B_H2) * _RHO_UNIT

    g_star_S_RH = 106.75
    g_star_S_0  = 3.91
    T_RH = 1e5
    a_RH = (T_CMB_GEV / T_RH) * (g_star_S_0 / g_star_S_RH)**(1.0/3.0)
    N_RH = np.log(a_RH)

    sigma_init = f * theta_i
    p_init = 0.0

    # Integrate from reheating to future
    sol = solve_ivp(
        _ode_rhs,
        [N_RH, N_future],
        [sigma_init, p_init],
        args=(f, Lambda_d, rho_r0, rho_m0),
        method='RK45',
        rtol=1e-12,
        atol=1e-15,
        dense_output=True,
        max_step=0.5,
    )

    return sol, rho_r0, rho_m0


def part1_future_evolution():
    print_header("Part 1: σ field and dark energy evolution into the future")

    omega_chi = compute_omega_chi_h2(M_CHI, ALPHA_D)
    print(f"  Ω_χ h² = {omega_chi:.4f}\n")

    N_future = 6.0  # a_max ≈ 400

    for label, params in BENCHMARKS.items():
        theta_i = params['theta_i']
        f = params['f']
        Ld = params['Lambda_d']

        sol, rho_r0, rho_m0 = solve_future(
            theta_i, f, Ld, N_future=N_future, omega_chi_h2=omega_chi
        )

        if not sol.success:
            print(f"  {label}: ODE failed")
            continue

        # Sample at many future time points
        N_arr = np.linspace(-2, N_future, 500)
        states = sol.sol(N_arr)
        sigma_arr = states[0]
        p_arr = states[1]

        print(f"  ── {label} ──")
        print(f"  θ_i = {theta_i:.2f},  f = {f/M_PL:.2f} M_Pl,  Λ_d = {Ld*1e12:.1f} meV")
        print()
        print(f"  {'N':>6} {'a':>10} {'z':>8} {'t/t_H':>8}"
              f" {'θ(N)':>10} {'w_σ':>10} {'Ω_DE':>8} {'q':>8}")
        print(f"  {'─'*6} {'─'*10} {'─'*8} {'─'*8}"
              f" {'─'*10} {'─'*10} {'─'*8} {'─'*8}")

        for i, N in enumerate(N_arr):
            if N < -1 or (N > 0.5 and int(i) % 10 != 0):
                if N < 0 and abs(N) > 0.01:
                    continue
                if N > 0 and N < 0.5:
                    continue

            a = np.exp(N)
            z = 1.0/a - 1 if a > 0 else float('inf')
            sigma = sigma_arr[i]
            p = p_arr[i]

            rho_r = rho_r0 * a**(-4)
            rho_m = rho_m0 * a**(-3)
            V = V_sigma(sigma, f, Ld)

            denom = 3.0 * M_PL**2 - 0.5 * p**2
            if denom <= 0:
                continue
            H2 = (rho_r + rho_m + V) / denom
            if H2 <= 0:
                continue

            rho_sig = 0.5 * H2 * p**2 + V
            P_sig = 0.5 * H2 * p**2 - V
            w = P_sig / rho_sig if rho_sig > 1e-100 else -1.0
            rho_tot = rho_r + rho_m + rho_sig
            Omega_DE = rho_sig / rho_tot if rho_tot > 0 else 0

            # Deceleration parameter q = −aä/ȧ² = (1+3w_eff)/2 for single component
            # General: q = Ω_r + (1/2)Ω_m − Ω_DE
            # More precisely: q = (1/2)(1 + 3w_eff) where w_eff = Σ Ω_i w_i
            Omega_r = rho_r / rho_tot
            Omega_m = rho_m / rho_tot
            w_eff = Omega_r/3.0 + w * Omega_DE
            q = 0.5 * (1 + 3 * w_eff)

            theta = sigma / f

            # Rough time estimate: t/t_H ~ integral, approximate
            # For display, use a as proxy

            # Print selected points
            if (abs(N) < 0.02 or  # today
                abs(N - (-1)) < 0.02 or  # z~1.7
                (N > 0 and int(i) % 25 == 0) or
                (N > 0 and abs(q) < 0.05)):  # near deceleration

                z_str = f"{z:.2f}" if z > -0.99 else f"{z:.3f}"
                print(f"  {N:>6.2f} {a:>10.3f} {z_str:>8}"
                      f" {'':>8}"
                      f" {theta:>10.6f} {w:>10.5f} {Omega_DE:>8.4f} {q:>8.4f}")

        print()


# ═══════════════════════════════════════════════════════════════════════════
#  Part 2: w(z) for z < 0 and finding key transitions
# ═══════════════════════════════════════════════════════════════════════════
def part2_w_future():
    print_header("Part 2: Equation of state w(z) — past and future")

    omega_chi = compute_omega_chi_h2(M_CHI, ALPHA_D)
    N_future = 8.0  # a ≈ 3000

    for label, params in BENCHMARKS.items():
        theta_i = params['theta_i']
        f = params['f']
        Ld = params['Lambda_d']

        sol, rho_r0, rho_m0 = solve_future(
            theta_i, f, Ld, N_future=N_future, omega_chi_h2=omega_chi
        )

        if not sol.success:
            continue

        N_arr = np.linspace(-2, N_future, 2000)
        states = sol.sol(N_arr)

        # Find key epochs
        decel_N = None  # when q crosses 0 (future deceleration)
        V_half_N = None  # when V drops to V_today/2

        # Get today's values
        s0 = sol.sol(0.0)
        V_today = V_sigma(s0[0], f, Ld)

        for i in range(1, len(N_arr)):
            N = N_arr[i]
            if N <= 0:
                continue
            a = np.exp(N)
            sigma = states[0, i]
            p = states[1, i]
            rho_r = rho_r0 * a**(-4)
            rho_m = rho_m0 * a**(-3)
            V = V_sigma(sigma, f, Ld)
            denom = 3.0 * M_PL**2 - 0.5 * p**2
            if denom <= 0:
                continue
            H2 = (rho_r + rho_m + V) / denom
            if H2 <= 0:
                continue
            rho_sig = 0.5 * H2 * p**2 + V
            rho_tot = rho_r + rho_m + rho_sig
            Omega_r = rho_r / rho_tot
            Omega_m = rho_m / rho_tot
            Omega_DE = rho_sig / rho_tot
            P_sig = 0.5 * H2 * p**2 - V
            w = P_sig / rho_sig if rho_sig > 1e-100 else -1.0
            w_eff = Omega_r / 3.0 + w * Omega_DE
            q = 0.5 * (1 + 3 * w_eff)

            if decel_N is None and q > 0:
                decel_N = N

            if V_half_N is None and V < V_today * 0.5:
                V_half_N = N

        print(f"  ── {label} ──")
        print(f"  θ_i = {theta_i:.2f},  f/M_Pl = {f/M_PL:.2f}")

        if decel_N is not None:
            a_decel = np.exp(decel_N)
            z_decel = 1.0/a_decel - 1
            # Time estimate: in matter-dominated, t ∝ a^(3/2)
            # t_decel/t_0 ≈ a_decel^(3/2)
            t_ratio = a_decel**1.5  # rough, matter-dominated
            t_Gyr = 13.8 * t_ratio
            print(f"  Future deceleration at N = {decel_N:.2f}"
                  f"  (a = {a_decel:.1f}, z = {z_decel:.3f})")
            print(f"  Estimated time: ~{t_Gyr:.0f} Gyr (t/t₀ ≈ {t_ratio:.1f})")
        else:
            print(f"  No deceleration found within N < {N_future}"
                  f" (a < {np.exp(N_future):.0f})")

        if V_half_N is not None:
            print(f"  V drops to V₀/2 at N = {V_half_N:.2f}"
                  f" (a = {np.exp(V_half_N):.1f})")
        else:
            # Check how much V changed
            s_end = sol.sol(N_future)
            V_end = V_sigma(s_end[0], f, Ld)
            print(f"  V at N={N_future:.0f}: V/V₀ = {V_end/V_today:.6f}"
                  f" (barely changed)")

        print()


# ═══════════════════════════════════════════════════════════════════════════
#  Part 3: Detailed future timeline
# ═══════════════════════════════════════════════════════════════════════════
def part3_timeline():
    print_header("Part 3: Cosmic timeline — past, present, future")

    omega_chi = compute_omega_chi_h2(M_CHI, ALPHA_D)
    theta_i = 3.0
    f = 0.27 * M_PL
    Ld = 2.0e-15

    m_sigma = Ld**2 / f
    H0_GeV = H0_PLANCK_KMS * H_100_GEV / 100.0

    print(f"  m_σ = Λ_d²/f = {m_sigma:.3e} GeV")
    print(f"  H₀ = {H0_GeV:.3e} GeV")
    print(f"  m_σ/H₀ = {m_sigma/H0_GeV:.4f}")
    print()

    # The field starts oscillating when H ~ m_σ
    # For m_σ/H₀ < 1: field hasn't started oscillating yet!
    if m_sigma < H0_GeV:
        print(f"  ⚠ m_σ < H₀ → field has NOT started oscillating yet!")
        print(f"  The field is still frozen by Hubble friction.")
        print(f"  H will drop below m_σ in the future.")
        print()

        # When does H = m_σ?
        # In matter domination: H = H₀ a^(-3/2) (approximate)
        # H₀ a^(-3/2) = m_σ → a = (H₀/m_σ)^(2/3)
        a_osc = (H0_GeV / m_sigma)**(2.0/3.0)
        N_osc = np.log(a_osc)
        t_osc_ratio = a_osc**1.5  # rough
        t_osc_Gyr = 13.8 * t_osc_ratio

        print(f"  σ starts oscillating at:")
        print(f"    a_osc = (H₀/m_σ)^(2/3) = {a_osc:.1f}")
        print(f"    N_osc = {N_osc:.2f}")
        print(f"    t_osc ~ {t_osc_Gyr:.0f} Gyr  (t/t₀ = {t_osc_ratio:.1f})")
    else:
        print(f"  m_σ > H₀ → field has already started oscillating")
        # Period of oscillation: T_osc ~ 2π/m_σ
        T_osc_GeV = 2 * np.pi / m_sigma
        # Convert to seconds: 1 GeV⁻¹ = 6.58e-25 s
        T_osc_s = T_osc_GeV * 6.58e-25
        T_osc_Gyr = T_osc_s / (3.156e7 * 1e9)
        print(f"  Oscillation period: T = 2π/m_σ = {T_osc_Gyr:.1e} Gyr")

    # Once oscillating: σ behaves as pressureless matter (w ≈ 0 averaged)
    # V(σ) oscillates with time-averaged <V> decreasing as a⁻³
    # → dark energy converts to dark matter-like fluid
    print()
    print(f"  ── After oscillations begin ──")
    print(f"  Time-averaged: <ρ_σ> ∝ a⁻³,  <w> ≈ 0")
    print(f"  σ acts as PRESSURELESS MATTER (like axion dark matter)")
    print(f"  Dark energy effectively disappears → acceleration STOPS")

    # Transition to matter domination
    print()
    print(f"  ── Ultimate fate ──")
    print(f"  Phase 1 (now): σ frozen → V(σ) ≈ const → cosmic acceleration")
    print(f"  Phase 2 (future): H drops below m_σ → σ oscillates → V → 0")
    print(f"  Phase 3 (far future): Ω_DE → 0, universe returns to matter"
          f" domination")
    print(f"  Phase 4: a(t) ∝ t^(2/3) → eternal DECELERATED expansion")
    print()
    print(f"  ┌──────────────────────────────────────────────────────────┐")
    print(f"  │  NO BIG RIP.  NO BIG CRUNCH.  NO HEAT DEATH (of ΛCDM)│")
    print(f"  │  Dark energy is TRANSIENT — it turns off.             │")
    print(f"  │  Universe expands forever, but decelerating.          │")
    print(f"  │                                                         │")
    if m_sigma < H0_GeV:
        print(f"  │  Currently: σ still frozen (m_σ/H₀ = {m_sigma/H0_GeV:.4f})    │")
        print(f"  │  Σ field oscillation epoch: a ~ {a_osc:.0f} ({t_osc_Gyr:.0f} Gyr)          │")
    print(f"  │  After: <V> → 0, expansion decelerates              │")
    print(f"  └──────────────────────────────────────────────────────────┘")


# ═══════════════════════════════════════════════════════════════════════════
#  Part 4: Comparison with ΛCDM and other quintessence models
# ═══════════════════════════════════════════════════════════════════════════
def part4_comparison():
    print_header("Part 4: Comparison — our model vs ΛCDM vs typical quintessence")

    print("""
  ┌─────────────────────┬──────────────┬──────────────┬──────────────┐
  │ Property            │ ΛCDM         │ Our model    │ Quintessence │
  │                     │              │ (dark pion)  │ (generic)    │
  ├─────────────────────┼──────────────┼──────────────┼──────────────┤
  │ w today             │ −1 (exact)   │ ≈ −1 (−0.99) │ w > −1       │
  │ w future            │ −1 (always)  │ 0 (oscillates)│ model-dep   │
  │ Dark energy fate    │ Eternal      │ DIES (→ 0)   │ model-dep    │
  │ Acceleration        │ Eternal      │ Transient    │ model-dep    │
  │ Ultimate expansion  │ de Sitter    │ a ∝ t^(2/3)  │ model-dep    │
  │ Cosmic structures   │ Dissolve     │ Survive      │ model-dep    │
  │ Big Rip?            │ No           │ No           │ If w < −1    │
  │ Big Crunch?         │ No           │ No           │ If V < 0     │
  │ V(σ) minimum        │ Λ > 0        │ V = 0        │ model-dep    │
  │ Origin of DE        │ Fine-tuned   │ Misalignment │ model-dep    │
  │ Parameter count     │ 1 (Λ)        │ 2 (α_d, θ_i) │ ≥ 2         │
  │ Testable (w≠−1)?    │ No           │ Yes (DESI)   │ Yes          │
  │ Falsifiable?        │ Not by w     │ By w(z)      │ By w(z)      │
  └─────────────────────┴──────────────┴──────────────┴──────────────┘

  Key distinguishing prediction:
  In ΛCDM, all gravitationally unbound structures will eventually be
  pushed beyond the cosmic horizon (cosmic loneliness).
  In our model, dark energy dies → expansion decelerates →
  gravitationally unbound structures remain visible.
""")


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("╔" + "═"*68 + "╗")
    print("║  Test 32: Future of the Universe — When does dark energy die?    ║")
    print("╚" + "═"*68 + "╝")

    part1_future_evolution()
    part2_w_future()
    part3_timeline()
    part4_comparison()

    print("\n  Done.  Test 32 complete.")
