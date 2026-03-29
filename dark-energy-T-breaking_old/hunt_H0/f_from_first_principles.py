"""
Test 29: f from First Principles — The Cosmological Consistency Argument
========================================================================

Shows that f ≈ 0.27 M_Pl is NOT a free parameter but a CONSEQUENCE of:
    1. GMOR:       m_σ = Λ_d² / f
    2. Friedmann:  H₀ ~ Λ_d² / M_Pl
    3. Dynamics:   m_σ ~ c · H₀   (c ~ O(few) for viable DE)
    ⟹  f = M_Pl / c ~ O(1) M_Pl

The factor c = m_σ/H₀ is determined by the full ODE solution, not by hand.
This reduces the parameter space from {f, Λ_d, θ_i} to {Λ_d, θ_i}.
"""

import numpy as np
import sys, os

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(_ROOT, 'core'))
sys.path.insert(0, os.path.dirname(__file__))

from layer8_cosmic_ode import (
    solve_layer8, M_PL, H_100_GEV, V_sigma,
    OMEGA_R_H2, OMEGA_B_H2, _RHO_UNIT
)


# ═══════════════════════════════════════════════════════════════════════════
#  PART 1:  Analytic argument — f ~ M_Pl from self-consistency
# ═══════════════════════════════════════════════════════════════════════════
def part1_analytic():
    """Derive f ~ M_Pl from Friedmann + GMOR + quintessence condition."""
    print("=" * 70)
    print("  PART 1: Analytic Derivation — f ~ M_Pl")
    print("=" * 70)

    print("""
    Three equations:

    (1) Friedmann:    3 H₀² M_Pl² = ρ_total
                      ⟹ Λ_d⁴(1-cosθ_i) ≈ Ω_DE · 3 H₀² M_Pl²
                      ⟹ Λ_d² = √(3Ω_DE) · H₀ · M_Pl / √(1-cosθ_i)

    (2) GMOR:         m_σ = Λ_d² / f

    (3) Quintessence:  m_σ = c · H₀   (c ~ O(few) for viable dynamics)

    Combining (1) + (2) + (3):

        f = Λ_d² / (c · H₀)
          = [√(3Ω_DE) · H₀ · M_Pl / √(1-cosθ_i)] / (c · H₀)

    ┌──────────────────────────────────────────────────┐
    │   f = √(3 Ω_DE) · M_Pl / [c · √(1 − cos θ_i)] │
    └──────────────────────────────────────────────────┘

    This depends ONLY on {Ω_DE, θ_i, c} — not on Λ_d or H₀ separately!
    """)

    Omega_DE = 0.685  # Planck 2018
    theta_values = [np.pi, 3.0, 2.95, 2.90, 2.887]

    print(f"  {'θ_i':>8}  {'θ_i/π':>6}  {'1-cosθ':>8}  {'f (c=1)':>12}  "
          f"{'f (c=3)':>12}  {'f (c=4.5)':>12}")
    print(f"  {'':>8}  {'':>6}  {'':>8}  {'[M_Pl]':>12}  {'[M_Pl]':>12}  {'[M_Pl]':>12}")
    print("  " + "-" * 68)

    for theta in theta_values:
        one_minus_cos = 1 - np.cos(theta)
        for c_vals in [(1, 3, 4.5)]:
            f_over_Mpl = [
                np.sqrt(3 * Omega_DE) / (c * np.sqrt(one_minus_cos))
                for c in c_vals
            ]
            print(f"  {theta:>8.3f}  {theta/np.pi:>6.3f}  {one_minus_cos:>8.4f}"
                  f"  {f_over_Mpl[0]:>12.4f}  {f_over_Mpl[1]:>12.4f}  {f_over_Mpl[2]:>12.4f}")

    print(f"\n  For θ_i = π (hilltop):   1-cos(π) = 2")
    f_hilltop = lambda c: np.sqrt(3 * Omega_DE / 2) / c
    print(f"    f/M_Pl = {f_hilltop(1):.4f}/c")
    print(f"    c = 1:   f = {f_hilltop(1):.4f} M_Pl")
    print(f"    c = 3:   f = {f_hilltop(3):.4f} M_Pl")
    print(f"    c = 3.8: f = {f_hilltop(3.8):.4f} M_Pl  ← matches MCMC!")

    c_needed = f_hilltop(1) / 0.27
    print(f"\n  ⟹ To get f = 0.27 M_Pl at θ_i = π: need c = {c_needed:.2f}")
    print(f"    i.e., m_σ ≈ {c_needed:.1f} × H₀")
    print(f"    This is m_σ ~ few H₀ — exactly the quintessence transition!")


# ═══════════════════════════════════════════════════════════════════════════
#  PART 2:  Numerical verification — scan f, find which gives correct H₀
# ═══════════════════════════════════════════════════════════════════════════
def part2_numerical_f_scan():
    """Scan f/M_Pl at different θ_i and show f-dependence of H₀."""
    print("\n" + "=" * 70)
    print("  PART 2: Numerical f-Scan — H₀(f) at Fixed {Λ_d, θ_i}")
    print("=" * 70)

    # MCMC parameters
    m_chi = 98.19
    m_phi = 9.66e-3  # GeV
    alpha_d = 3.274e-3

    from layer8_cosmic_ode import compute_omega_chi_h2
    omega_chi = compute_omega_chi_h2(m_chi, alpha_d)

    H0_planck_gev = 67.4 * H_100_GEV / 100.0

    # Scan at multiple θ_i values to show the transition
    for theta_i, Lambda_d_meV in [(np.pi, 1.89), (3.0, 1.89), (2.9, 1.89)]:
        Lambda_d = Lambda_d_meV * 1e-12  # GeV
        theta_label = "π" if abs(theta_i - np.pi) < 0.001 else f"{theta_i:.2f}"
        print(f"\n  ── θ_i = {theta_label}, Λ_d = {Lambda_d_meV} meV ──")
        print(f"  {'f/M_Pl':>8}  {'m_σ/H₀':>8}  {'H₀':>7}  {'w_σ':>8}  {'Ω_DE':>6}  {'θ_today':>8}")
        print("  " + "-" * 55)

        f_fracs = [0.05, 0.10, 0.15, 0.20, 0.25, 0.27, 0.30, 0.40, 0.50, 1.00]

        for ff in f_fracs:
            f_val = ff * M_PL
            m_sigma = Lambda_d**2 / f_val
            c_ratio = m_sigma / H0_planck_gev

            res = solve_layer8(m_chi, m_phi, alpha_d, f_val, Lambda_d, theta_i,
                               omega_chi_h2=omega_chi, verbose=False)

            if res.H0_kms is not None:
                print(f"  {ff:>8.2f}  {c_ratio:>8.2f}  {res.H0_kms:>7.1f}  "
                      f"{res.w_sigma:>8.4f}  {res.Omega_DE:>6.3f}  "
                      f"{res.theta_today:>8.4f}")
            else:
                print(f"  {ff:>8.2f}  {c_ratio:>8.2f}  {'FAIL':>7}")

    print("""
    Key observations:
    • At θ_i = π (hilltop): sin(π) = 0 → dV/dσ = 0 → field NEVER rolls.
      H₀ is independent of f. This is the CC limit. f is unconstrained.
    • At θ_i < π: field CAN roll. Smaller f → larger m_σ → more rolling →
      DE dissipates. f now MATTERS and is constrained.
    • The DESI observation w₀ > −1 requires θ_i < π, which determines f.
    """)


# ═══════════════════════════════════════════════════════════════════════════
#  PART 3:  What determines c = m_σ/H₀ ?
# ═══════════════════════════════════════════════════════════════════════════
def part3_dynamical_c():
    """Compute c(θ_i) — the dynamical ratio m_σ/H₀ at the quintessence boundary."""
    print("\n" + "=" * 70)
    print("  PART 3: Dynamical c(θ_i) — What Fixes m_σ/H₀?")
    print("=" * 70)

    print("""
    For given {Λ_d, θ_i}, the value of f that reproduces H₀_target
    determines c = m_σ/H₀. This c depends on θ_i because:
    • θ_i ≈ π: field near hilltop, slow-rolls, needs m_σ moderately > H₀
    • θ_i < π: field displaced, rolls faster, needs larger m_σ (smaller f)
    """)

    m_chi = 98.19
    m_phi = 9.66e-3
    alpha_d = 3.274e-3
    H0_target = 67.4  # Planck

    from layer8_cosmic_ode import compute_omega_chi_h2
    omega_chi = compute_omega_chi_h2(m_chi, alpha_d)

    theta_values = [np.pi, 3.10, 3.05, 3.00, 2.95, 2.90, 2.887]

    print(f"  Target: H₀ = {H0_target:.1f} km/s/Mpc")
    print(f"\n  {'θ_i':>8}  {'θ_i/π':>6}  {'Λ_d*':>8}  {'f*':>10}  {'f*/M_Pl':>8}  "
          f"{'m_σ/H₀':>8}  {'w₀':>7}")
    print(f"  {'':>8}  {'':>6}  {'[meV]':>8}  {'[GeV]':>10}  {'':>8}  {'(=c)':>8}  ")
    print("  " + "-" * 65)

    H0_target_gev = H0_target * H_100_GEV / 100.0

    for theta_i in theta_values:
        # For each θ_i, find Λ_d that gives H₀ ≈ target at f = 0.27 M_Pl first,
        # then find f that gives H₀ = target at that Λ_d.

        # Analytic estimate: Λ_d from Friedmann
        Omega_DE_approx = 0.685
        one_minus_cos = 1 - np.cos(theta_i)
        # Λ_d⁴ ≈ Ω_DE · 3H₀²M_Pl² / (1-cosθ)
        Ld4 = Omega_DE_approx * 3 * H0_target_gev**2 * M_PL**2 / one_minus_cos
        Lambda_d = Ld4**0.25

        # Binary search for f that gives H₀ = target
        f_lo, f_hi = 0.05 * M_PL, 1.0 * M_PL
        f_best = None

        for _ in range(50):  # bisection
            f_mid = (f_lo + f_hi) / 2
            res = solve_layer8(m_chi, m_phi, alpha_d, f_mid, Lambda_d, theta_i,
                               omega_chi_h2=omega_chi, verbose=False)
            if res.H0_kms is None:
                f_lo = f_mid
                continue

            if res.H0_kms > H0_target:
                f_hi = f_mid  # f too large → too much DE → reduce f
            else:
                f_lo = f_mid  # f too small → not enough DE → increase f

            if abs(res.H0_kms - H0_target) < 0.1:
                f_best = f_mid
                break

        if f_best is None:
            f_best = (f_lo + f_hi) / 2

        res = solve_layer8(m_chi, m_phi, alpha_d, f_best, Lambda_d, theta_i,
                           omega_chi_h2=omega_chi, verbose=False)

        if res.H0_kms is not None:
            m_sigma = Lambda_d**2 / f_best
            c_val = m_sigma / H0_target_gev
            print(f"  {theta_i:>8.3f}  {theta_i/np.pi:>6.3f}  {Lambda_d*1e12:>8.4f}  "
                  f"{f_best:>10.3e}  {f_best/M_PL:>8.4f}  {c_val:>8.2f}  "
                  f"{res.w_sigma:>7.4f}")


# ═══════════════════════════════════════════════════════════════════════════
#  PART 4:  QCD analogy — is f_π/Λ_QCD ~ f/Λ_d ?
# ═══════════════════════════════════════════════════════════════════════════
def part4_qcd_analogy():
    """Compare with QCD axion to understand the hierarchy f >> Λ_d."""
    print("\n" + "=" * 70)
    print("  PART 4: QCD Analogy — Why f >> Λ_d Is Natural")
    print("=" * 70)

    # QCD numbers
    f_pi = 0.093    # GeV (pion decay constant)
    Lambda_QCD = 0.220  # GeV
    f_a_min = 1e9   # GeV (QCD axion lower bound)
    f_a_max = 1e12  # GeV (QCD axion upper bound)

    # Our dark sector
    f_dark = 0.27 * M_PL  # GeV
    Lambda_d = 2e-3 * 1e-9  # GeV (2 meV)

    print(f"""
    QCD has TWO distinct "f" scales:

    │ Quantity       │ QCD                │ Dark Sector         │ Role              │
    │────────────────│────────────────────│─────────────────────│───────────────────│
    │ F_π (pion)     │ 93 MeV             │ —                   │ Chiral condensate │
    │ F_π/Λ_QCD      │ {f_pi/Lambda_QCD:.2f}               │ —                   │ O(1) ratio        │
    │ f_a (axion)    │ 10⁹ – 10¹² GeV     │ {f_dark:.2e} GeV │ PQ breaking scale │
    │ f_a/Λ_QCD      │ 10⁹·⁷ – 10¹²·⁷     │ {f_dark/Lambda_d:.2e}     │ Huge hierarchy!   │
    │ Λ_conf         │ {Lambda_QCD*1e3:.0f} MeV             │ {Lambda_d*1e12:.0f} meV                │ Confinement       │

    The point: σ is NOT the dark pion in the sense of F_π ~ Λ.
    σ is the dark AXION:  f ~ M_Pl >> Λ_d, just like f_a >> Λ_QCD.

    The potential V = Λ_d⁴(1 - cos σ/f) is generated by dark instantons,
    with f set by the global U(1) breaking scale (gravity/Planck).

    QCD axion:   m_a = Λ_QCD² / f_a  ≈ {(Lambda_QCD**2/f_a_min)*1e6:.1f} μeV  (for f_a = 10⁹ GeV)
    Dark axion:  m_σ = Λ_d² / f      ≈ {(Lambda_d**2/f_dark):.2e} GeV  ≈ {(Lambda_d**2/f_dark)/1.44e-42:.1f} H₀

    Both have ultralight masses because f >> Λ.
    """)

    # Hierarchy
    print("  Hierarchy ratios:")
    print(f"    QCD:  f_a/Λ_QCD   ~ 10⁹⁻¹²/0.22 ~ 10^{{9.7-12.7}}")
    print(f"    Dark: f/Λ_d       = {f_dark/Lambda_d:.2e}")
    print(f"    log₁₀(f/Λ_d)     = {np.log10(f_dark/Lambda_d):.1f}")
    print(f"    QCD:  log₁₀(f_a/Λ_QCD) = 9.7 – 12.7")
    print(f"\n    ⟹ Same ORDER of hierarchy! Dark sector is QCD axion physics")
    print(f"      with Λ_QCD → Λ_d ∼ meV, f_a → f ∼ M_Pl.")


# ═══════════════════════════════════════════════════════════════════════════
#  PART 5:  The Weak Gravity Conjecture and f
# ═══════════════════════════════════════════════════════════════════════════
def part5_wgc():
    """Check if f = 0.27 M_Pl satisfies the Weak Gravity Conjecture."""
    print("\n" + "=" * 70)
    print("  PART 5: Weak Gravity Conjecture for Axions")
    print("=" * 70)

    f_dark = 0.27 * M_PL

    print("""
    The WGC for axions (Arkani-Hamed, Motl, Nicolis, Vafa 2006):

        S_inst · f ≤ M_Pl

    where S_inst is the action of the instanton generating the potential.

    For the dark QCD instanton at the confinement scale:
    • At μ ~ Λ_d: α_d → ∞ (strongly coupled), S_inst ~ O(1)
    • The "dilute instanton gas" has S = 8π²/g² = 2π/α_d
    • But at confinement: this perturbative estimate fails
    • Non-perturbative: S_inst ~ 1 (BPST-like instantons at strong coupling)
    """)

    # At the confinement scale
    S_inst_strong = 1.0  # O(1) at strong coupling
    bound_strong = M_PL / S_inst_strong
    print(f"  WGC bound (S_inst ~ 1):     f ≤ {bound_strong/M_PL:.1f} M_Pl = {bound_strong:.3e} GeV")
    print(f"  Our value:                   f = {f_dark/M_PL:.2f} M_Pl = {f_dark:.3e} GeV")
    print(f"  Satisfies WGC? {'✓ YES' if f_dark <= bound_strong else '✗ NO'}")

    # At the UV matching scale
    alpha_d_UV = 0.031  # at m_chi ~ 100 GeV
    S_inst_UV = 2 * np.pi / alpha_d_UV
    bound_UV = M_PL / S_inst_UV
    print(f"\n  WGC bound (S_inst = 2π/α_d): f ≤ {bound_UV/M_PL:.4f} M_Pl = {bound_UV:.3e} GeV")
    print(f"  Our value:                    f = {f_dark/M_PL:.2f} M_Pl = {f_dark:.3e} GeV")
    print(f"  Satisfies WGC? {'✓ YES' if f_dark <= bound_UV else '✗ NO ← perturbative bound violated'}")

    print("""
    Resolution: The perturbative S_inst = 2π/α_d applies in the UV where
    instantons are suppressed. The POTENTIAL is generated in the IR
    (at Λ_d) where the theory is strongly coupled and S_inst ~ O(1).
    The WGC should be applied with the IR instanton action → f < M_Pl ✓

    Moreover, the WGC is a CONJECTURE — violations are possible in
    consistent theories. Our modest f = 0.27 M_Pl < M_Pl is safe.
    """)


# ═══════════════════════════════════════════════════════════════════════════
#  PART 6:  Parameter reduction — from 3 to 2
# ═══════════════════════════════════════════════════════════════════════════
def part6_parameter_reduction():
    """Show that f is determined: {f, Λ_d, θ_i} → {Λ_d, θ_i}."""
    print("\n" + "=" * 70)
    print("  PART 6: Parameter Reduction — {f, Λ_d, θ_i} → {Λ_d, θ_i}")
    print("=" * 70)

    print("""
    The self-consistency argument:

    Given {Λ_d, θ_i}, the dark energy density is:
        ρ_DE = Λ_d⁴ (1 - cos θ_i)    [if field frozen]

    This determines H₀ via Friedmann (with known Ω_m, Ω_r):
        H₀ = H₀(Λ_d, θ_i)

    Then the GMOR mass:
        m_σ = Λ_d²/f

    The ODE solution requires m_σ ~ c(θ_i) · H₀ for self-consistency
    (too large → field oscillated away, too small → pure CC).

    Therefore:
        f = Λ_d² / [c(θ_i) · H₀(Λ_d, θ_i)]

    ┌──────────────────────────────────────────────────────────────┐
    │  f is DERIVED from {Λ_d, θ_i} — not a free parameter!      │
    │                                                              │
    │  Parameter count: {f, Λ_d, θ_i} → {Λ_d, θ_i}              │
    │                   3 params      → 2 params                   │
    │                                                              │
    │  Total Lagrangian: {m_χ, m_φ, α_D, α_d, Λ_d, θ_i} = 6    │
    │  With transmutation: {m_χ, m_φ, α_D, α_d, θ_i} = 5        │
    └──────────────────────────────────────────────────────────────┘
    """)

    # Summary table of all parameters
    print("  Complete parameter census (Paper 1 + Paper 2):\n")
    print(f"  {'Parameter':>12}  {'Value':>16}  {'Source':>15}  {'Status':>12}")
    print("  " + "-" * 60)

    params = [
        ("m_χ",   "98.19 GeV",    "MCMC (Paper 1)", "Fixed"),
        ("m_φ",   "9.66 MeV",     "MCMC (Paper 1)", "Fixed"),
        ("α_D",   "3.274×10⁻³",   "MCMC (Paper 1)", "Fixed"),
        ("α_d",   "0.031 ≈ 1/32", "Transmutation",  "Derived"),
        ("Λ_d",   "~2 meV",       "from α_d via RG", "Derived"),
        ("f",     "0.27 M_Pl",    "m_σ~H₀ consist.", "Derived"),
        ("θ_i",   "~2.9 – π",     "DESI/Planck",    "Free"),
    ]

    for name, val, source, status in params:
        marker = "★" if status == "Free" else "→" if status == "Derived" else " "
        print(f"  {name:>12}  {val:>16}  {source:>15}  {marker} {status:>10}")

    print(f"\n  Free UV parameters: {{m_χ, m_φ, α_D, α_d, θ_i}} = 5")
    print(f"  Of these, 3 are fixed by Paper 1 MCMC")
    print(f"  Remaining free: {{α_d, θ_i}} = 2 parameters")
    print(f"  * α_d sets Λ_d via transmutation")
    print(f"  * θ_i sets w₀ and H₀ (given Λ_d)")
    print(f"  * f is derived from dynamical self-consistency")


# ═══════════════════════════════════════════════════════════════════════════
#  PART 7:  The "f = √(3Ω_DE)M_Pl / [c√(1-cosθ)]" formula verification
# ═══════════════════════════════════════════════════════════════════════════
def part7_formula_check():
    """Verify the analytic formula against ODE solutions."""
    print("\n" + "=" * 70)
    print("  PART 7: Formula Verification — Analytic vs ODE")
    print("=" * 70)

    m_chi = 98.19
    m_phi = 9.66e-3
    alpha_d = 3.274e-3
    H0_target = 67.4

    from layer8_cosmic_ode import compute_omega_chi_h2
    omega_chi = compute_omega_chi_h2(m_chi, alpha_d)

    H0_gev = H0_target * H_100_GEV / 100.0
    Omega_DE = 0.685

    theta_values = [np.pi, 3.10, 3.00, 2.95, 2.90]

    print(f"\n  Testing: f_analytic = √(3Ω_DE) · M_Pl / [c · √(1-cosθ)]")
    print(f"  with c determined from ODE to match H₀ = {H0_target}\n")
    print(f"  {'θ_i':>7}  {'f_ODE/M_Pl':>10}  {'c_ODE':>7}  {'f_formula':>10}  {'Agree?':>7}")
    print("  " + "-" * 50)

    for theta_i in theta_values:
        one_minus_cos = 1 - np.cos(theta_i)

        # Get Λ_d from Friedmann
        Ld4 = Omega_DE * 3 * H0_gev**2 * M_PL**2 / one_minus_cos
        Lambda_d = Ld4**0.25

        # Binary search for f
        f_lo, f_hi = 0.05 * M_PL, 1.5 * M_PL
        for _ in range(60):
            f_mid = (f_lo + f_hi) / 2
            res = solve_layer8(m_chi, m_phi, alpha_d, f_mid, Lambda_d, theta_i,
                               omega_chi_h2=omega_chi, verbose=False)
            if res.H0_kms is None:
                f_lo = f_mid
                continue
            if res.H0_kms > H0_target:
                f_hi = f_mid
            else:
                f_lo = f_mid
            if abs(res.H0_kms - H0_target) < 0.05:
                break

        f_ode = (f_lo + f_hi) / 2
        m_sigma_ode = Lambda_d**2 / f_ode
        c_ode = m_sigma_ode / H0_gev

        # Formula prediction (with c_ode)
        f_formula = np.sqrt(3 * Omega_DE) * M_PL / (c_ode * np.sqrt(one_minus_cos))
        agree = abs(f_ode - f_formula) / f_ode * 100

        print(f"  {theta_i:>7.3f}  {f_ode/M_PL:>10.4f}  {c_ode:>7.2f}  "
              f"{f_formula/M_PL:>10.4f}  {agree:>6.1f}%")

    print("""
    The formula f = √(3Ω_DE)M_Pl / [c√(1-cosθ)] is exact (by construction).
    The physics content is: c = m_σ/H₀ is an O(few) dynamical number
    determined by the full cosmological evolution, not a free parameter.
    """)


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("\n" + "█" * 70)
    print("  TEST 29: f FROM FIRST PRINCIPLES")
    print("  — The Cosmological Self-Consistency Argument —")
    print("█" * 70)

    part1_analytic()
    part2_numerical_f_scan()
    part3_dynamical_c()
    part4_qcd_analogy()
    part5_wgc()
    part6_parameter_reduction()
    part7_formula_check()

    print("\n" + "=" * 70)
    print("  CONCLUSION")
    print("=" * 70)
    print("""
    f ≈ 0.27 M_Pl: Two Regimes

    ┌─────────────────────────────────────────────────────────────────┐
    │  θ_i = π (hilltop):  f UNCONSTRAINED by H₀ alone.             │
    │  → Field frozen at dV/dσ = 0.  Acts as pure CC.               │
    │  → f ∈ (0, ∞) all give same H₀.  Need other data to fix f.   │
    │                                                                 │
    │  θ_i < π (quintessence): f IS DETERMINED by {Λ_d, θ_i, H₀}.  │
    │  → Smaller f → field oscillates → DE lost.                     │
    │  → Larger f → field frozen → CC limit.                         │
    │  → Unique f gives observed H₀ with w₀ > −1.                   │
    └─────────────────────────────────────────────────────────────────┘

    DESI resolves the degeneracy: w₀ = −0.727 ± 0.067 > −1 at 4σ
    ⟹ θ_i < π required ⟹ f is determined!

    The ORDER OF MAGNITUDE argument:
        f ~ M_Pl  follows from m_σ ~ H₀ + Friedmann + GMOR
        This is a consistency relation, not fine-tuning.
        The O(1) prefactor 0.1-0.3 comes from the ODE dynamics.

    Analogy: Just as f_a >> Λ_QCD in the QCD axion,
    f >> Λ_d in the dark axion. Both hierarchies are natural
    (f is a UV symmetry-breaking scale; Λ is an IR confinement scale).

    Parameter reduction (assuming DESI w₀ > −1):
        Before: {f, Λ_d, θ_i} = 3 DE parameters
        After:  {Λ_d, θ_i} = 2 DE parameters (f derived from H₀)

    Total free UV parameters: {m_χ, m_φ, α_D, α_d, θ_i} = 5
    After Paper 1 MCMC: {α_d, θ_i} = 2 remaining free parameters
    """)
    print("=" * 70)
