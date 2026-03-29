"""
RG Transmutation:  α_D → Λ_d → m_σ → H₀
==========================================

Path 2 from the hunt-H₀ roadmap.

Physics:  dimensional transmutation in dark SU(2)
    Λ_d = μ · exp( −2π / (b₀ · α_d(μ)) )

where b₀ is the 1-loop beta-function coefficient of the dark gauge group.

For SU(2)_d with N_f Majorana fermions in the fundamental:
    b₀ = (11/3)·C₂(G) − (2/3)·T(R)·N_f
       = (11/3)·2 − (2/3)·(1/2)·N_f
       = 22/3 − N_f/3

1-loop RG running (from μ₁ to μ₂):
    1/α(μ₂) = 1/α(μ₁) + (b₀/(2π)) · ln(μ₂/μ₁)

Two key questions:
    Q1: Given α_D(m_χ) from MCMC, what Λ_d does transmutation predict?
    Q2: What α_D(m_χ) would we need to get Λ_d ∼ 2 meV?
    Q3: Can the RG running between m_φ (9.66 MeV) and m_χ (98 GeV)
        reconcile any tension?
"""

import numpy as np
import math
import sys
import os

# ── Constants (same as layer8_cosmic_ode.py) ─────────────────────────────
M_PL       = 2.435e18       # reduced Planck mass [GeV]
H_100_GEV  = 2.1332e-42     # 100 km/s/Mpc [GeV]
H0_PLANCK  = 67.4 * H_100_GEV / 100.0  # Planck H₀ [GeV]

# ── MCMC best-fit parameters ────────────────────────────────────────────
M_CHI      = 98.19           # GeV
M_PHI      = 9.66e-3         # GeV (dark mediator)
ALPHA_D    = 3.274e-3        # dark fine-structure constant (MCMC)

# ── Dark QCD sector ─────────────────────────────────────────────────────
LAMBDA_D_TARGET = 2.0e-12    # GeV  (2 meV — needed for dark energy)
F_AXION         = 0.27 * M_PL  # dark pion decay constant


# ═══════════════════════════════════════════════════════════════════════════
#  Beta-function coefficient
# ═══════════════════════════════════════════════════════════════════════════
def beta_coeff_SU2(N_f_majorana):
    """
    1-loop beta function coefficient for SU(2) with N_f Majorana fundamentals.
    
    β(α) = −b₀ α² / (2π)
    
    b₀ = (11/3)·C₂(G) − (2/3)·T(R)·N_f
    
    SU(2): C₂(G) = 2, T(fund) = 1/2
    Each Majorana = 1/2 Dirac, so coefficient is (2/3)·(1/2)·N_f for N_f
    Majorana fermions, where the standard formula uses Dirac fermions.
    
    More precisely: for Majorana fermions, each contributes (1/3)·T(R) to
    the coefficient, while Dirac contributes (4/3)·T(R).  A Majorana
    fundamental of SU(2) contributes (1/3)·(1/2) = 1/6 per species.
    
    Actually let's be precise:
    General formula: b₀ = (11/3)C₂(G) - (4/3)T(R)·n_D - (2/3)T(R)·n_M
    where n_D = number of Dirac fermions, n_M = number of Majorana fermions.
    
    For SU(2) with N_f Majorana quarks (no Dirac):
    b₀ = (11/3)·2 - (2/3)·(1/2)·N_f = 22/3 - N_f/3
    """
    C2_G = 2.0       # SU(2) Casimir
    T_R  = 0.5       # fundamental rep
    return (11.0/3.0) * C2_G - (2.0/3.0) * T_R * N_f_majorana


def asymptotic_freedom_limit(N_f_majorana=None):
    """Maximum N_f for asymptotic freedom (b₀ > 0)."""
    # b₀ = 22/3 - N_f/3 > 0  →  N_f < 22
    if N_f_majorana is not None:
        return beta_coeff_SU2(N_f_majorana) > 0
    return 22  # max Majorana species for SU(2)_d asymptotic freedom


# ═══════════════════════════════════════════════════════════════════════════
#  1-loop RG running
# ═══════════════════════════════════════════════════════════════════════════
def alpha_running(alpha_mu1, mu1, mu2, b0):
    """
    1-loop RG evolution:  1/α(μ₂) = 1/α(μ₁) + (b₀/(2π))·ln(μ₂/μ₁)
    
    Returns α(μ₂). Returns None if Landau pole encountered.
    """
    inv_alpha = 1.0/alpha_mu1 + (b0 / (2*math.pi)) * math.log(mu2 / mu1)
    if inv_alpha <= 0:
        return None  # Landau pole (shouldn't happen for AF theory going UP)
    return 1.0 / inv_alpha


def transmutation_scale(alpha_mu, mu, b0):
    """
    Dimensional transmutation:
        Λ = μ · exp( −2π / (b₀ · α(μ)) )
    
    This is RG-invariant (independent of choice of μ, to 1-loop accuracy).
    """
    exponent = -2.0 * math.pi / (b0 * alpha_mu)
    # Use log to check for underflow
    log_Lambda = math.log(mu) + exponent
    return math.exp(log_Lambda)


def alpha_needed_for_Lambda(Lambda_d, mu, b0):
    """
    Invert transmutation: what α(μ) gives a desired Λ_d?
    
        α(μ) = 2π / (b₀ · ln(μ/Λ_d))
    """
    return 2.0 * math.pi / (b0 * math.log(mu / Lambda_d))


# ═══════════════════════════════════════════════════════════════════════════
#  Main analysis
# ═══════════════════════════════════════════════════════════════════════════
def main():
    N_f = 3  # 3 Majorana dark quarks (from model: χ₁, χ₂, χ₃ under A₄)
    b0 = beta_coeff_SU2(N_f)
    
    print("=" * 78)
    print("  Path 2: RG Transmutation — Λ_d from α_D")
    print("=" * 78)
    
    # ── Dark gauge group parameters ──────────────────────────────────────
    print(f"\n  Dark gauge group:  SU(2)_d")
    print(f"  Majorana quarks:   N_f = {N_f}  (A₄ triplet: χ₁, χ₂, χ₃)")
    print(f"  b₀ = (11/3)·2 − (2/3)·(1/2)·{N_f} = {b0:.4f}  ({22-N_f}/3)")
    print(f"  Asymptotic freedom: {'Yes ✓' if b0 > 0 else 'No ✗'}")
    print(f"  (Max N_f for AF: {asymptotic_freedom_limit()})")
    
    # ── Input parameters ─────────────────────────────────────────────────
    print(f"\n  ── MCMC best-fit ──")
    print(f"  m_χ   = {M_CHI:.2f} GeV")
    print(f"  m_φ   = {M_PHI*1e3:.2f} MeV")
    print(f"  α_D   = {ALPHA_D:.4e}")
    print(f"  1/α_D = {1/ALPHA_D:.2f}")
    
    print(f"\n  ── Target ──")
    print(f"  Λ_d (needed) = {LAMBDA_D_TARGET:.1e} GeV  = {LAMBDA_D_TARGET*1e12:.1f} meV")
    print(f"  f             = {F_AXION:.3e} GeV  = {F_AXION/M_PL:.2f} M_Pl")
    print(f"  m_σ = Λ_d²/f  = {LAMBDA_D_TARGET**2/F_AXION:.3e} GeV")
    
    # ═════════════════════════════════════════════════════════════════════
    #  Analysis 1: Forward — α_D(m_χ) → Λ_d
    # ═════════════════════════════════════════════════════════════════════
    print(f"\n{'='*78}")
    print(f"  Analysis 1:  Forward transmutation at μ = m_χ")
    print(f"{'='*78}")
    
    Lambda_from_mchi = transmutation_scale(ALPHA_D, M_CHI, b0)
    print(f"\n  Λ_d = m_χ · exp(−2π/(b₀·α_D)) ")
    print(f"      = {M_CHI} · exp(−2π/({b0:.4f} × {ALPHA_D:.4e}))")
    exponent_val = -2*math.pi/(b0 * ALPHA_D)
    print(f"      = {M_CHI} · exp({exponent_val:.1f})")
    print(f"      = {Lambda_from_mchi:.3e} GeV")
    print(f"\n  ⚠ log₁₀(Λ_d/GeV) = {math.log10(Lambda_from_mchi):.1f}")
    print(f"  ⚠ This is absurdly small: ~10^{math.log10(Lambda_from_mchi):.0f} GeV")
    print(f"  ⚠ Target was: Λ_d = 2×10⁻¹² GeV (2 meV)")
    print(f"\n  ❌ MCMC α_D = {ALPHA_D:.4e} is too small for transmutation!")
    
    # ═════════════════════════════════════════════════════════════════════
    #  Analysis 2: Inverse — what α_D(m_χ) gives Λ_d = 2 meV?
    # ═════════════════════════════════════════════════════════════════════
    print(f"\n{'='*78}")
    print(f"  Analysis 2:  Inverse — what α_D(m_χ) gives Λ_d = 2 meV?")
    print(f"{'='*78}")
    
    alpha_needed = alpha_needed_for_Lambda(LAMBDA_D_TARGET, M_CHI, b0)
    print(f"\n  α_D(m_χ) needed = 2π / (b₀ · ln(m_χ/Λ_d))")
    print(f"                   = 2π / ({b0:.4f} · ln({M_CHI}/{LAMBDA_D_TARGET:.1e}))")
    print(f"                   = 2π / ({b0:.4f} · {math.log(M_CHI/LAMBDA_D_TARGET):.2f})")
    print(f"                   = {alpha_needed:.5f}")
    print(f"                   ≈ 1/{1/alpha_needed:.1f}")
    
    ratio = alpha_needed / ALPHA_D
    print(f"\n  Comparison with MCMC:")
    print(f"    α_D(MCMC)    = {ALPHA_D:.5f}     (1/{1/ALPHA_D:.1f})")
    print(f"    α_D(needed)  = {alpha_needed:.5f}   (1/{1/alpha_needed:.1f})")
    print(f"    ratio        = {ratio:.2f}×")
    print(f"\n  Needed coupling is {ratio:.1f}× larger than MCMC value")
    
    # ═════════════════════════════════════════════════════════════════════
    #  Analysis 3: RG running from m_φ to m_χ
    #  
    #  Could α_D(m_φ) be larger? In the broken phase below m_φ, the gauge
    #  boson is massive and the running changes. But above m_φ (in the
    #  unbroken SU(2)_d phase), the coupling runs with the full b₀.
    #
    #  The MCMC α_D enters at the m_φ scale (mediator exchange), so:
    #     α_D = α_d(m_φ)
    #  Then running UP from m_φ to M_Pl:
    #     1/α_d(M_Pl) = 1/α_d(m_φ) + (b₀/2π)·ln(M_Pl/m_φ)
    # ═════════════════════════════════════════════════════════════════════
    print(f"\n{'='*78}")
    print(f"  Analysis 3:  RG running from m_φ to M_Pl")
    print(f"{'='*78}")
    
    # The scale hierarchy: m_φ (9.66 MeV) → m_χ (98 GeV) → M_Pl (2.44×10¹⁸ GeV)
    # Possibility: α_D is defined at m_φ scale
    alpha_mpl = alpha_running(ALPHA_D, M_PHI, M_PL, b0)
    Lambda_from_mpl = transmutation_scale(alpha_mpl, M_PL, b0) if alpha_mpl else None
    
    print(f"\n  Running α_d from m_φ = {M_PHI*1e3:.2f} MeV to M_Pl = {M_PL:.3e} GeV:")
    print(f"    α_d(m_φ) = {ALPHA_D:.5f}  (= α_D from MCMC)")
    print(f"    1/α_d(m_φ) = {1/ALPHA_D:.2f}")
    log_ratio = math.log(M_PL / M_PHI)
    delta_inv = (b0 / (2*math.pi)) * log_ratio
    print(f"    Δ(1/α) = (b₀/2π)·ln(M_Pl/m_φ) = ({b0:.4f}/2π)·{log_ratio:.2f} = {delta_inv:.2f}")
    if alpha_mpl:
        print(f"    1/α_d(M_Pl) = {1/ALPHA_D:.2f} + {delta_inv:.2f} = {1/alpha_mpl:.2f}")
        print(f"    α_d(M_Pl) = {alpha_mpl:.6e}")
        print(f"\n  Transmutation from M_Pl:")
        print(f"    Λ_d = M_Pl · exp(−2π/(b₀·α_d(M_Pl)))")
        exp_val = -2*math.pi/(b0 * alpha_mpl)
        print(f"        = {M_PL:.3e} · exp({exp_val:.1f})")
        if Lambda_from_mpl and Lambda_from_mpl > 0:
            print(f"        = {Lambda_from_mpl:.3e} GeV")
            print(f"        = {Lambda_from_mpl*1e12:.3e} meV")
        else:
            print(f"        ≈ 10^{exp_val/math.log(10) + math.log10(M_PL):.0f} GeV")
    
    print(f"\n  Note: transmutation is RG-invariant at 1-loop.")
    print(f"  Λ_d(from m_χ) = Λ_d(from m_φ) = Λ_d(from M_Pl) by definition.")
    
    # ═════════════════════════════════════════════════════════════════════
    #  Analysis 4: What α_d(m_φ) gives Λ_d = 2 meV? (at m_φ scale)
    # ═════════════════════════════════════════════════════════════════════
    print(f"\n{'='*78}")
    print(f"  Analysis 4:  What α_d(m_φ) gives Λ_d = 2 meV?")
    print(f"{'='*78}")
    
    alpha_needed_mphi = alpha_needed_for_Lambda(LAMBDA_D_TARGET, M_PHI, b0)
    print(f"\n  α_d(m_φ) needed = 2π / (b₀ · ln(m_φ/Λ_d))")
    print(f"                   = 2π / ({b0:.4f} · {math.log(M_PHI/LAMBDA_D_TARGET):.2f})")
    print(f"                   = {alpha_needed_mphi:.5f}")
    print(f"                   ≈ 1/{1/alpha_needed_mphi:.1f}")
    
    ratio_mphi = alpha_needed_mphi / ALPHA_D
    print(f"\n  Ratio to MCMC:  {ratio_mphi:.2f}×")
    
    # ═════════════════════════════════════════════════════════════════════
    #  Analysis 5: Scan — what if α_D is NOT the SU(2)_d gauge coupling?
    #  What if the SIDM coupling g_D is a Yukawa coupling, and the dark
    #  gauge coupling g_d is separate?
    # ═════════════════════════════════════════════════════════════════════
    print(f"\n{'='*78}")
    print(f"  Analysis 5:  Scan over α_d — from MCMC value to needed value")
    print(f"{'='*78}")
    
    print(f"\n  If α_D(MCMC) ≠ α_d(dark gauge), the gauge coupling is independent.")
    print(f"  What α_d gives what Λ_d?\n")
    
    print(f"  {'α_d':>10s}  {'1/α_d':>8s}  {'Λ_d [GeV]':>12s}  {'Λ_d [meV]':>12s}  {'m_σ [GeV]':>12s}  {'m_σ/H₀':>8s}")
    print("  " + "−"*72)
    
    # Scan at μ = m_χ (characteristic dark sector scale)
    for alpha_test in [0.001, 0.002, 0.003, ALPHA_D, 0.005, 0.010, 
                       0.015, 0.020, 0.025, 0.030, alpha_needed, 
                       0.040, 0.050, 0.060, 0.080, 0.100, 0.118, 0.15, 0.20]:
        Lambda = transmutation_scale(alpha_test, M_CHI, b0)
        Lambda_meV = Lambda * 1e12
        m_sigma = Lambda**2 / F_AXION if Lambda > 0 else 0
        m_over_H0 = m_sigma / H0_PLANCK if m_sigma > 0 and H0_PLANCK > 0 else 0
        
        marker = ""
        if abs(alpha_test - ALPHA_D) < 1e-6:
            marker = " ← MCMC"
        elif abs(alpha_test - alpha_needed) / alpha_needed < 0.01:
            marker = " ← needed for Λ_d=2meV"
        elif abs(alpha_test - 0.118) < 0.001:
            marker = " ← QCD α_s(M_Z)!"
        
        if Lambda_meV > 1e-100:
            print(f"  {alpha_test:10.5f}  {1/alpha_test:8.1f}  {Lambda:12.3e}  {Lambda_meV:12.3e}  {m_sigma:12.3e}  {m_over_H0:8.2f}{marker}")
        else:
            print(f"  {alpha_test:10.5f}  {1/alpha_test:8.1f}  {'≈ 0':>12s}  {'≈ 0':>12s}  {'≈ 0':>12s}  {'≈ 0':>8s}{marker}")
    
    # ═════════════════════════════════════════════════════════════════════
    #  Analysis 6: The QCD analogy
    # ═════════════════════════════════════════════════════════════════════
    print(f"\n{'='*78}")
    print(f"  Analysis 6:  QCD analogy — sanity check")
    print(f"{'='*78}")
    
    # QCD: SU(3), N_f=5 active at M_Z
    b0_qcd = (11.0/3.0) * 3.0 - (4.0/3.0) * 0.5 * 5  # SU(3), Dirac
    # b0_qcd = 11 - 10/3 = 23/3
    alpha_s_MZ = 0.1179
    M_Z = 91.2  # GeV
    Lambda_QCD = transmutation_scale(alpha_s_MZ, M_Z, b0_qcd)
    print(f"\n  QCD: SU(3), N_f=5 Dirac at M_Z")
    print(f"    b₀(QCD) = {b0_qcd:.4f} (= 23/3)")
    print(f"    α_s(M_Z) = {alpha_s_MZ}")
    print(f"    Λ_QCD = M_Z · exp(−2π/(b₀·α_s)) = {Lambda_QCD:.3f} GeV = {Lambda_QCD*1e3:.0f} MeV")
    print(f"    (Experimental: Λ_QCD ≈ 220 MeV — 1-loop estimate is O(1) correct)")
    
    print(f"\n  Lesson: QCD works because α_s(M_Z) ≈ 0.12 is O(0.1).")
    print(f"  Our α_D(MCMC) ≈ 0.003 is O(0.001) — 40× smaller.")
    print(f"  Transmutation is exponentially sensitive to α.")
    
    # ═════════════════════════════════════════════════════════════════════
    #  Summary
    # ═════════════════════════════════════════════════════════════════════
    print(f"\n{'='*78}")
    print(f"  SUMMARY")
    print(f"{'='*78}")
    
    print(f"""
  1. MCMC α_D = {ALPHA_D:.4e} is the SIDM coupling at low energy.
     Forward transmutation gives Λ_d ≈ 10^{math.log10(Lambda_from_mchi):.0f} GeV — essentially zero.
     
  2. To get Λ_d = 2 meV, we need α_d(m_χ) ≈ {alpha_needed:.4f} (1/{1/alpha_needed:.0f}),
     which is {ratio:.1f}× the MCMC value.
     
  3. RG running doesn't help: Λ_d is RG-invariant (by definition).
     Running α from m_φ to M_Pl just redistributes the same information.
     
  4. KEY INSIGHT: α_D (SIDM Yukawa coupling) ≠ α_d (dark gauge coupling).
     
     In our model, χ couples to φ via Yukawa: g_D φ χ̄χ
     But φ is a dark Higgs, and the dark gauge coupling g_d is SEPARATE.
     
     α_D = g_D²/(4π) — Yukawa coupling → SIDM σ_T
     α_d = g_d²/(4π) — gauge coupling → confinement → Λ_d
     
     These are DIFFERENT couplings unless there's a specific relation
     (e.g., g_D = g_d from gauge-Yukawa unification).
     
  5. RESOLUTION: α_d (dark gauge coupling) is a SEPARATE parameter,
     currently hidden in the model. It must satisfy:
     
         α_d(m_χ) ≈ {alpha_needed:.4f}  (≈ 1/{1/alpha_needed:.0f})
     
     to produce Λ_d = 2 meV via transmutation.
     This is NOT fine-tuned — it's analogous to α_s(M_Z) = 0.118 in QCD,
     which gives Λ_QCD ≈ 220 MeV.
     
  6. Open question: Is α_d = {alpha_needed:.4f} consistent with other
     constraints (perturbativity, SIDM, BBN)?
     
     Perturbativity: α_d < 1 ✓ (easily)
     Dark pion mass: m_σ = Λ_d²/f = {LAMBDA_D_TARGET**2/F_AXION:.3e} GeV 
                     → m_σ/H₀ ≈ {LAMBDA_D_TARGET**2/F_AXION/H0_PLANCK:.1f} ✓
""")
    
    print(f"  ╔══════════════════════════════════════════════════════════════╗")
    print(f"  ║  CONCLUSION:  Path 2 transmutation works IF:               ║")
    print(f"  ║                                                            ║")
    print(f"  ║    α_d(dark gauge) ≈ {alpha_needed:.4f} ≈ 1/{1/alpha_needed:.0f}  (separate from α_D)    ║")
    print(f"  ║                                                            ║")
    print(f"  ║  This REDUCES free parameters: Λ_d is no longer free,      ║")
    print(f"  ║  it's derived from α_d via dimensional transmutation.      ║")
    print(f"  ║                                                            ║")
    print(f"  ║  Parameter count: f, Λ_d, θ_i → f, α_d, θ_i              ║")
    print(f"  ║  (same count, but α_d is more fundamental than Λ_d)        ║")
    print(f"  ╚══════════════════════════════════════════════════════════════╝")


if __name__ == "__main__":
    main()
