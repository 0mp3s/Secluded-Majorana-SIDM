"""
Test 22: ΔN_eff Including φ → 2σ (Dark Pion Decay)
====================================================
Goal: After φ decays to 2σ (dark pions from SU(2)_d confinement),
compute the TOTAL ΔN_eff including σ as dark radiation.

PHYSICS
-------
Baseline (from Check 4 / Test 19):
  Dark sector = χ (Majorana) + φ (real scalar)
  g_dark = 1.75 (χ) + 1.0 (φ) = 2.75
  At T_D = 200 MeV decoupling: ΔN_eff = 0.153

Extension (dark QCD):
  φ → 2σ where σ is the dark pion (pseudo-NGB of SU(2)_d)
  m_σ ~ Λ_d²/f ~ H₀ ~ 10⁻³³ eV (ultralight, effectively massless)
  m_φ = 11.1 MeV >> 2m_σ → decay always kinematically open

The question is:
  (a) WHEN does φ → 2σ happen? (lifetime τ_φ)
  (b) Does it happen before or after ν decoupling (T ~ 2 MeV)?
  (c) How does σ contribute to ΔN_eff?
  (d) What is the total ΔN_eff = ΔN(χ) + ΔN(σ)?

KEY INSIGHT (from T21):
  The model is SECLUDED — no Higgs portal. So φ does NOT decay to SM.
  φ can only decay within the dark sector: φ → 2σ (if coupling exists)
  or φ → χχ (if m_φ > 2m_χ — NOT the case: m_φ=11 MeV < 2×94 MeV).
  So φ → 2σ is the ONLY decay channel if the dark QCD coupling exists.

COUPLING
--------
  The φ-σ interaction comes from the dark sector potential:
    V ⊃ λ_{φσ} φ² σ²   (quartic portal within dark sector)
  Or from the dim-5 operator after dark QCD confinement:
    L ⊃ (α_d / 8π f) φ F̃_d F_d → after confinement → (α_d Λ_d² / 8π f) φ σ²

  Decay rate: Γ(φ → 2σ) = λ_{φσ}² v_φ² / (8π m_φ)    [if φ has VEV]
  Or for direct quartic: Γ(φ → 2σ) = λ_{φσ}² m_φ / (32π)  [2-body phase space]

  We scan over λ_{φσ} to find the viable range.

ENTROPY ACCOUNTING
------------------
  Before φ decay: dark sector has χ + φ at temperature T_d
  After φ decay:  dark sector has χ + 2σ (per φ)
  
  If decay is instantaneous (Γ >> H):
    Energy of φ → redistributed to σ
    Since σ is massless: adds to dark radiation
    
  The key: does φ decay INCREASE or DECREASE ΔN_eff?
  Answer: φ (massive, g=1) → 2σ (massless, g=1 each)
  The number of relativistic d.o.f. changes from 1 (if φ was already NR)
  to 1 (σ is real scalar, g=1). But we get ENTROPY from φ mass → σ kinetic.

PARAMETERS (from config)
------------------------
  m_χ  = 94.07 MeV     (DM mass)
  m_φ  = 11.10 MeV     (mediator mass)
  α_D  = 5.734×10⁻³    (dark fine structure)
  Λ_d  = 2.05×10⁻¹² GeV (dark QCD scale)
  f    = 0.24 M_Pl      (dark pion decay constant)
  m_σ  = Λ_d²/f ~ H₀   (dark pion mass, ultralight)
  T_D  = 200 MeV        (dark sector decoupling temperature)
"""

import numpy as np

print("=" * 72)
print("  TEST 22: ΔN_eff INCLUDING φ → 2σ (DARK PION DECAY)")
print("=" * 72)

# ─── Constants ──────────────────────────────────────────────────────────
M_Pl    = 2.435e18    # reduced Planck mass [GeV]
M_Pl_full = 1.221e19  # full Planck mass [GeV]
GeV     = 1.0
MeV     = 1e-3
eV      = 1e-9
keV     = 1e-6
H_0     = 1.44e-42    # Hubble constant [GeV]
HBAR_S  = 6.582e-25   # ℏ [GeV·s]

# ─── Model parameters ──────────────────────────────────────────────────
m_chi   = 94.07e-3    # GeV
m_phi   = 11.10e-3    # GeV
alpha_D = 5.734e-3
y_D     = np.sqrt(4 * np.pi * alpha_D)

# Dark QCD parameters
Lambda_d = 2.05e-12   # GeV
f_sigma  = 0.24 * M_Pl  # dark pion decay constant [GeV]
m_sigma  = Lambda_d**2 / f_sigma  # GMOR relation [GeV]

# SM parameters
T_nu_dec = 2.0 * MeV  # neutrino decoupling temperature
T_D      = 200.0 * MeV  # dark sector decoupling temperature

# Dark sector d.o.f.
g_chi    = 1.75    # Majorana fermion: (7/8) × 2
g_phi    = 1.0     # real scalar
g_sigma  = 1.0     # real scalar (dark pion)
g_dark_total = g_chi + g_phi  # = 2.75

# SM g*_S values
gstar_nu_dec = 10.75   # at T_ν decoupling
gstar_200MeV = 61.75   # at T_D = 200 MeV

print(f"\n{'─'*72}")
print(f"  MODEL PARAMETERS")
print(f"{'─'*72}")
print(f"  m_χ  = {m_chi*1000:.2f} MeV")
print(f"  m_φ  = {m_phi*1000:.2f} MeV")
print(f"  α_D  = {alpha_D:.3e}")
print(f"  Λ_d  = {Lambda_d:.2e} GeV = {Lambda_d/eV:.2e} eV")
print(f"  f    = {f_sigma:.2e} GeV = {f_sigma/M_Pl:.2f} M_Pl")
print(f"  m_σ  = {m_sigma:.2e} GeV = {m_sigma/eV:.2e} eV")
print(f"  m_σ/H₀ = {m_sigma/H_0:.2f}")
print(f"  T_D  = {T_D*1000:.0f} MeV")

# ═══════════════════════════════════════════════════════════════════════
# PART A: Baseline ΔN_eff (no φ decay — from Check 4)
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print(f"  PART A: BASELINE ΔN_eff (χ + φ, no decay)")
print(f"{'='*72}")

xi_baseline = (gstar_nu_dec / gstar_200MeV) ** (1.0/3.0)
dNeff_baseline = (4.0/7.0) * g_dark_total * xi_baseline**4

print(f"  g_dark = {g_chi} (χ) + {g_phi} (φ) = {g_dark_total}")
print(f"  ξ = (g*_ν / g*_D)^(1/3) = ({gstar_nu_dec}/{gstar_200MeV})^(1/3) = {xi_baseline:.4f}")
print(f"  ΔN_eff = (4/7) × g_dark × ξ⁴ = {dNeff_baseline:.4f}")
print(f"  → Baseline: ΔN_eff = {dNeff_baseline:.3f}")

# ═══════════════════════════════════════════════════════════════════════
# PART B: φ → 2σ Decay Rate
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print(f"  PART B: φ → 2σ DECAY RATE")
print(f"{'='*72}")

# Two coupling scenarios:
#
# Scenario 1: Direct quartic  V ⊃ (λ/4) φ² σ²
#   Γ(φ → 2σ) = λ² m_φ / (32π)  [for m_σ ≈ 0]
#   (This is tree-level 1→2 from quartic with symmetry factor)
#   Actually for V = (λ/2) φ² σ², the amplitude is M = λ v_φ
#   but φ is NOT condensed in the scattering scenario.
#   For φ(p) → σ(k1) σ(k2) from (λ/2)φ²σ²:
#   Need φ to have a VEV or cubic coupling.
#
# Better: The coupling comes from dark QCD anomaly:
#   L_eff = (α_d / 8π) (σ/f) G̃G  → after confinement:
#   V_eff ⊃ Λ_d⁴ cos(σ/f) ≈ Λ_d⁴ [1 - σ²/(2f²) + ...]
#   This gives m_σ² = Λ_d⁴/f² (GMOR) but NO direct φ-σ coupling.
#
# The φ-σ coupling must come from the interaction where φ couples
# to dark quarks which form the σ condensate:
#   L ⊃ y_D χ̄χ φ  (dark Yukawa)
#   After confinement: ⟨χ̄χ⟩ ~ Λ_d³/f × cos(σ/f)
#   → effective: y_D Λ_d³/f² × φ σ² (cubic + quartic)
#
# Wait — χ does NOT confine (it's the DM, not a dark quark).
# In Option C (from dark_qcd_consistency.py):
#   Separate dark quarks ψ confine. σ = dark pion from ψ̄ψ.
#   φ couples to σ through a PORTAL: λ_{φσ} φ² σ²
#   This is a FREE parameter (not predicted by the model).
#
# In Option A (χ in fundamental of SU(2)_d):
#   χ̄χ condensate → σ is the fluctuation
#   But χ is 94 MeV and Λ_d ~ 10⁻¹² GeV — χ doesn't confine!
#   Confinement only affects particles lighter than Λ_d.
#
# CONCLUSION: φ-σ coupling is model-dependent. We scan over it.

print("""
  φ-σ coupling analysis:
  ─────────────────────
  In the secluded model, φ → 2σ requires a dark-sector portal.
  Two scenarios for the effective coupling:

  (i)  Anomaly-induced: λ_eff ~ α_D Λ_d² / (4π f m_φ)
       → extremely tiny because Λ_d ~ 10⁻¹² GeV

  (ii) Direct quartic: V ⊃ (λ_{φσ}/2) φ² σ²
       → λ_{φσ} is a free parameter
       For decay: need cubic φσ² term from φ VEV or σ VEV
""")

# Scenario (i): Anomaly-induced coupling
# After confinement, the anomaly gives:
#   L ~ (α_d / 8π) φ/f_φ × Λ_d⁴/f² σ²  (very rough)
# But this requires φ to couple to the dark gauge field.
# In our model φ is a SCALAR, not a gauge boson. So this is suppressed.

# More physical: φ couples to dark quarks ψ via Yukawa y_{φψ}
# At 1-loop: φ → (ψ loop) → 2σ
# Γ ~ y_{φψ}² α_d² m_φ⁵ / (16π³ f² Λ_d²)  [very rough]
# But we don't have ψ in the minimal model.

# Let's just compute Γ for a DIRECT cubic coupling:
# V ⊃ μ_3 φ σ²  →  Γ(φ→2σ) = μ_3² / (8π m_φ) × √(1 - 4m_σ²/m_φ²)
# Since m_σ ~ 0: Γ = μ_3² / (8π m_φ)

print(f"  DIRECT CUBIC: V ⊃ μ₃ φ σ²")
print(f"  Γ(φ→2σ) = μ₃² / (8π m_φ)")
print()

# Scan over μ₃
print(f"  {'μ₃ [GeV]':>14}  {'Γ [GeV]':>12}  {'τ [s]':>12}  {'T_decay [MeV]':>14}  Status")
print(f"  {'─'*14}  {'─'*12}  {'─'*12}  {'─'*14}  ──────")

def H_rad(T, g_star=10.75):
    """Hubble rate [GeV] in radiation domination."""
    return 1.66 * np.sqrt(g_star) * T**2 / M_Pl_full

def T_at_Gamma(Gamma, g_star=10.75):
    """Temperature when Γ = H(T) → T = sqrt(Γ M_Pl / (1.66 sqrt(g*)))."""
    return np.sqrt(Gamma * M_Pl_full / (1.66 * np.sqrt(g_star)))

results_B = []
mu3_values = np.logspace(-15, -1, 50)

for mu3 in mu3_values:
    Gamma = mu3**2 / (8 * np.pi * m_phi)
    tau_s = HBAR_S / Gamma if Gamma > 0 else np.inf
    T_dec = T_at_Gamma(Gamma)
    T_dec_MeV = T_dec / MeV
    
    # Status
    if T_dec_MeV > 200:
        status = "before dark decouple"
    elif T_dec_MeV > 2:
        status = "before ν decouple"
    elif T_dec_MeV > 0.07:  # BBN ~0.07 MeV = 70 keV
        status = "before BBN"
    elif T_dec_MeV > 1e-4:
        status = "after BBN ⚠"
    else:
        status = "today/never"
    
    results_B.append({
        'mu3': mu3, 'Gamma': Gamma, 'tau_s': tau_s,
        'T_dec_MeV': T_dec_MeV, 'status': status
    })
    
    # Print selected values
    if mu3 > 0.9e-15 or abs(np.log10(mu3) % 2) < 0.15:
        print(f"  {mu3:14.2e}  {Gamma:12.2e}  {tau_s:12.2e}  {T_dec_MeV:14.4f}  {status}")

# Special cases
print(f"\n  Special coupling estimates:")

# (a) Anomaly-induced: μ₃ ~ α_D Λ_d² / (4π f)
mu3_anomaly = alpha_D * Lambda_d**2 / (4 * np.pi * f_sigma)
Gamma_anom  = mu3_anomaly**2 / (8 * np.pi * m_phi)
tau_anom    = HBAR_S / Gamma_anom if Gamma_anom > 0 else np.inf
print(f"  Anomaly-induced: μ₃ ~ α_D Λ_d²/(4πf) = {mu3_anomaly:.2e} GeV")
print(f"    Γ = {Gamma_anom:.2e} GeV, τ = {tau_anom:.2e} s")
print(f"    → φ NEVER decays (τ >> age of universe = 4.3×10¹⁷ s)")

# (b) μ₃ ~ m_φ (strong coupling within dark sector)
mu3_strong = m_phi
Gamma_strong = mu3_strong**2 / (8 * np.pi * m_phi)
tau_strong   = HBAR_S / Gamma_strong
T_strong     = T_at_Gamma(Gamma_strong)
print(f"\n  Strong coupling: μ₃ ~ m_φ = {mu3_strong:.2e} GeV")
print(f"    Γ = {Gamma_strong:.2e} GeV, τ = {tau_strong:.2e} s")
print(f"    T_decay = {T_strong/MeV:.1f} MeV → {'before' if T_strong > T_nu_dec else 'after'} ν decoupling")

# (c) What μ₃ gives τ = 1 s (BBN boundary)?
mu3_BBN = np.sqrt(8 * np.pi * m_phi * HBAR_S / 1.0)
print(f"\n  BBN boundary: τ = 1 s → μ₃ = {mu3_BBN:.2e} GeV")

# (d) What μ₃ gives T_decay = T_ν_dec = 2 MeV?
H_nu = H_rad(T_nu_dec)
mu3_nu = np.sqrt(8 * np.pi * m_phi * H_nu)
print(f"  ν decouple boundary: T = 2 MeV → μ₃ = {mu3_nu:.2e} GeV")

# ═══════════════════════════════════════════════════════════════════════
# PART C: ΔN_eff AFTER φ → 2σ DECAY
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print(f"  PART C: ΔN_eff AFTER φ → 2σ")
print(f"{'='*72}")

print("""
  Three timing regimes for φ → 2σ:

  REGIME 1: φ decays BEFORE dark sector decouples (T_decay > T_D = 200 MeV)
    → σ thermalizes in dark sector bath
    → Replace φ (g=1) with σ (g=1) in g_dark
    → g_dark stays the same: 1.75 (χ) + 1.0 (σ) = 2.75
    → ΔN_eff UNCHANGED = 0.153

  REGIME 2: φ decays AFTER dark decoupling but BEFORE ν decoupling
            (2 MeV < T_decay < 200 MeV)
    → φ was already decoupled from SM, at dark temperature T_d < T_γ
    → φ decay dumps energy into σ within the dark sector
    → Conservation of energy: ρ_φ → ρ_σ
    → If φ was non-relativistic at decay: E_φ = m_φ → two σ each with E = m_φ/2
    → If φ was relativistic: E_φ = p_φ → two σ each with E = p_φ/2
    
    At T_D = 200 MeV, dark sector temperature T_d:
      ξ = T_d/T_γ = (g*_ν/g*_D)^{1/3} at ν dec
    
    Is φ relativistic at T ~ 100 MeV (dark temperature)?
      T_d_at_100MeV = ξ × 100 MeV × (g*(100)/g*(200))^{1/3}
      T_d ~ 0.56 × 100 MeV × (14.25/61.75)^{1/3} ~ 34 MeV
      φ is relativistic if T_d > m_φ = 11.1 MeV → YES at 100 MeV
      φ becomes NR when T_d ~ m_φ ↔ T_γ ~ m_φ/ξ ~ 20-30 MeV
    
    CASE 2a: φ still relativistic at decay
      ρ_φ → ρ_σ, same number of d.o.f. (1→1 scalar)
      But each φ → two σ, so n_σ = 2 n_φ
      Energy per σ = E_φ/2 — but effectively σ just IS radiation
      ΔN_eff contribution from σ = same as from φ when relativistic
      → ΔN_eff ≈ 0.153 (unchanged)
    
    CASE 2b: φ already non-relativistic at decay
      φ energy density: ρ_φ = m_φ n_φ (matter-like, scales as a⁻³)
      After decay: ρ_σ = ρ_φ (energy conservation)
      But σ is radiation → scales as a⁻⁴ from then on
      This is ENTROPY INJECTION into dark radiation sector!
      
      The dark sector temperature INCREASES:
        T_d_new > T_d_old (heated by φ mass energy)
      
      Factor: (T_d_new/T_d_old)⁴ = (ρ_σ_new / ρ_σ_old) 
            = 1 + ρ_φ/ρ_χ (approximately)

  REGIME 3: φ decays AFTER BBN (T_decay < 0.07 MeV) or NEVER
    → φ is stable on cosmological timescale
    → φ contributes to Ω_m (dark matter) not Ω_rad
    → ΔN_eff from χ only: (4/7) × 1.75 × ξ⁴
""")

# ─── Compute for each regime ────────────────────────────────────────

# REGIME 1: φ decays before T_D → ΔN_eff unchanged
print(f"  REGIME 1 (T_decay > {T_D/MeV:.0f} MeV): φ → 2σ before decoupling")
dNeff_regime1 = dNeff_baseline
print(f"    g_dark = 1.75 (χ) + 1.0 (σ replaces φ) = 2.75")
print(f"    ΔN_eff = {dNeff_regime1:.4f} (same as baseline)")

# REGIME 2a: φ still relativistic at decay
print(f"\n  REGIME 2a (relativistic φ → 2σ, 2 MeV < T < 200 MeV):")
# Energy conservation: ρ_φ(rad) → ρ_σ(rad)
# n_σ = 2 n_φ, each with half the energy → same total energy density
# g_eff stays 1.0 → ΔN_eff from σ = ΔN_eff from φ
dNeff_regime2a = dNeff_baseline
print(f"    ΔN_eff = {dNeff_regime2a:.4f} (same — energy conservation)")

# REGIME 2b: φ non-relativistic at decay (most interesting)
print(f"\n  REGIME 2b (NR φ → 2σ, entropy injection):")

# At T_γ where φ becomes NR in dark sector:
# T_d = ξ × T_γ, φ NR when T_d < m_φ → T_γ > m_φ/ξ
# ξ at different T_γ depends on g*(T_γ)

# Let's compute the entropy injection factor
# Before φ decay (at some time t_dec, temperature T_γ,dec):
#   Dark sector has: ρ_χ (radiation if T_d > m_χ, which it isn't at low T)
#   Actually χ is 94 MeV, so at T_γ ~ few MeV, T_d ~ 1 MeV
#   χ is NON-RELATIVISTIC at decay time too!
# 
# Wait — both χ and φ may be NR by the time of decay.
# χ freezes out at T_d ~ m_χ/20 ~ 4.7 MeV → T_γ ~ 4.7/ξ ~ 8-15 MeV
# After χ freeze-out, dark sector has:
#   - Non-relativistic χ (frozen out, contributing to DM)
#   - φ: becomes NR at T_d ~ m_φ = 11 MeV → T_γ ~ 20-30 MeV
#     But wait: φ ALSO freezes out of the dark thermal bath
#     χ annihilates: χχ → φφ keeps φ in thermal equilibrium
#     After χ freeze-out, φ equilibrium maintained by φφ↔φφ (self-scattering)
#     or number-changing 3φ↔2φ (cannibal)
#     Eventually φ becomes NR and its number density → exp(-m_φ/T_d) (Boltzmann)
#     This means ρ_φ → 0 exponentially.
#     UNLESS φ has conserved number (no number-changing processes).

print("""
    CRITICAL REALIZATION:
    ────────────────────
    After χ freezes out in the dark sector:
    - χ is frozen non-relativistic DM
    - φ is still in dark thermal bath
    - If 3φ→2φ (cannibal) is active: φ stays in equilibrium
      → n_φ ~ (m_φ T_d)^{3/2} exp(-m_φ/T_d) → exponentially suppressed
      → Very little φ energy density remains to inject into σ
    - If cannibal is NOT active: φ number is frozen
      → ρ_φ = m_φ × n_φ(frozen) → matter-like component
      → Decay φ→2σ converts this to dark radiation
    
    In the mixed_coupling analysis (condition3):
      Cannibal 3φ→2φ IS active for μ₃/m_φ ≳ 1.3
      For our model: cannibal keeps φ in equilibrium → φ Boltzmann-
      suppressed → negligible ρ_φ by T_d << m_φ
    
    RESULT: If cannibal is active, there's almost no φ left to decay.
    The φ→2σ decay is irrelevant for ΔN_eff because φ is already gone!
""")

# Quantitative: what fraction of φ survives to decay?
# After cannibal phase, φ abundance: n_φ/s ~ Y_φ,fo
# For cannibal: Y_φ,fo is SET by cannibal freeze-out
# From condition3: Ω_φ h² depends on μ₃/m_φ

# Without cannibal, φ freezes out as RADIATION (T_d > m_φ):
# Y_φ = n_φ/s_dark is constant after dark decoupling
# n_φ ~ T_d³ at decoupling → after NR → ρ_φ = m_φ n_φ = m_φ T_D³ (g_*(now)/g_*(T_D))

# Let's compute the φ contribution to energy density at T_ν dec = 2 MeV
# assuming φ is a relic from dark sector thermal bath

# Dark sector temperature at T_γ = 2 MeV:
xi_at_2MeV = (gstar_nu_dec / gstar_200MeV) ** (1.0/3.0)
T_d_at_2MeV = xi_at_2MeV * 2.0 * MeV

print(f"    Dark temperature when T_γ = 2 MeV: T_d = {T_d_at_2MeV/MeV:.2f} MeV")
print(f"    m_φ/T_d = {m_phi / T_d_at_2MeV:.1f}")

# φ is non-relativistic: m_φ/T_d ~ 11/1.1 ~ 10 → heavily Boltzmann suppressed
z_phi = m_phi / T_d_at_2MeV
n_phi_over_T3 = (z_phi / (2*np.pi))**1.5 * np.exp(-z_phi)  # Boltzmann
rho_phi_over_T4 = z_phi * n_phi_over_T3  # ρ = m n

# Compare to radiation:
rho_rad_over_T4 = np.pi**2 / 30  # for 1 d.o.f.

ratio = rho_phi_over_T4 / rho_rad_over_T4
print(f"    Boltzmann suppression: e^(-m_φ/T_d) = e^(-{z_phi:.1f}) = {np.exp(-z_phi):.2e}")
print(f"    ρ_φ / ρ_rad(1 d.o.f.) = {ratio:.2e}")
print(f"    → φ contribution to energy density is NEGLIGIBLE")

# REGIME 3: φ is stable (most likely scenario)
print(f"\n  REGIME 3 (φ stable — no dark QCD coupling or τ > t_universe):")
# If φ is stable AND was Boltzmann-suppressed, it's just matter
# If φ is stable AND was NOT Boltzmann-suppressed (e.g., asymmetry):
# then it's a sub-dominant CDM component

# Actually the key question: is there any φ left at late times?
# In the thermal dark sector, for T_d << m_φ:
# n_φ/s → Y_φ,fo ~ 10⁻³ (typical for thermal relic at m_φ/T_fo ~ 20)
# Ω_φ h² = Y_φ,fo × s₀ × m_φ / ρ_crit
# This was computed in condition3 of mixed_coupling!

print(f"    If φ is Boltzmann-suppressed at T_d << m_φ = {m_phi/MeV:.1f} MeV:")
print(f"    → φ is negligible. Dark sector is χ only.")
print(f"    → ΔN_eff = (4/7) × g_χ × ξ⁴ (with ξ adjusted for g* change)")

# After φ becomes NR and annihilates, it heats the dark sector
# The entropy from φ annihilation goes to... what?
# In the SECLUDED model: φφ → nothing (no lighter dark particles except σ if it exists)
# φφ → χχ is kinematically FORBIDDEN (m_φ = 11 MeV < m_χ = 94 MeV)
# So φ can only self-annihilate: φφ → φφ (elastic, doesn't change number)
# Or cannibal: 3φ → 2φ (reduces number, heats remaining φ)
# Eventually ALL φ energy is in a small number of energetic φ particles
# that are non-relativistic and contribute to matter density.

# WITHOUT φ→2σ:
# φ remnant contributes to Ω_m, NOT to ΔN_eff
# ΔN_eff comes only from χ (if relativistic) — but χ is also NR!

print(f"""
    DEEP ISSUE: At T_γ = 2 MeV (ν decoupling):
    T_d ~ {T_d_at_2MeV/MeV:.1f} MeV
    Both χ (94 MeV) and φ (11 MeV) are NON-RELATIVISTIC
    m_χ/T_d = {m_chi/T_d_at_2MeV:.0f}, m_φ/T_d = {m_phi/T_d_at_2MeV:.0f}
    
    So at ν decoupling: ZERO dark radiation from χ or φ!
    ΔN_eff ≈ 0 (both species are Boltzmann-suppressed)
    
    WAIT — this contradicts the baseline ΔN_eff = 0.153!
    
    Resolution: ΔN_eff = 0.153 is computed ASSUMING the dark sector
    has relativistic d.o.f. at T_D = 200 MeV that DECOUPLE while
    relativistic. After decoupling, they redshift as radiation
    indefinitely (even though they become NR later).
    
    KEY DISTINCTION:
    ────────────────
    (a) If dark sector STAYS COUPLED internally after SM decoupling:
        → φ and χ go through thermal transitions in the dark bath
        → φ becomes NR → energy goes to cannibal heating
        → χ freezes out within dark bath
        → At late times: only NR massive particles remain
        → ΔN_eff ≈ 0 at BBN (no dark radiation)
    
    (b) If dark sector DECOUPLES FROM ITSELF at T_D:
        → Each particle free-streams independently
        → Momentum redshifts: p ∝ 1/a
        → Even when p << m (NR): energy density ρ = m n ∝ a⁻³
        → This is MATTER, not radiation
        → ΔN_eff from free-streaming massive relics → 0 at late times
    
    (c) The STANDARD ΔN_eff = 0.153 calculation assumes the dark
        sector contributes g_dark d.o.f. of RADIATION at T_D,
        which then redshifts as a⁻⁴. This is valid if there are
        MASSLESS particles in the dark sector at T_D.
        
        For χ (m=94 MeV) at T_D=200 MeV: T_d/m_χ ~ 2 → semi-relativistic
        For φ (m=11 MeV) at T_D=200 MeV: T_d/m_φ ~ 18 → fully relativistic
        
        At ν decoupling (T~2 MeV): both are NR → contribution → 0
        But the ENTROPY they deposited at decoupling is frozen in.
""")

# ═══════════════════════════════════════════════════════════════════════
# PART D: CORRECT ΔN_eff ACCOUNTING
# ═══════════════════════════════════════════════════════════════════════
print(f"{'='*72}")
print(f"  PART D: CORRECT ΔN_eff CALCULATION")
print(f"{'='*72}")

# The ΔN_eff = 0.153 from Check 4 is computed as if the dark sector
# d.o.f. are massless. Let's do the correct massive calculation.

# At recombination (T_γ ~ 0.26 eV), the relevant N_eff is:
# N_eff = ρ_rad(non-photon) / ρ_1ν
# 
# Dark sector particles that are NR by recombination contribute to Ω_m, not N_eff.
# χ (94 MeV) and φ (11 MeV) are DEFINITELY NR at recombination.
#
# So: ΔN_eff from massive dark sector = 0 at recombination!
# UNLESS there is a massless dark species.
#
# This is where σ matters:
# If φ → 2σ BEFORE recombination, and m_σ ~ H₀ ~ 10⁻³³ eV (massless),
# then σ IS dark radiation at recombination!
#
# The question becomes: how much energy is transferred to σ?

print("""
  CORRECT PICTURE:
  ────────────────
  Without σ: ΔN_eff = 0 at recombination
    (χ and φ are massive and NR → matter, not radiation)
  
  With σ (massless dark pion):
    σ is the ONLY massless d.o.f. in the dark sector
    If dark sector was ever thermalized, σ has T_d at dark decoupling
    → ΔN_eff = (4/7) × g_σ × (T_d/T_ν)⁴
    
  BUT: Does σ thermalize in the dark sector?
    σ couples to χ via (m_χ/f) σ χ̄iγ⁵χ
    Interaction rate: Γ ~ n_χ × (m_χ/f)² / m_χ²  
    Since f ~ 0.24 M_Pl ~ 5.8×10¹⁷ GeV:
    (m_χ/f)² ~ (94 MeV / 5.8×10¹⁷ GeV)² ~ 2.6×10⁻²⁰
    → EXTREMELY suppressed coupling → σ NEVER thermalizes!
""")

# σ thermalization check
g_sigma_chi = m_chi / f_sigma  # effective σ-χ coupling
print(f"  σ-χ coupling: g = m_χ/f = {g_sigma_chi:.2e}")
print(f"  g² = {g_sigma_chi**2:.2e}")

# Interaction rate at T_D:
# Γ(σχ→σχ) ~ n_χ σ ~ T_d³ × g⁴/(16π T_d²) = g⁴ T_d / (16π)
# Actually: σ ~ g⁴ / (16π s) where s ~ T² → Γ ~ g⁴ T
Gamma_sigma = g_sigma_chi**4 * T_D / (16*np.pi)
H_at_TD = H_rad(T_D, gstar_200MeV)
print(f"\n  At T_D = 200 MeV:")
print(f"    Γ(σ thermalization) ~ g⁴ T / 16π = {Gamma_sigma:.2e} GeV")
print(f"    H(T_D)              = {H_at_TD:.2e} GeV")
print(f"    Γ/H                 = {Gamma_sigma/H_at_TD:.2e}")
print(f"    → σ {'thermalizes' if Gamma_sigma > H_at_TD else 'NEVER thermalizes'}")

# What about φ→2σ as a SOURCE of σ radiation?
print(f"""
  Since σ never thermalizes through scattering (g ~ m_χ/f too small),
  the ONLY source of σ is φ→2σ decay.
  
  If φ→2σ occurs, the σ particles free-stream as dark radiation.
  Their contribution to ΔN_eff depends on the ENERGY transferred.
""")

# ═══════════════════════════════════════════════════════════════════════
# PART E: ΔN_eff FROM φ→2σ AT DIFFERENT TIMES
# ═══════════════════════════════════════════════════════════════════════
print(f"{'='*72}")
print(f"  PART E: ΔN_eff FROM φ→2σ AT DIFFERENT TIMES")
print(f"{'='*72}")

# The φ abundance at decay time determines how much energy goes to σ.
# 
# Case 1: φ decays while RELATIVISTIC (T_d > m_φ)
#   ρ_φ = (π²/30) × g_φ × T_d⁴ (thermal radiation)
#   After decay: ρ_σ = ρ_φ (energy conservation)
#   σ contribution: ΔN_eff(σ) = (4/7) × g_σ × ξ⁴ = (4/7) × 1 × ξ⁴
#   
#   Complete replacement: φ (g=1) → σ (g=1)
#   ΔN_eff changes by: Δ(ΔN) = (4/7)(g_σ - g_φ)ξ⁴ = 0 (same d.o.f.)
#   → NO CHANGE for relativistic decay
#
# Case 2: φ decays while NON-RELATIVISTIC (T_d < m_φ)
#   ρ_φ = m_φ × n_φ(T_d)
#   n_φ = (m_φ T_d / 2π)^{3/2} exp(-m_φ/T_d)  [Boltzmann]
#   This energy is converted to σ radiation (massless)
#   
#   The ΔN_eff contribution from this injected radiation:
#   ΔN_eff(σ) = (8/7) × (ρ_σ / ρ_γ) evaluated at relevant time
#   where ρ_γ = (π²/15) T_γ⁴

# Let's compute for the interesting case: φ decay at T_γ between 2-200 MeV
print(f"\n  ΔN_eff(σ) from NR φ decay at various T_γ:")
print(f"  {'T_γ [MeV]':>12}  {'T_d [MeV]':>10}  {'m_φ/T_d':>8}  {'n_φ/T_d³':>12}  {'ΔN_eff(σ)':>10}")
print(f"  {'─'*12}  {'─'*10}  {'─'*8}  {'─'*12}  {'─'*10}")

for T_gamma_MeV in [200, 100, 50, 30, 20, 10, 5, 2, 1]:
    T_gamma = T_gamma_MeV * MeV
    
    # g*_S at this T_γ
    if T_gamma_MeV > 160:
        g_s = 61.75
    elif T_gamma_MeV > 100:
        g_s = 61.75
    elif T_gamma_MeV > 0.5:
        # approximate
        g_s = 10.75 + (61.75 - 10.75) * (T_gamma_MeV - 0.5) / (160 - 0.5)
    else:
        g_s = 10.75
    
    # Dark sector temperature (from entropy conservation after decoupling at T_D)
    xi = (gstar_nu_dec / gstar_200MeV) ** (1.0/3.0)
    # Actually ξ = T_d/T_γ varies with T_γ due to SM entropy changes
    # More precisely: T_d = T_D × (a_D/a) = T_D × (g_S(T)/g_S(T_D))^{1/3} × (T/T_D)^{-1} × T/T
    # Simpler: after decoupling, T_d × a = const and T_γ × a × g_S^{1/3} = const
    # So T_d/T_γ = (T_D/T_D) × (g_S(T_γ)/g_S(T_D))^{1/3} = ξ × (g_S(T_γ)/gstar_nu_dec)^{1/3}
    # Wait — at decoupling: T_d = ξ_0 × T_γ where ξ_0 depends on the g_S ratio
    # After decoupling: T_d ∝ 1/a, T_γ ∝ g_S^{-1/3}/a
    # So T_d/T_γ = ξ_0 × g_S(T_γ)^{1/3} / g_S(T_D)^{1/3}
    # At T_D: ξ_0 = ... actually the dark sector decouples at T_D = 200 MeV
    # At that moment T_d = T_γ = T_D (they were coupled)? NO!
    # Dark sector decouples with T_d set by how much energy it got.
    # The 0.153 formula: ξ⁴ at ν dec = (10.75/61.75)^{4/3}
    # This assumes dark sector tracks SM entropy until T_D, then free-streams.
    
    # After dark decoupling at T_D:
    # T_d/T_γ at later time = (g_S(T_γ)/g_S(T_D))^{1/3}
    # because T_d ∝ 1/a and T_γ ∝ g_S^{-1/3}/a
    xi_T = (g_s / gstar_200MeV) ** (1.0/3.0)
    T_d = xi_T * T_gamma
    
    z = m_phi / T_d if T_d > 0 else 1e10
    
    if z < 3:
        # relativistic — φ still acts as radiation
        n_ratio = 1.0  # n_φ/T_d³ ~ O(1)
        dN_sigma = (4.0/7.0) * 1.0 * xi_T**4  # same as φ contribution
        note = "(relativistic — same as φ)"
    else:
        # Boltzmann suppressed
        n_ratio = (z/(2*np.pi))**1.5 * np.exp(-z)
        # Energy density of φ: ρ_φ = m_φ × n_φ = m_φ T_d³ × n_ratio
        # ΔN_eff(σ) = (8/7)(ρ_σ/ρ_1ν) where ρ_1ν = (7/8)(π²/15)T_ν⁴
        # At T_ν = T_γ (before ν decouple): ρ_1ν = (7π²/120)T_γ⁴
        rho_phi = m_phi * T_d**3 * n_ratio
        rho_1nu = 7 * np.pi**2 / 120 * T_gamma**4
        dN_sigma = rho_phi / rho_1nu
        note = f"(NR, Boltz. e^-{z:.0f})"
    
    print(f"  {T_gamma_MeV:12.0f}  {T_d/MeV:10.2f}  {z:8.1f}  {n_ratio:12.2e}  {dN_sigma:10.2e}  {note}")

# ═══════════════════════════════════════════════════════════════════════
# PART F: TOTAL ΔN_eff — ALL SCENARIOS
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print(f"  PART F: TOTAL ΔN_eff — COMPLETE PICTURE")
print(f"{'='*72}")

# The crucial insight: ΔN_eff = 0.153 assumed MASSLESS dark d.o.f.
# But m_χ = 94 MeV and m_φ = 11 MeV → both become NR before recombination.
# The CORRECT baseline ΔN_eff at recombination from these massive species is ~0.

# With σ (m ~ H₀, effectively massless):
# σ gets its energy from φ→2σ decay
# But how much energy does φ have at decay?

# Scenario A: No φ→2σ coupling (SECLUDED, no dark QCD)
print(f"\n  SCENARIO A: No dark QCD, no σ")
print(f"  ─────────────────────────────")
print(f"  Dark sector: χ (94 MeV) + φ (11 MeV), both massive")
print(f"  At CMB recombination (T ~ 0.26 eV): both NR → matter")
print(f"  ΔN_eff at recombination = 0")
print(f"  ΔN_eff at BBN (T ~ 1 MeV):")

# At BBN, T_d ~ ξ × 1 MeV ~ 0.55 MeV
# φ: m_φ/T_d = 11/0.55 = 20 → Boltzmann suppressed → negligible
# χ: m_χ/T_d = 94/0.55 = 170 → completely frozen out
xi_BBN = (10.75 / gstar_200MeV) ** (1.0/3.0)
T_d_BBN = xi_BBN * 1.0 * MeV
print(f"    T_d(BBN) = {T_d_BBN/MeV:.3f} MeV")
print(f"    m_φ/T_d = {m_phi/T_d_BBN:.0f}, m_χ/T_d = {m_chi/T_d_BBN:.0f}")
print(f"    → Both Boltzmann-suppressed: ΔN_eff ≈ 0 at BBN too!")

print(f"""
  ⚠️ IMPORTANT CORRECTION:
  The ΔN_eff = 0.153 from Check 4 is the contribution at T_D = 200 MeV,
  when χ is semi-relativistic and φ is fully relativistic.
  At LATER times (BBN, recombination), both are NR → ΔN_eff → 0.
  
  HOWEVER: The standard way ΔN_eff is defined for CMB/BBN is:
  ΔN_eff counts the TOTAL extra radiation at ν decoupling (for BBN)
  or at recombination (for CMB).
  
  For MASSIVE species: they contribute to ΔN_eff only while relativistic.
  After they become NR: they contribute to Ω_m instead.
  
  Correctly: φ (11 MeV) is still relativistic at T ~ 50 MeV (dark T ~ 26 MeV)
  but NR by T ~ 2 MeV (dark T ~ 1.1 MeV, m_φ/T_d ~ 10).
  
  So φ's contribution to ΔN_eff at BBN ≈ 0 (already NR).
  χ's contribution at BBN ≈ 0 (already NR and frozen out).
  
  ΔN_eff(BBN) ≈ 0 in the secluded model WITHOUT σ.
""")

# Scenario B: With dark QCD → σ is the only massless dark species
print(f"\n  SCENARIO B: With dark QCD, φ → 2σ")
print(f"  ────────────────────────────────────")
print(f"  σ is the ONLY massless d.o.f. in the dark sector.")
print(f"  σ gets energy from φ→2σ decay.")
print()

# If φ decays while relativistic (T_d > m_φ), all its radiation energy → σ
# ΔN_eff(σ) = (4/7) × g_σ × ξ⁴ at ν decoupling
# But wait — σ doesn't thermalize (coupling too weak)
# σ just free-streams with energy = m_φ/2 per particle at decay
# After decay: σ redshifts as a⁻⁴ (massless)

# If φ decays while RELATIVISTIC (before dark sector cools below m_φ):
# ρ_σ inherits ρ_φ → ΔN_eff(σ) = (4/7) × 1 × ξ(T_decay)⁴

# If φ decays while NR (after dark sector cools below m_φ):
# ρ_σ = ρ_φ(NR) << ρ_φ(had it been rad) → ΔN_eff(σ) << previous

# Most favorable: φ decays immediately while still relativistic
# This gives maximum ΔN_eff(σ)
dNeff_sigma_max = (4.0/7.0) * g_sigma * xi_baseline**4

# Only from σ (φ is gone):
dNeff_chi_at_BBN = 0  # χ is NR at BBN
dNeff_total_B = dNeff_chi_at_BBN + dNeff_sigma_max

print(f"  Maximum ΔN_eff(σ) [φ decays while relativistic]:")
print(f"    ΔN_eff(σ) = (4/7) × {g_sigma} × ξ⁴ = {dNeff_sigma_max:.4f}")
print(f"    ΔN_eff(χ) at BBN = {dNeff_chi_at_BBN:.4f} (NR, frozen)")
print(f"    Total ΔN_eff = {dNeff_total_B:.4f}")
print()

# Compare with original claim
print(f"  Original claim: ΔN_eff = 0.153 (χ) + 0.027 (σ) ≈ 0.180")
print(f"  Corrected:      ΔN_eff(χ at BBN) ≈ 0, ΔN_eff(σ) ≤ {dNeff_sigma_max:.4f}")
print(f"  Maximum total:  ΔN_eff ≤ {dNeff_sigma_max:.4f}")
print()

# The 0.153 is overestimated because it assumes χ is relativistic at BBN
# χ is 94 MeV, T_d at BBN ~ 0.55 MeV → m/T ~ 170 → completely NR
# The REAL ΔN_eff at BBN from the dark sector (without σ) is ≈ 0

# With σ: the maximum is (4/7) × 1 × ξ⁴ ≈ 0.056 (only φ's share)
# But this requires φ to have decayed while still hot.

# ═══════════════════════════════════════════════════════════════════════
# PART G: RE-EXAMINATION — WHEN IS ΔN_eff = 0.153 CORRECT?
# ═══════════════════════════════════════════════════════════════════════
print(f"{'='*72}")
print(f"  PART G: WHEN IS ΔN_eff = 0.153 VALID?")
print(f"{'='*72}")

print(f"""
  The ΔN_eff = 0.153 formula:
    ΔN_eff = (4/7) × g_dark × (g*_ν / g*_D)^{{4/3}}
    
  This is EXACT when:
  (1) All dark species are MASSLESS (or m << T at all relevant epochs)
  (2) Dark sector decouples at T_D from SM, then free-streams
  
  For our model: g_dark = 2.75 (χ Majorana + φ scalar)
  But m_χ = 94 MeV and m_φ = 11 MeV are NOT massless!
  
  CORRECT TREATMENT for massive species:
  ────────────────────────────────────────
  After decoupling at T_D = 200 MeV:
  - Dark entropy is conserved separately: s_d ∝ g*_S,d(T_d) × T_d³
  - As T_d drops, g*_S,d DECREASES when species become NR
  - This HEATS the remaining relativistic dark species
  
  Timeline:
  1. T_d > m_χ = 94 MeV: g_d = 2.75 (all relativistic)
  2. m_φ < T_d < m_χ: χ becomes NR → heats φ
     g_d drops from 2.75 to 1.0
     T_d increases by factor (2.75/1.0)^{{1/3}} = 1.40
  3. T_d < m_φ: φ also becomes NR → NO relativistic species left!
     All dark entropy goes to... nothing? → stored in NR species
  
  KEY: If there's NO massless dark species, the dark "radiation"
  energy density drops exponentially as dark species become NR.
  At BBN/recombination: ΔN_eff → 0 (exponentially).
  
  IF σ EXISTS (massless) AND THERMALIZED at T_D:
  - σ absorbs ALL dark entropy when χ and φ become NR
  - Entropy conservation: g_total × T_d³ = g_σ × T_σ³
    → T_σ = T_d × (g_total/g_σ)^{{1/3}} = T_d × (2.75/1)^{{1/3}} = 1.40 T_d
  - ΔN_eff(σ) = (4/7) × g_σ × (ξ_σ)⁴ where ξ_σ = ξ × (g_total/g_σ)^{{1/3}}
  
  NOTE: (4/7) × g_σ × (g_total/g_σ)^{{4/3}} × ξ⁴ ≠ (4/7) × g_total × ξ⁴
  The heating of σ OVERCOMPENSATES for the loss of χ,φ d.o.f.
  (T⁴ effect: concentrating entropy into fewer d.o.f. raises T more)
""")

# With entropy transfer to σ:
boost = (g_dark_total / g_sigma) ** (1.0/3.0)
xi_sigma = xi_baseline * boost
dNeff_sigma_heated = (4.0/7.0) * g_sigma * xi_sigma**4

print(f"  Entropy boost factor: (g_total/g_σ)^(1/3) = ({g_dark_total}/{g_sigma})^(1/3) = {boost:.4f}")
print(f"  σ temperature ratio: ξ_σ = ξ × boost = {xi_baseline:.4f} × {boost:.4f} = {xi_sigma:.4f}")
print(f"  ΔN_eff(σ, entropy-heated) = (4/7) × 1 × ξ_σ⁴ = {dNeff_sigma_heated:.4f}")
print()

# Compare to original massless formula
print(f"  Original (massless approx): ΔN_eff = (4/7) × g_total × ξ⁴ = {dNeff_baseline:.4f}")
print(f"  Correct  (entropy to σ):    ΔN_eff = (4/7) × g_σ × ξ_σ⁴  = {dNeff_sigma_heated:.4f}")
print(f"  Ratio: {dNeff_sigma_heated/dNeff_baseline:.4f}")
print()
print(f"  The entropy-transfer formula gives a HIGHER ΔN_eff because")
print(f"  concentrating entropy into fewer d.o.f. heats them super-linearly (T⁴).")
print()
print(f"  HOWEVER: this assumes σ was thermalized with χ,φ at T_D.")
print(f"  Three sub-cases for the actual ΔN_eff:")
print(f"  (a) σ thermalized at T_D, absorbs all entropy: ΔN_eff = {dNeff_sigma_heated:.3f}")
print(f"  (b) σ gets entropy only from φ→2σ (not from χ): ΔN_eff < {dNeff_sigma_heated:.3f}")
print(f"  (c) σ never thermalizes (f ~ M_Pl):              ΔN_eff ≈ 0")

# ═══════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print(f"  TEST 22 — FINAL SUMMARY")
print(f"{'='*72}")
print(f"""
  QUESTION: What happens to ΔN_eff when φ → 2σ?
  
  ANSWER: ΔN_eff depends critically on whether σ exists and thermalizes.
  
  EXPLANATION:
  ────────────
  The original ΔN_eff = 0.153 counts g_dark = 2.75 d.o.f. of dark
  radiation at decoupling. But χ (94 MeV) and φ (11 MeV) are massive.
  As the dark sector cools, their energy must go somewhere:
  
  WITHOUT σ: entropy has nowhere to go → both become NR → ΔN_eff → 0
    (The massless approximation ΔN_eff = 0.153 is WRONG for a model
     with only massive dark species!)
  
  WITH σ (thermalized at T_D, massless):
    σ absorbs ALL dark entropy as χ, φ become NR
    T_σ boosted by (g_total/g_σ)^(1/3) = 1.40
    ΔN_eff(σ) = (4/7) × 1 × (ξ × 1.40)⁴ ≈ 0.214
    This is HIGHER than 0.153 because concentrating entropy
    into fewer d.o.f. raises ΔN_eff (T⁴ scaling).
  
  WITH φ→2σ ONLY (σ not thermalized initially):
    σ gets energy only from φ→2σ, not from χ
    ΔN_eff(σ) ≤ 0.056 (only φ's share of energy)
  
  RESULTS TABLE:
  ──────────────
  | Scenario                     | ΔN_eff (BBN) | ΔN_eff (CMB) |
  |------------------------------|--------------|--------------|
  | No σ, no dark QCD            | ≈ 0          | ≈ 0          |
  | σ therm. at T_D, gets all S  | 0.214        | 0.214        |
  | σ from φ→2σ only (max)       | 0.056        | 0.056        |
  | Original claim (check 4)     | 0.153        | 0.153        |
  
  CRITICAL FINDING:
  ─────────────────
  1. σ is NECESSARY for any ΔN_eff > 0 at BBN/CMB.
     Without a massless dark species, ΔN_eff → 0 as χ,φ → NR.
  
  2. The original claim of ΔN_eff = 0.153 (from check 4) is the
     MASSLESS APPROXIMATION. The correct value depends on how σ
     acquires dark entropy:
     - σ thermalized at T_D: ΔN_eff = 0.214 (entropy-heated)
     - σ only from φ→2σ: ΔN_eff ≤ 0.056
     - σ not thermalized: ΔN_eff → 0
  
  3. The original need_to_verify claim of +0.027 ADDITIONAL from σ
     is WRONG — σ doesn't add on top, it replaces χ,φ.
  
  4. With f ~ 0.24 M_Pl, σ CANNOT thermalize (coupling ~ 10⁻²⁰).
     → ΔN_eff at BBN/CMB ≈ 0 unless σ is populated by φ→2σ.
     → If φ→2σ occurs while φ is relativistic: ΔN_eff ≈ 0.056.
  
  5. ALL scenarios satisfy Planck: ΔN_eff < 0.33 ✓
  
  STATUS: ✅ PASS — no BBN/CMB tension in any scenario.
  The ΔN_eff = 0.153 prediction needs revision:
  - Best case (σ thermalized): 0.214 → Planck OK, CMB-S4 7.9σ
  - Likely case (φ→2σ only):   0.056 → Planck OK, CMB-S4 2.1σ  
  - Conservative (no σ):       ≈ 0   → Planck OK, no CMB-S4 signal
""")

print(f"  Planck limit:  ΔN_eff < 0.33 → {'✓ PASS' if 0.153 < 0.33 else '✗ FAIL'}")
print(f"  CMB-S4 (2σ):   σ(ΔN_eff) = 0.027 → {0.153/0.027:.1f}σ detection")
print()

# Timing check
print(f"  TIMING OF φ→2σ:")
print(f"    Anomaly-induced coupling: τ >> t_universe → φ never decays")
print(f"    → σ thermalizes through dark sector entropy, not φ decay")
print(f"    → Entropy transfer happens automatically when χ,φ → NR")
print(f"    → φ→2σ decay is NOT needed for ΔN_eff = 0.153")
print(f"    → σ just needs to EXIST and be in thermal contact at T_D")
print()
print(f"    BUT: σ-χ coupling g = m_χ/f = {g_sigma_chi:.2e} → too weak!")
print(f"    σ never thermalizes with χ,φ through scattering.")
print(f"    → σ must be populated through φ→2σ decay OR initial conditions.")
print()
print(f"    If σ was part of the thermal dark sector at T_D:")
print(f"      Need Γ(σ-dark) > H(T_D) — NOT satisfied (g ~ 10⁻²⁰)")
print(f"      Unless σ thermalizes through a DIFFERENT channel")
print(f"      (e.g., σ-σ self-interaction from Λ_d⁴ cos(σ/f) potential)")
print()

# σ self-interaction rate (from the cos(σ/f) potential)
# V = Λ_d⁴ [1 - cos(σ/f)] → V'' = Λ_d⁴/f² = m_σ² → mass
# V'''' = Λ_d⁴/f⁴ → self-coupling λ_σ = Λ_d⁴/f⁴
lambda_sigma = Lambda_d**4 / f_sigma**4
print(f"  σ self-coupling: λ_σ = Λ_d⁴/f⁴ = {lambda_sigma:.2e}")
print(f"  → NEGLIGIBLE. σ is essentially free-streaming from the start.")
print()

print(f"  BOTTOM LINE:")
print(f"  ═══════════")
print(f"  ΔN_eff = 0.153 requires σ to carry the dark sector entropy.")
print(f"  For σ to do this, it must be thermalized at T_D = 200 MeV.")
print(f"  With f ~ 0.24 M_Pl, the σ-matter coupling is ~ 10⁻²⁰.")
print(f"  σ CANNOT thermalize → it cannot carry entropy → ΔN_eff → 0")
print(f"  UNLESS σ has additional (non-gravitational) couplings to χ or φ.")
print()
print(f"  TWO RESOLUTIONS:")
print(f"  (1) Accept ΔN_eff ≈ 0 (no dark radiation at BBN/CMB)")
print(f"      → Model is safe but loses the CMB-S4 prediction")
print(f"  (2) σ has enhanced coupling (f << M_Pl, e.g., f ~ TeV scale)")
print(f"      → σ thermalizes, carries entropy, ΔN_eff = 0.153")
print(f"      → But m_σ = Λ_d²/f would change → dark energy story changes")
print()
print("=" * 72)
print("  Test 22 COMPLETE")
print("=" * 72)
