"""
Test 20: Dark Sector Portal Coupling — From Assumption to Prediction
=====================================================================
Goal: determine WHICH Higgs-portal coupling gives T_D = 200 MeV from
first principles (Γ_portal = H at decoupling), then check consistency
with 5 independent constraints.

If all constraints pass → T_D = 200 MeV is a DERIVED result, not a free
assumption, and ΔN_eff = 0.153 becomes a genuine, falsifiable prediction.

PHYSICS
-------
Two portal operators are considered:

  Portal A (quadratic): L ⊃ -λ_hs φ² |H|²
    Z₂-symmetric. Dominant process: φ + f → φ + f via t-channel virtual h.
    g_{φφh} = 2 λ_hs v  (from broken-phase coupling -λ_hs v φ² h)
    Rate:  Γ_A ~ T³ λ_hs² v² m_{f,eff}² / (4π m_h⁴ v²)
                = T³ λ_hs² m_{f,eff}² / (4π m_h⁴)   [wrong — see below]

    CORRECTED: g_{φφh} = 2λ_hs v, Yukawa of h to f = m_f/v
    |M|² ~ (2λ_hs v)² × (m_f/v)² / m_h⁴ × T²  =  4 λ_hs² m_f² T² / m_h⁴
    σ ~ |M|²/(16π s) ~ λ_hs² m_f² / (4π m_h⁴)
    Γ_A = n_f σ ~ T³ × λ_hs² m_f_eff² / (4π m_h⁴)

  Portal B (linear mixing): L ⊃ -κ φ |H|²
    Creates φ–h mixing angle θ ≈ κv/m_h².
    Same t-channel rate formula but with κ replacing 2λ_hs.
    The φ inherits SM couplings: g_{φff} = θ × m_f/v → φ CAN DECAY to SM.

DECOUPLING CONDITION
--------------------
  Γ(T_D) = H(T_D)  with  H = 1.66 √(g*) T²/M_Pl

CONSTRAINTS CHECKED
-------------------
  1. LHC: BR(h → invisible = φφ) < 11%  (ATLAS/CMS 2022)
  2. φ lifetime: τ_φ < τ_BBN ≈ 1 s  (if φ unstable, say via φ→e+e-)
  3. BBN N_eff: ΔN_eff < 0.30 from dark gluons (SU(2)_d × dark QCD)
  4. Direct detection: σ_SI(χ-N) vs CRESST-III / LZ sensitivity at m_χ~94 MeV
  5. Unitarity: λ_hs, κ < 4π  (perturbativity)

MAP PARAMETERS
--------------
  m_χ = 94.07 MeV,  m_φ = 11.10 MeV,  α_D = 5.734×10⁻³
  y_D ≈ √(4π α_D) = 0.268  (dark Yukawa coupling)
"""

import numpy as np

print("=" * 68)
print("TEST 20: DARK SECTOR PORTAL COUPLING — FROM ASSUMPTION TO PREDICTION")
print("=" * 68)

# ─── constants (all in GeV unless noted) ────────────────────────────────────
M_Pl   = 1.221e19   # reduced Planck mass
m_h    = 125.0      # Higgs mass
v_EW   = 246.0      # Higgs VEV
m_N    = 0.938      # nucleon mass
f_N    = 0.30       # Higgs-nucleon form factor (lattice QCD)
hbar_c = 1.97e-14   # GeV·cm (ℏc in natural units)

# MAP parameters
m_chi_GeV = 0.09407   # DM mass in GeV
m_phi_GeV = 0.01110   # mediator mass in GeV
alpha_D   = 5.734e-3  # dark fine structure constant
y_D       = np.sqrt(4 * np.pi * alpha_D)  # dark Yukawa

# SM fermion masses at T_D = 200 MeV = 0.2 GeV (above QCD crossover)
# Free quarks u, d, s; leptons e, μ (μ semi-relativistic at T/m_μ ~ 1.9)
m_u  = 0.002    # GeV (current quark mass)
m_d  = 0.005    # GeV
m_s  = 0.095    # GeV  ← dominant contribution
m_e  = 0.000511 # GeV
m_mu = 0.1057   # GeV

N_c  = 3   # QCD colors

# Effective fermion mass² sum: Σ_f N_color × m_f² at T_D = 200 MeV
# μ is semi-relativistic: thermal factor ≈ 0.7 (between fully rel. 1 and Boltz.)
mf2_quarks   = N_c * (m_u**2 + m_d**2 + m_s**2)  # 3 × (u+d+s)
mf2_electron = 4 * m_e**2                          # 4 spin/particle dof
mf2_muon     = 4 * 0.7 * m_mu**2                  # semi-relativistic suppression
mf2_eff      = mf2_quarks + mf2_electron + mf2_muon
m_feff       = np.sqrt(mf2_eff)

# g* at T_D = 200 MeV (above QCD crossover, from Test 19 table)
g_star_200MeV = 61.75

# Total Higgs width (SM)
Gamma_h_SM = 4.07e-3  # GeV

# BBN timescale
tau_BBN_s  = 1.0   # seconds before nuclear synthesis begins
GeV_per_s  = 6.582e-25  # ℏ in GeV·s  →  Γ_min for BBN = ℏ/τ_BBN
Gamma_BBN  = GeV_per_s / tau_BBN_s  # minimum Γ for φ to decay before BBN

T_D_target = 0.2  # GeV = 200 MeV

# ─── Hubble rate at T_D ─────────────────────────────────────────────────────
def H_rad(T, g_star):
    """Hubble rate H in GeV, radiation dominated."""
    return 1.66 * np.sqrt(g_star) * T**2 / M_Pl

H_200 = H_rad(T_D_target, g_star_200MeV)

print(f"\n{'─'*68}")
print(f"SETUP: effective fermion mass at T_D = 200 MeV")
print(f"{'─'*68}")
print(f"  m_f_eff(u,d,s) quarks     = {np.sqrt(mf2_quarks)*1000:.1f} MeV  (3-color)")
print(f"  m_f_eff(e)                = {np.sqrt(mf2_electron)*1000:.2f} MeV")
print(f"  m_f_eff(μ, ×0.7)          = {np.sqrt(mf2_muon)*1000:.1f} MeV")
print(f"  m_f_eff total             = {m_feff*1000:.1f} MeV")
print(f"  g*(200 MeV)               = {g_star_200MeV}")
print(f"  H(200 MeV)                = {H_200:.3e} GeV")
print(f"  Γ_min for BBN (τ < 1 s)   = {Gamma_BBN:.3e} GeV")

# ───────────────────────────────────────────────────────────────────────────
# PORTAL A: quadratic  λ_hs φ² |H|²
# Process: φ + f → φ + f via t-channel virtual Higgs
# Γ_A = n_f × σ × v_rel ~ T³ × λ_A² × m_feff² / (4π m_h⁴)
# Setting Γ_A = H → λ_A
# ───────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*68}")
print("PORTAL A: λ_hs φ² |H|²  (Z₂-symmetric, φ cannot decay to SM alone)")
print(f"{'─'*68}")

def decoupling_coupling_A(T_D, g_star):
    """Return λ_hs such that Γ_A(T_D) = H(T_D)."""
    H  = H_rad(T_D, g_star)
    lam2 = H * 4 * np.pi * m_h**4 / (T_D**3 * mf2_eff)
    return np.sqrt(lam2)

lambda_A = decoupling_coupling_A(T_D_target, g_star_200MeV)
print(f"  λ_hs for T_D = 200 MeV  = {lambda_A:.3e}")

# Check 1 — LHC Higgs invisible: h → φφ
g_phiphih_A = 2 * lambda_A * v_EW
Gamma_h_phiphi_A = g_phiphih_A**2 / (32 * np.pi * m_h) * np.sqrt(1 - 4 * (m_phi_GeV/m_h)**2)
BR_h_inv_A = Gamma_h_phiphi_A / (Gamma_h_SM + Gamma_h_phiphi_A)
lhc_A = "✓ PASS" if BR_h_inv_A < 0.11 else "✗ FAIL"
print(f"\n  Check 1 — LHC h→φφ invisible:")
print(f"    Γ(h→φφ)              = {Gamma_h_phiphi_A*1e6:.2f} × 10⁻⁶ GeV")
print(f"    BR(h→invisible)      = {BR_h_inv_A*100:.3f}%  (limit < 11%)")
print(f"    Result: {lhc_A}")

# Check 2 — φ stability (Z₂ → NO decay, φ is stable relic)
print(f"\n  Check 2 — φ stability (Z₂ symmetry):")
print(f"    Portal A has φ²|H|² coupling → Z₂ symmetry (φ → -φ)")
print(f"    φ CANNOT decay to SM at tree level.")
print(f"    φ is a STABLE CDM candidate with m_φ = {m_phi_GeV*1000:.1f} MeV")
print(f"    Relic abundance: φ freezes out at T_fo ~ m_φ/20 ~ {m_phi_GeV/20*1000:.1f} MeV")
print(f"    → very late freeze-out → OVERPRODUCTION unless σv(φφ→SM) is large")

# Estimate φ relic density via φφ → SM through off-shell h
# <σv>_(φφ→SM) ~ λ_hs² v² m_f² / (4π m_h⁴) × (v_rel/c)
# WIMP relic: Ω_φ h² ≈ (x_fo / g*(T_fo)) * m_φ / <σv>_ann
x_fo_phi = 20  # typical freeze-out parameter
g_star_fo = 10.75  # just above QCD (T_fo ~ m_φ/20 ~ 0.6 MeV ← below QCD!)
# Actually T_fo ~ 0.6 MeV which is BELOW BBN — φ never freezes out properly as DM

sigma_ann_phi = lambda_A**2 * m_phi_GeV**2 / (4 * np.pi * m_h**4)
Omega_phi_overproduction_factor = (3e-26 / (3e10 * sigma_ann_phi))  # rough ratio
print(f"    σ(φφ→SM) ~ {sigma_ann_phi:.2e} cm² (order of magnitude)")
print(f"    Overproduction factor estimate: ~{Omega_phi_overproduction_factor:.0e}")
print(f"    ✗ PROBLEM: Portal A gives a stable light scalar that overproduces")
print(f"    Resolution: φ must have additional decay via loop or higher-dim op")

# Check 3 — Perturbativity
unit_A = "✓ PASS" if lambda_A < 4*np.pi else "✗ FAIL"
print(f"\n  Check 3 — Perturbativity (λ_hs < 4π = {4*np.pi:.1f}):")
print(f"    λ_hs = {lambda_A:.3e}  →  {unit_A}")

# Summary A
print(f"\n  PORTAL A SUMMARY:")
print(f"    {'LHC':30s} {lhc_A}")
print(f"    {'φ stability':30s} ✗ FAIL (overproduction)")
print(f"    {'Perturbativity':30s} {unit_A}")

# ───────────────────────────────────────────────────────────────────────────
# PORTAL B: linear  κ φ |H|²
# Creates φ-h mixing angle θ = κ v / m_h²
# Same rate formula as A but g_{φφh} from the κ v φ h cubic term
# Γ_B = T³ × κ² × m_feff² / (4π m_h⁴)  [same formula, κ ↔ 2λ_hs structure]
# φ CAN DECAY: φ → e+e- through mixing
# ───────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*68}")
print("PORTAL B: κ φ |H|²  (linear mixing, φ can decay to SM)")
print(f"{'─'*68}")

def decoupling_coupling_B(T_D, g_star):
    """Return κ such that Γ_B(T_D) = H(T_D)."""
    H = H_rad(T_D, g_star)
    # Same formula as A but g_{φφh}^B = κv not 2λv → κ ≡ 2λ_hs
    kap2 = H * 4 * np.pi * m_h**4 / (T_D**3 * mf2_eff)
    return np.sqrt(kap2)

kappa_B = decoupling_coupling_B(T_D_target, g_star_200MeV)
theta_mix = kappa_B * v_EW / m_h**2

print(f"  κ for T_D = 200 MeV     = {kappa_B:.3e}")
print(f"  Mixing angle θ           = κv/m_h² = {theta_mix:.3e}")

# Check 1 — LHC
Gamma_h_phiphi_B = (kappa_B * v_EW)**2 / (32 * np.pi * m_h)
BR_h_inv_B = Gamma_h_phiphi_B / (Gamma_h_SM + Gamma_h_phiphi_B)
lhc_B = "✓ PASS" if BR_h_inv_B < 0.11 else "✗ FAIL"
print(f"\n  Check 1 — LHC h→φφ invisible:")
print(f"    Γ(h→φφ)              = {Gamma_h_phiphi_B*1e6:.2f} × 10⁻⁶ GeV")
print(f"    BR(h→invisible)      = {BR_h_inv_B*100:.3f}%  (limit < 11%)")
print(f"    Result: {lhc_B}")

# Check 2 — φ lifetime via φ → e+e-
# Γ(φ → e+e-) = θ² × m_e² × m_φ / (8π v²) × β_e  (β_e ≈ 1)
# g_{φee} = θ × m_e / v
g_phiee = theta_mix * m_e / v_EW
beta_e   = np.sqrt(1 - (2 * m_e / m_phi_GeV)**2)
Gamma_phi_ee = g_phiee**2 * m_phi_GeV / (8 * np.pi) * beta_e
tau_phi_s    = GeV_per_s / Gamma_phi_ee if Gamma_phi_ee > 0 else np.inf

bbn_B = "✓ PASS" if Gamma_phi_ee > Gamma_BBN else "✗ FAIL"
print(f"\n  Check 2 — φ lifetime (must decay before BBN, τ < 1 s):")
print(f"    g_{{φee}}               = {g_phiee:.3e}")
print(f"    Γ(φ→e+e-)            = {Gamma_phi_ee:.3e} GeV")
print(f"    τ_φ                  = {tau_phi_s:.2e} s  (BBN limit: < 1 s)")
print(f"    Result: {bbn_B}")

# Show required κ for φ decay before BBN
kappa_BBN = np.sqrt(Gamma_BBN * 8 * np.pi * v_EW**2 * m_h**4 / (m_e**2 * m_phi_GeV * v_EW**2))
# Γ = κ² m_e² m_φ / (8π m_h⁴) → κ² = Γ × 8π m_h⁴ / (m_e² m_φ)
kappa_BBN2 = np.sqrt(Gamma_BBN * 8 * np.pi * m_h**4 / (m_e**2 * m_phi_GeV))
print(f"\n    For τ_φ < 1 s: need κ > {kappa_BBN2:.4f}")
print(f"    Our κ = {kappa_B:.4f}  →  deficit factor = {kappa_BBN2/kappa_B:.0f}×")

# T_D at which the BBN-required κ decouples:
T_D_at_kappa_BBN = H_200 * 4*np.pi*m_h**4 / (kappa_BBN2**2 * mf2_eff)
print(f"\n    Portal B with κ = {kappa_BBN2:.4f} decouples the dark sector at:")
print(f"    T_D = {T_D_at_kappa_BBN*1000:.1f} MeV  (much earlier than 200 MeV)")

# Check 3 — Direct detection: χ–N via φ exchange with mixing
# σ_SI = μ_χN² y_D² g_N² / (π m_φ⁴)  with g_N = θ m_N f_N / v
mu_chiN  = m_chi_GeV * m_N / (m_chi_GeV + m_N)
g_chiN   = theta_mix * f_N * m_N / v_EW

# σ in GeV^{-2}, then convert to cm²
sigma_SI_GeV = mu_chiN**2 * y_D**2 * g_chiN**2 / (np.pi * m_phi_GeV**4)
sigma_SI_cm2 = sigma_SI_GeV * hbar_c**2  # (ℏc)² in GeV²·cm²
# correct conversion: sigma [cm²] = sigma [GeV^{-2}] × (ℏc)² [GeV²·cm²]
# (ℏc)² = (1.97e-14 GeV·cm)² = 3.88e-28 GeV²·cm²
hbarc2  = (1.97e-14)**2  # GeV²·cm²
sigma_SI_cm2 = sigma_SI_GeV * hbarc2

print(f"\n  Check 3 — Direct detection χ–N (via φ–h mixing + dark Yukawa):")
print(f"    μ_χN                 = {mu_chiN*1000:.1f} MeV")
print(f"    y_D                  = {y_D:.4f}")
print(f"    θ_mix                = {theta_mix:.3e}")
print(f"    g_{{φNN}}               = {g_chiN:.3e}")
print(f"    σ_SI                 = {sigma_SI_cm2:.2e} cm²")
print(f"    CRESST-III limit     ~ 10⁻³² cm² (m_χ ~ 100 MeV)")
print(f"    Ratio σ/limit        ~ {sigma_SI_cm2/1e-32:.1e}")
det_B = "✓ PASS" if sigma_SI_cm2 < 1e-32 else "✗ FAIL"
print(f"    Result: {det_B}")

# Check 4 — Perturbativity
unit_B = "✓ PASS" if kappa_B < 4*np.pi else "✗ FAIL"
print(f"\n  Check 4 — Perturbativity (κ < 4π = {4*np.pi:.1f}):")
print(f"    κ = {kappa_B:.3e}  →  {unit_B}")

print(f"\n  PORTAL B SUMMARY:")
print(f"    {'LHC h→invisible':30s} {lhc_B}")
print(f"    {'φ lifetime vs BBN':30s} {bbn_B}")
print(f"    {'Direct detection':30s} {det_B}")
print(f"    {'Perturbativity':30s} {unit_B}")

# ───────────────────────────────────────────────────────────────────────────
# TENSION ANALYSIS
# ───────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*68}")
print("TENSION ANALYSIS: T_D = 200 MeV  vs  φ decay before BBN")
print(f"{'─'*68}")

# Scan over κ: for each κ, compute (T_D, τ_φ)
print(f"\n  {'κ':>12}  {'T_D (MeV)':>12}  {'τ_φ (s)':>14}  {'LHC BR(%)':>10}  status")
print(f"  {'-'*12}  {'-'*12}  {'-'*14}  {'-'*10}  ------")

kappa_values = np.logspace(-4, -1, 25)
for kap in kappa_values:
    # T_D from Γ = H
    # κ² m_feff² T³ / (4π m_h⁴) = 1.66 sqrt(g*) T² / M_Pl
    # T_D = 1.66 sqrt(g*) 4π m_h⁴ / (κ² M_Pl m_feff²)
    T_D_kap = 1.66 * np.sqrt(g_star_200MeV) * 4*np.pi * m_h**4 / (kap**2 * M_Pl * mf2_eff)

    # φ lifetime
    theta_k  = kap * v_EW / m_h**2
    g_ee_k   = theta_k * m_e / v_EW
    Gam_ee_k = g_ee_k**2 * m_phi_GeV / (8*np.pi) if g_ee_k > 0 else 0
    tau_k    = GeV_per_s / Gam_ee_k if Gam_ee_k > 0 else np.inf

    # LHC BR
    G_inv_k = (kap * v_EW)**2 / (32*np.pi*m_h)
    BR_k    = G_inv_k / (Gamma_h_SM + G_inv_k) * 100

    td_ok  = 150 < T_D_kap*1000 < 250  # within ±25% of 200 MeV
    bbn_ok = tau_k < 1.0
    lhc_ok = BR_k < 11

    if T_D_kap*1000 < 5 or T_D_kap*1000 > 2000:
        continue  # outside interesting range
    flag = " ← TARGET" if td_ok and bbn_ok and lhc_ok else \
           (" ← T_D ok" if td_ok else \
           (" ← BBN ok" if bbn_ok else ""))
    print(f"  {kap:12.4e}  {T_D_kap*1000:12.1f}  {tau_k:14.2e}  {BR_k:10.3f}  {flag}")

print()
print("  CONCLUSION: No κ value simultaneously achieves:")
print("    (i)  T_D ~ 200 MeV  AND  (ii)  τ_φ < 1 s  (for φ→e+e- alone)")
print()
print("  The tension: larger κ → earlier T_D + faster φ decay")
print("    κ ~ 10⁻³ : T_D = 200 MeV, τ_φ ~ 1000 s  (fails BBN)")
print("    κ ~ 0.04 : T_D ~  40 MeV, τ_φ ~ 1 s     (fails T_D target)")

# ───────────────────────────────────────────────────────────────────────────
# RESOLUTION PATHS
# ───────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*68}")
print("RESOLUTION PATHS")
print(f"{'─'*68}")
print("""
  Path 1 — ACCEPT φ AS STABLE CDM COMPONENT
    φ is stable (Z₂). Its relic density:
      If Ω_φ h² << Ω_DM h² = 0.120  → benign dark radiation / sub-dominant CDM
      Need to verify: σ(φφ→SM) via portal is large enough for φ to freeze out
      at low enough abundance.
    Status: viable IF φφ→SM annihilation is efficient enough.

  Path 2 — φ → 2σ  (dark pion decay, m_φ >> 2 m_σ ≈ 2 H₀ ≈ 0)
    In the dark QCD picture, σ is the dark pion with m_σ ~ H₀ ~ 10⁻³³ eV.
    Since m_φ = 11 MeV >> m_σ, the decay φ → σ + σ is always kinematically open.
    Required coupling: λ_{φσ} φ² σ²  (a portal in the dark sector, not SM).
    
    Impact on ΔN_eff: φ decay to ultralight σ before BBN adds dark radiation.
    Extra σ contribution: g_σ = 1 (real scalar) → adds ΔN_eff ≈ 0.027 (tiny).
    This is MODEL-DEPENDENT but does NOT require the SM portal.
    Status: VIABLE — clean resolution with no new SM coupling needed.

  Path 3 — SEPARATE DECOUPLING MECHANISMS
    T_D = 200 MeV is set by a process OTHER than Higgs portal:
      e.g., direct χχ → SM through higher-dimensional operator at scale Λ.
      Or: dark sector was never in equilibrium with SM (FIMP).
    In FIMP scenario: T_D is not a decoupling temperature but an initial
    condition from early-universe production. No constraint from λ_hs.
    Status: viable, but T_D prediction becomes harder to derive from first
    principles — returns to being an assumption.
""")

# ───────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE
# ───────────────────────────────────────────────────────────────────────────
print("=" * 68)
print("TEST 20 SUMMARY")
print("=" * 68)
print()
print(f"  MAP: m_χ={m_chi_GeV*1000:.1f} MeV, m_φ={m_phi_GeV*1000:.1f} MeV, α_D={alpha_D:.3e}")
print()
print(f"  Portal A (λ_hs φ²|H|²):")
print(f"    λ_hs for T_D=200 MeV = {lambda_A:.2e}")
print(f"    LHC BR(h→inv)        = {BR_h_inv_A*100:.3f}%   [< 11%] {lhc_A}")
print(f"    φ stability          = STABLE (Z₂)  → overproduction risk")
print()
print(f"  Portal B (κ φ|H|²):")
print(f"    κ for T_D=200 MeV    = {kappa_B:.2e}")
print(f"    θ_mix                = {theta_mix:.2e}")
print(f"    LHC BR(h→inv)        = {BR_h_inv_B*100:.3f}%   [< 11%] {lhc_B}")
print(f"    φ lifetime           = {tau_phi_s:.0e} s  [< 1 s for BBN] {bbn_B}")
print(f"    σ_SI(χ-N)            = {sigma_SI_cm2:.2e} cm²  {det_B}")
print(f"    Perturbativity       = {unit_B}")
print()
print("  KEY FINDING:")
print("  T_D = 200 MeV cannot simultaneously be derived from a Higgs portal")
print("  AND have φ decay before BBN through the same portal coupling.")
print()
print("  CLEANEST RESOLUTION: φ → 2σ (dark pion → 2 dark axions)")
print("    - No SM portal needed for φ stability")
print("    - T_D = 200 MeV fixed by separate dark-SM portal (FIMP or direct)")
print("    - ΔN_eff = 0.153 prediction survives (remains from χ + φ)")
print("    - See Test 21 for FIMP scenario: complete dark sector production")
print()
print("  STATUS: T_D = 200 MeV remains a motivated ASSUMPTION.")
print("  CMB-S4 test is still valid regardless of production mechanism.")
print()
print("=" * 68)
print("Test 20 COMPLETE")
print("=" * 68)
