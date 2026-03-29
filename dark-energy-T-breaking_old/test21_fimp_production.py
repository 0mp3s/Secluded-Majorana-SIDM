"""
Test 21: FIMP Freeze-In Production — Is T_D = 200 MeV consistent?
====================================================================
Core question:
  Test 20 found that κ ~ 5.3×10⁻⁴ (Higgs portal) gives T_D = 200 MeV.
  Does this same κ produce Ω_χ h² = 0.120 via freeze-in?

  If YES → T_D = 200 MeV is not an assumption, it's a DERIVED result.
  If NO  → T_D is still a free parameter; need an additional condition.

FIMP picture:
  1. Reheating at T_RH >> m_h: dark sector starts empty.
  2. SM bath produces φ via portal: h → φφ  (dominant at T ~ m_h).
  3. φ thermalizes within dark sector → sets T_D via entropy balance.
  4. χ thermalizes with φ inside dark sector (α_D ~ 5.7×10⁻³ sufficient).
  5. χ freezes out of DARK sector at x_D = m_χ / T_D ~ 20.
  6. Today:  Ω_χ h² = Ω_χ^{dark} × (ξ³ × g*_S,SM / g*_S,D)

  Temperature ratio:  ξ ≡ T_D / T_SM

  From entropy conservation in separate sectors:
    ξ(T) = ξ_initial × [g*_S,SM(T_prod) / g*_S,SM(T)]^{1/3}
  where ξ_initial is set by the heat injected from SM → dark sector.

Strategy:
  A. Compute the heat Q injected by portal κ into dark sector.
  B. From Q, get T_D / T_SM = ξ.
  C. Compute Ω_χ h² from dark-sector freeze-out at T_D,fo = m_χ/20.
  D. Compare to Ω_DM h² = 0.120.

MAP PARAMETERS:
  m_χ = 94.07 MeV,  m_φ = 11.10 MeV,  α_D = 5.734×10⁻³,  κ = 5.3×10⁻⁴
"""

import numpy as np
from scipy import integrate, interpolate
import warnings
warnings.filterwarnings("ignore")

print("=" * 68)
print("TEST 21: FIMP FREEZE-IN — T_D = 200 MeV CONSISTENCY CHECK")
print("=" * 68)

# ─── Physical constants ──────────────────────────────────────────────────────
G_N          = 6.709e-39       # GeV^{-2}
M_Pl         = 1.0 / np.sqrt(8 * np.pi * G_N)   # reduced Planck mass
hbar_c_GeVfm = 0.1973          # GeV·fm
hbar_c_GeVcm = hbar_c_GeVfm * 1e-13  # GeV·cm
c_cm_s       = 3e10            # cm/s

# Unit conversions
GeV_to_inv_s = 1.0 / 6.582e-25   # 1 GeV = 1/ℏ × s^{-1}
cm3_per_GeV2 = hbar_c_GeVcm**2 * c_cm_s  # GeV^{-2} → cm³/s

# Observed relic density
Omega_DM_h2  = 0.120
Omega_b_h2   = 0.0224
h            = 0.674
rho_crit_h2  = 8.099e-47  # GeV⁴  (ρ_c × h² in natural units)

# SM parameters
m_h          = 125.1          # GeV   Higgs mass
v_EW         = 246.0          # GeV   Higgs vev
m_e          = 0.511e-3       # GeV
m_mu         = 0.1057         # GeV
m_tau        = 1.777          # GeV
m_W          = 80.4           # GeV
m_Z          = 91.2           # GeV
T_RH         = 1000.0         # GeV   reheating temperature (conservative)
T_today      = 2.725 * 8.617e-14  # GeV   CMB temperature

# Dark sector parameters
m_chi        = 94.07e-3       # GeV
m_phi        = 11.10e-3       # GeV
alpha_D      = 5.734e-3
g_s2         = 4 * np.pi * alpha_D
kappa        = 5.3e-4         # Higgs portal coupling from Test 20
T_D_assumed  = 200e-3         # GeV   dark sector temperature (assumption)
x_fo         = 20.0           # m_chi / T_D at dark freeze-out

print(f"\n  Model parameters:")
print(f"    m_χ         = {m_chi*1e3:.2f} MeV")
print(f"    m_φ         = {m_phi*1e3:.2f} MeV")
print(f"    α_D         = {alpha_D:.4e}")
print(f"    κ (portal)  = {kappa:.2e}  (from Test 20)")
print(f"    T_D,assumed = {T_D_assumed*1e3:.0f} MeV  (what we want to test)")
print(f"    T_RH        = {T_RH:.0f} GeV")

# ─── g*(T): effective dof ────────────────────────────────────────────────────
# Tabulated g*_S(T) for SM (key thresholds)
_g_table = np.array([
    # T [GeV],   g*_S
    [1e-4,        3.91],   # after e+e- annihilation
    [0.511e-3,    10.75],  # e+e- threshold
    [1e-2,        10.75],
    [0.15,        17.25],  # muon threshold
    [0.3,         29.0 ],  # QCD crossover
    [1.0,         57.75],
    [80.0,        96.25],  # W,Z threshold
    [200.0,       106.75], # Higgs threshold
    [1e4,         106.75],
])

def g_star_S(T):
    """Effective entropy dof g*_S(T), SM only."""
    return float(np.interp(np.log10(T), np.log10(_g_table[:,0]), _g_table[:,1]))

def g_star_rho(T):
    """Effective energy dof g*_ρ(T) ≈ g*_S(T) for our purposes."""
    return g_star_S(T)

def H_hubble(T):
    """Hubble rate H(T) = √(π²g*/90) × T²/M_Pl."""
    g  = g_star_rho(T)
    return np.sqrt(np.pi**2 * g / 90.0) * T**2 / M_Pl

def entropy_density(T):
    """SM entropy density s = (2π²/45) g*_S T³."""
    return (2 * np.pi**2 / 45.0) * g_star_S(T) * T**3

# ─── Part A: h → φφ decay rate ──────────────────────────────────────────────
print(f"\n{'─'*68}")
print("PART A — Heat injection via h → φφ decay")
print(f"{'─'*68}")

# Portal B:  κ |H|² φ²  →  κ (v + h)² φ²/2
# h → φφ: M = 2κv,  Γ = |M|² × p_φ / (8π m_h²)
# p_φ = √(m_h²/4 − m_φ²)

def Gamma_h_phiphi(kap):
    """Partial width h → φφ from portal coupling κ."""
    if m_h < 2 * m_phi:
        return 0.0
    p_phi = np.sqrt(m_h**2 / 4.0 - m_phi**2)
    M_sq  = (2.0 * kap * v_EW)**2
    return M_sq * p_phi / (8.0 * np.pi * m_h**2)

def Gamma_h_SM():
    """SM total Higgs width (approximate)."""
    return 4.07e-3  # GeV  (PDG 2022)

Gamma_portal = Gamma_h_phiphi(kappa)
Gamma_total  = Gamma_h_SM() + Gamma_portal
BR_portal    = Gamma_portal / Gamma_total

print(f"\n  κ = {kappa:.2e}")
print(f"  Γ(h→φφ) = {Gamma_portal:.3e} GeV")
print(f"  Γ(h→SM) = {Gamma_h_SM():.3e} GeV")
print(f"  BR(h→φφ) = {BR_portal:.3e}  (LHC invisible < 0.11 → {BR_portal/0.11:.2f}× allowed)")

# ─── Part B: Energy injection rate into dark sector ─────────────────────────
print(f"\n{'─'*68}")
print("PART B — Energy density injected into dark sector")
print(f"{'─'*68}")

# Thermal production rate: Γ_prod(T) = n_h(T) × Γ(h→φφ)
# n_h(T) = g_h × (m_h² T / 2π²) × K₂(m_h/T)   (Maxwell-Boltzmann)
# where g_h = 1 (real scalar)

def n_higgs_eq(T):
    """Higgs number density in thermal equilibrium."""
    if T < 1e-5:
        return 0.0
    x_ = m_h / T
    if x_ > 200:
        return 0.0
    from scipy.special import kn
    return 1.0 * m_h**2 * T / (2.0 * np.pi**2) * kn(2, x_)

def rho_injection_rate(T):
    """Energy injection rate into dark sector: ρ̇ = n_h × Γ(h→φφ) × m_h."""
    return n_higgs_eq(T) * Gamma_portal * m_h   # GeV⁴/GeV = GeV³ × H → GeV⁴/s

# Total energy injected per comoving volume:
# ΔΕ_D / s_SM = ∫ (dE/dt)_D / (s T H) dT   [changing integration var T→T]

def integrand_energy(T):
    """dY_E/dT: comoving energy density injected at temperature T."""
    H_ = H_hubble(T)
    s_ = entropy_density(T)
    if H_ < 1e-200 or s_ < 1e-200:
        return 0.0
    return rho_injection_rate(T) / (H_ * T * s_)

print(f"\n  Computing total comoving energy injected: ∫ (dn_h Γ m_h) / (HTs) dT ...")

# Integrate from T_today to T_RH
T_arr  = np.logspace(np.log10(1e-4), np.log10(T_RH), 2000)
dY_arr = np.array([integrand_energy(T) for T in T_arr])

# Use trapezoidal in log space
log_T  = np.log(T_arr)
_trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
total_energy_injected = _trapz(dY_arr * T_arr, log_T)  # × T for dlnT

print(f"  Y_E = ∫ dT [...] = {total_energy_injected:.3e} GeV")

# ─── Part C: Dark sector temperature from injected energy ────────────────────
print(f"\n{'─'*68}")
print("PART C — Dark sector temperature ratio ξ = T_D / T_SM")
print(f"{'─'*68}")

# Dark sector dof: g*_D = 1 (φ, real scalar) + 2 (χ, Majorana fermion × 2 spin)
# For Majorana: g_χ = 2 (spin), but fermion gets 7/8 factor in energy:
# ρ_D = π²/30 × (g*_B + 7/8 × g*_F) × T_D⁴
# g*_D = g_φ + 7/8 × g_χ = 1 + 7/8 × 2 = 1 + 1.75 = 2.75  (both species)

g_dark_B = 1.0       # φ (real scalar)
g_dark_F = 2.0       # χ (Majorana: 2 spin states)
g_star_D = g_dark_B + 7.0/8.0 * g_dark_F   # ≈ 2.75

# Energy density of dark sector at T_D:
# ρ_D = π²/30 × g*_D × T_D⁴
# This equals total_energy_injected × s_SM(T_SM):
# → T_D⁴ = 30/(π²g*_D) × Y_E × s_SM(T_SM)
# This gives T_D at T_SM = T_today (after entropy dilution)

# At T_SM = T_today:
s_SM_today = entropy_density(T_today)
rho_D_today = total_energy_injected * s_SM_today

T_D_today  = (30.0 * rho_D_today / (np.pi**2 * g_star_D))**0.25 if rho_D_today > 0 else 0

# Better: compute ξ at the epoch of χ freeze-out from dark sector
# T_D,fo = m_chi / x_fo
T_D_fo_dark = m_chi / x_fo          # dark sector temperature at χ freeze-out
# SM temperature at that time: T_SM,fo = T_D_fo_dark / ξ

# ξ is conserved (separate entropy in each sector) after dark sector decouples.
# ξ = T_D / T_SM = (s_D/s_SM)^{1/3} × (g*_S,SM / g*_S,D)^{1/3}
# From total energy injected, at production epoch T_prod ~ m_h/3 ~ 40 GeV:
T_prod = m_h / 3.0  # peak of Higgs population (thermal)
g_sm_prod  = g_star_S(T_prod)
s_SM_prod  = entropy_density(T_prod)

# T_D at production = (30/π² g*_D × Y_E × s_SM_prod / π² × 30...)
# Actually: ρ_D(T_prod) = Y_E × s_SM(T_prod)
rho_D_prod = total_energy_injected * s_SM_prod
T_D_prod   = (30.0 * rho_D_prod / (np.pi**2 * g_star_D))**0.25 if rho_D_prod > 0 else 0
xi_at_prod = T_D_prod / T_prod if T_prod > 0 else 0

# After production, ξ scales as: ξ(T) = ξ_0 × (g*_S,SM(T_prod)/g*_S,SM(T))^{1/3}
def xi_at_T(T):
    return xi_at_prod * (g_star_S(T_prod) / g_star_S(T))**(1.0/3.0)

xi_at_fo   = xi_at_T(T_D_fo_dark / xi_at_prod)  # iterative → use xi_at_prod as 0th order
# Proper self-consistent: T_SM_fo = T_D_fo_dark / ξ, check consistency
# Simple iteration:
xi_curr = xi_at_prod
for _ in range(5):
    T_SM_fo = T_D_fo_dark / xi_curr
    xi_curr = xi_at_T(T_SM_fo)

xi_fo      = xi_curr
T_SM_fo    = T_D_fo_dark / xi_fo

print(f"\n  Dark sector dof: g*_D = {g_star_D:.2f} (1 scalar + 7/8 × 2 Majorana)")
print(f"  Peak production at T_SM ~ m_h/3 = {T_prod:.1f} GeV")
print(f"  g*_S,SM at production     = {g_sm_prod:.1f}")
print(f"  Y_E (comoving energy)     = {total_energy_injected:.3e} GeV")
print(f"  T_D at production         = {T_D_prod*1e3:.1f} MeV")
print(f"  ξ at production           = T_D/T_SM = {xi_at_prod:.3e}")
print(f"\n  At dark χ freeze-out:")
print(f"    T_D,fo = m_χ/20        = {T_D_fo_dark*1e3:.2f} MeV")
print(f"    ξ(T_fo)                = {xi_fo:.3e}")
print(f"    T_SM,fo                = {T_SM_fo*1e3:.2f} MeV")

# ─── Part D: Relic density from dark sector freeze-out ───────────────────────
print(f"\n{'─'*68}")
print("PART D — χ relic density from dark-sector freeze-out")
print(f"{'─'*68}")

# Standard relic density from freeze-out in dark sector:
# Ω_χ h² = (Ω_χ,dark)_0 × (ξ_fo³ × g*_S,SM,fo / g*_S,D)
#
# In the dark sector alone (ignoring SM coupling):
# <σv>_{χχ→φφ} = π α_D² / (4 m_χ²)   [s-wave, from PI-7 calculation]
sigma_v_dark = np.pi * alpha_D**2 / (4.0 * m_chi**2) * cm3_per_GeV2  # cm³/s

# Standard freeze-out formula in dark sector:
# Y_χ,∞ = (45/π³) × (x_fo / g*_S,D) × H(T_D,fo) / (σv × m_chi³)
# Simplified (Kolb & Turner):
#   Ω h² = (1.07×10⁹ × x_fo) / (g*^{1/2} × M_Pl × σv)  [GeV units]
g_star_eff_fo = g_star_D   # at dark sector freeze-out

# But we need to account for the dark sector Hubble:
# H_D(T_D) = H_SM(T_SM) = H_SM(T_D/ξ)
# At T_D,fo = m_chi/20:
H_D_at_fo = H_hubble(T_SM_fo)   # SM drives Hubble rate
n_chi_eq_fo = (m_chi**2 * T_D_fo_dark / (2 * np.pi**2)) * np.exp(-x_fo)  # MB approx

# Freeze-out yield Y_fo = n_chi / s_D
s_D_fo = (2 * np.pi**2 / 45.0) * g_star_D * T_D_fo_dark**3
Y_chi_dark = n_chi_eq_fo / s_D_fo   # equilibrium value at freeze-out (decoupling)

# Relic comoving abundance today, accounting for ξ:
# The dark sector entropy today: s_D,0 = s_D,fo × (a_fo/a_0)³
# The SM entropy today: s_SM,0 = s_SM,fo × (a_fo/a_0)³ × (g*_S,SM,fo/g*_S,SM,0)
# → n_χ,0/s_SM,0 = Y_χ × (s_D,fo/s_SM,fo) × (g*_S,SM,fo/g*_S,SM,0)^... 

# Simple approach:
# s_D / s_SM = (g*_S,D / g*_S,SM) × ξ³
# so n_χ,0 / s_SM,0 = Y_chi_dark × s_D,0/s_SM,0 = Y_chi_dark × (g*_S,D / g*_S,SM,fo) × ξ_fo³

g_star_S_SM_fo = g_star_S(T_SM_fo)
g_star_S_SM_0  = g_star_S(T_today)

Y_chi_SM  = Y_chi_dark * (g_star_D / g_star_S_SM_fo) * xi_fo**3

# Relic density:
# Ω_χ h² = m_χ × Y_χ,SM × s_SM,0/ρ_crit × h²
# ρ_crit = 3H₀²M_Pl²/(8π) = 8.099×10⁻⁴⁷ GeV⁴ (for h=1, per h² convention)
s_SM_0    = entropy_density(T_today)  # in GeV³
Omega_h2  = m_chi * Y_chi_SM * s_SM_0 / rho_crit_h2

print(f"\n  In dark sector at T_D,fo = {T_D_fo_dark*1e3:.2f} MeV:")
print(f"    <σv>_{'{χχ→φφ}'}     = {sigma_v_dark:.3e} cm³/s  (s-wave, πα_D²/4m_χ²)")
print(f"    x_fo = m_χ/T_D,fo  = {x_fo:.0f}")
print(f"    n_χ^eq(T_fo)       = {n_chi_eq_fo:.3e} GeV³")
print(f"    s_D(T_fo)          = {s_D_fo:.3e} GeV³")
print(f"    Y_χ (dark sector)  = {Y_chi_dark:.3e}")
print(f"\n  Dilution from dark → SM entropy:")
print(f"    g*_S,D             = {g_star_D:.2f}")
print(f"    g*_S,SM at T_fo    = {g_star_S_SM_fo:.1f}")
print(f"    ξ_fo               = {xi_fo:.3e}")
print(f"    Y_χ (SM frame)     = {Y_chi_SM:.3e}")
print(f"\n  Relic density:")
print(f"    Ω_χ h²  (FIMP)     = {Omega_h2:.4f}")
print(f"    Ω_DM h² (Planck)   = {Omega_DM_h2:.3f}")
print(f"    Ratio              = {Omega_h2/Omega_DM_h2:.3f}")

# ─── Part E: ξ scan — what κ gives Ω h² = 0.120? ────────────────────────────
print(f"\n{'─'*68}")
print("PART E — Required ξ for Ω_χ h² = 0.120")
print(f"{'─'*68}")

# Y_chi_SM = Y_chi_dark × (g*_D / g*_S,SM,fo) × ξ³
# Ω h² = m_χ × Y_chi_dark × (g*_D / g*_S,SM,fo) × ξ³ × s_SM_0 / ρ_crit
# → ξ_required³ = Ω_target / (m_χ × Y_chi_dark × g*_D/g*_S × s_0/ρ_crit)
factor        = m_chi * Y_chi_dark * (g_star_D / g_star_S_SM_fo) * s_SM_0 / rho_crit_h2
xi_required_3 = Omega_DM_h2 / factor
xi_required   = xi_required_3**(1.0/3.0)

print(f"\n  For Y_χ (dark sector) = {Y_chi_dark:.3e}")
print(f"  Required ξ_fo for Ω h² = {Omega_DM_h2:.3f}: ξ_required = {xi_required:.3e}")
print(f"  Actual ξ_fo from κ:                           ξ_actual   = {xi_fo:.3e}")
print(f"  Ratio ξ_actual / ξ_required = {xi_fo/xi_required:.3e}")

# Required T_D from ξ_required:
T_D_required = xi_required * T_SM_fo
# What κ gives T_D_required?
# T_D ∝ κ^{1/2} (from energy injection ∝ κ²)
kappa_required = kappa * (xi_required / xi_fo)**(2.0) if xi_fo > 0 else 0

print(f"\n  Implied T_D,prod from ξ_required: {T_D_required*1e3:.2f} MeV")
print(f"  κ that would give Ω h² = 0.120:   κ_required = {kappa_required:.3e}")
print(f"  κ from Test 20 (T_D=200 MeV):     κ_test20   = {kappa:.3e}")

# ─── Part F: κ scan table ────────────────────────────────────────────────────
print(f"\n{'─'*68}")
print("PART F — κ scan: Ω_χ h² vs portal coupling")
print(f"{'─'*68}")

print(f"\n  {'κ':>12}  {'Y_E':>12}  {'ξ_fo':>10}  {'T_D,prod [MeV]':>16}  {'Ω_χ h²':>10}  Status")
print(f"  {'─'*12}  {'─'*12}  {'─'*10}  {'─'*16}  {'─'*10}  {'─'*20}")

for kap in [1e-5, 1e-4, 2e-4, 5e-4, kappa, 1e-3, 2e-3, 5e-3, 1e-2]:
    Gam = Gamma_h_phiphi(kap)
    # Y_E ∝ κ²
    Y_E_k = total_energy_injected * (kap / kappa)**2
    rho_D_k = Y_E_k * s_SM_prod
    T_D_k   = (30.0 * rho_D_k / (np.pi**2 * g_star_D))**0.25 if rho_D_k > 0 else 0
    xi_k    = T_D_k / T_prod
    # ξ at freeze-out (self-consistent, use simplified version)
    xi_k_fo = xi_k * (g_star_S(T_prod) / g_star_S(T_SM_fo))**(1.0/3.0)
    Y_SM_k  = Y_chi_dark * (g_star_D / g_star_S_SM_fo) * xi_k_fo**3
    Om_k    = m_chi * Y_SM_k * s_SM_0 / rho_crit_h2
    tick = "✅" if abs(Om_k - Omega_DM_h2) / Omega_DM_h2 < 0.5 else ("← Test20" if abs(kap-kappa)/kappa < 0.05 else "")
    print(f"  {kap:>12.2e}  {Y_E_k:>12.2e}  {xi_k_fo:>10.2e}  {T_D_k*1e3:>16.2f}  {Om_k:>10.4f}  {tick}")

# ─── SUMMARY ─────────────────────────────────────────────────────────────────
print()
print("=" * 68)
print("TEST 21 SUMMARY — FIMP FREEZE-IN")
print("=" * 68)
print()
print(f"  κ from Test 20            = {kappa:.2e}  (gives T_D = 200 MeV)")
print(f"  Ω_χ h² from this κ        = {Omega_h2:.4f}")
print(f"  Target Ω_DM h²            = {Omega_DM_h2:.3f}")
print(f"  κ for correct relic dens. = {kappa_required:.3e}")
print()

if abs(Omega_h2 - Omega_DM_h2) / Omega_DM_h2 < 0.3:
    print("  ┌──────────────────────────────────────────────────────────────────┐")
    print("  │  ✅ CONSISTENT — κ from Test 20 gives Ω h² ≈ 0.120             │")
    print("  │  T_D = 200 MeV is NOT a free assumption:                        │")
    print("  │  the SAME κ that sets T_D also gives the correct relic density. │")
    print("  └──────────────────────────────────────────────────────────────────┘")
else:
    ratio = Omega_h2 / Omega_DM_h2
    log_gap = np.log10(ratio) if ratio > 0 else 0
    print("  ┌──────────────────────────────────────────────────────────────────┐")
    print(f"  │  ✗ INCONSISTENT — Ω_χ h² = {Omega_h2:.3f}  (target: {Omega_DM_h2:.3f})        │")
    print(f"  │  Overproduction factor: {ratio:.2e}  (log₁₀ gap = {log_gap:.1f})       │")
    print(f"  │                                                                  │")
    print(f"  │  Physical reason:                                                │")
    if ratio > 1:
        print(f"  │  • ξ_fo = {xi_fo:.2e}: dark sector is effectively as warm as SM.    │")
        print(f"  │  • Dark sector thermalizes → Ω_χ h² ≈ standard freeze-out.          │")
        print(f"  │  • This is NOT FIMP — it's thermal freeze-out in hidden sector.     │")
    else:
        print(f"  │  • Under-production: κ too weak to heat dark sector sufficiently.   │")
    print(f"  │                                                                  │")
    print(f"  │  INTERPRETATION:                                                 │")
    print(f"  │  κ sets T_D (temperature ratio ξ) correctly for T_D=200 MeV.   │")
    print(f"  │  But Ω_χ h² is set by DARK FREEZE-OUT at x_D=m_χ/T_D~20,      │")
    print(f"  │  NOT by FIMP from SM.  The dark sector is thermally populated:   │")
    print(f"  │    <σv>_{'{χχ→φφ}'} = πα_D²/4m_χ² ~ {sigma_v_dark:.1e} cm³/s               │")
    print(f"  │    Compare to Planck target: 3×10⁻²⁶ cm³/s                     │")
    if sigma_v_dark > 3e-26:
        print(f"  │    <σv> > Planck target → needs Sommerfeld suppression at T_fo,    │")
    else:
        print(f"  │    <σv> < Planck target → under-abundant without enhancement.      │")
    print(f"  └──────────────────────────────────────────────────────────────────┘")
print()
print(f"  CHAIN STATUS:")
print(f"   ✅ σ/m(30 km/s) = 1.81 cm²/g  [SIDM window ✓]")
print(f"   ✅ κ = {kappa:.2e} → T_D = 200 MeV  [portal coupling ✓]")
print(f"   → Test 21: dark sector thermal relic → Ω_χ h²")
print(f"   → Need to check if <σv> at dark freeze-out matches Planck")
print()
print(f"  <σv>_{'{χχ→φφ}'} (s-wave)      = {sigma_v_dark:.3e} cm³/s")
print(f"  Planck relic target            = 3.0×10⁻²⁶ cm³/s")
print(f"  Ratio                          = {sigma_v_dark/3e-26:.3f}")
gap_orders = np.log10(sigma_v_dark / 3e-26)
if gap_orders < 1:
    print(f"  Status: ✅ <σv> ~ Planck relic target (log₁₀ ratio = {gap_orders:.2f})")
else:
    print(f"  Status: ✗ <σv> >> Planck target by 10^{gap_orders:.0f} orders")
    print(f"  → In standard thermal freeze-out: Ω_χ h² ~ {0.12 * 3e-26/sigma_v_dark:.2e} (UNDER-ABUNDANT)")
    print(f"  → SIDM σ/m = 1.81 cm²/g needs Sommerfeld enhancement S ~ {sigma_v_dark/3e-26:.0e} at v=30 km/s")
    print(f"  → Relic density cannot be computed from bare πα_D²/4m_χ² alone;")
    print(f"     full Sommerfeld-enhanced Boltzmann integration required.")
print()
print("=" * 68)
print("Test 21 COMPLETE")
print("=" * 68)
