"""
Dark Force Accumulation Test
============================
Omer's idea: each DM particle radiates a tiny φ/σ field.
In a uniform DM background, the fields accumulate.
Does the stored energy = ρ_Λ?

Physics:
  Static Poisson eq: (∇² - m²)φ = -g n_χ
  Uniform solution:  φ_0 = g n_χ / m²
  Energy density:    ρ = ½ m² φ_0² = g² n_χ² / (2 m²)

Compare to ρ_Λ for both φ (MeV mediator) and σ (dark axion).
"""
import numpy as np

# ============= Constants (natural units, GeV) ==============
M_Pl   = 2.435e18       # reduced Planck mass [GeV]
H_0    = 1.44e-42        # Hubble constant [GeV]  (67.4 km/s/Mpc)
rho_crit = 3 * H_0**2 * M_Pl**2 / (8 * np.pi)  # ~ 8.1e-47 GeV^4
rho_DM = 0.26 * rho_crit
rho_Lambda = 0.69 * rho_crit

hbar_c_cm = 1.9733e-14   # ℏc in GeV·cm
hbar_c_m  = hbar_c_cm * 1e-2  # GeV·m

print("="*65)
print("  DARK FORCE ACCUMULATION: φ and σ vs ρ_Λ")
print("="*65)
print(f"\nρ_crit  = {rho_crit:.2e} GeV⁴")
print(f"ρ_DM    = {rho_DM:.2e} GeV⁴")
print(f"ρ_Λ     = {rho_Lambda:.2e} GeV⁴")

# ============= MAP Benchmark Point ==============
m_chi  = 0.094    # GeV (94 MeV)
m_phi  = 1e-3     # GeV (1 MeV)
alpha  = 5.734e-3 # coupling (MAP value; corrected from 3.4e-3)
y      = np.sqrt(4 * np.pi * alpha)  # Yukawa coupling
theta  = np.arcsin(1/3)
y_s    = y * np.cos(theta)
y_p    = y * np.sin(theta)
f_sigma = 0.2 * M_Pl  # σ decay constant

n_chi = rho_DM / m_chi  # DM number density [GeV³]

print(f"\n--- MAP Benchmark ---")
print(f"m_χ = {m_chi*1e3:.0f} MeV,  m_φ = {m_phi*1e3:.0f} MeV")
print(f"α = {alpha:.4f},  y = {y:.4f}")
print(f"y_s = {y_s:.4f},  y_p = {y_p:.4f}")
print(f"n_χ = {n_chi:.2e} GeV³")

# ============= (1) φ field from DM background ==============
print("\n" + "="*65)
print("  (1) φ MEDIATOR — accumulated background field")
print("="*65)

# Static Poisson: m_φ² φ₀ = y_s n_χ
phi_0 = y_s * n_chi / m_phi**2
rho_phi = 0.5 * m_phi**2 * phi_0**2

range_phi = 1 / m_phi  # in GeV⁻¹
range_phi_m = range_phi * hbar_c_m  # in meters

print(f"\nφ₀ = y_s × n_χ / m_φ² = {phi_0:.2e} GeV")
print(f"ρ_φ = ½ m_φ² φ₀² = {rho_phi:.2e} GeV⁴")
print(f"ρ_φ / ρ_Λ = {rho_phi/rho_Lambda:.2e}")
print(f"Range of φ = 1/m_φ = {range_phi_m:.2e} m = {range_phi_m*1e15:.1f} fm")

# What m_φ would give ρ_Λ?
m_phi_needed = np.sqrt(y_s**2 * n_chi**2 / (2 * rho_Lambda))
range_needed = 1 / m_phi_needed * hbar_c_m

print(f"\n→ For ρ_φ = ρ_Λ, need m_φ = {m_phi_needed:.2e} GeV = {m_phi_needed*1e9:.2e} eV")
print(f"  Range = {range_needed:.2e} m = {range_needed/1e3:.1f} km = {range_needed/1.5e11:.2e} AU")

# ============= (2) σ field from DM background ==============
print("\n" + "="*65)
print("  (2) σ DARK AXION — accumulated background field")
print("="*65)

# σ couples to DM through CW (1-loop effective coupling)
# dm_χ/dσ ~ y² m_χ / (16π² f)
g_sigma = y**2 * m_chi / (16 * np.pi**2 * f_sigma)

# CW mass of σ
m_sigma_CW = 1e-21  # GeV (~ 10⁻¹² eV, from CW calculation)

sigma_0 = g_sigma * n_chi / m_sigma_CW**2
rho_sigma = 0.5 * m_sigma_CW**2 * sigma_0**2

range_sigma = 1 / m_sigma_CW * hbar_c_m

print(f"\ng_σ = y² m_χ / (16π² f) = {g_sigma:.2e}")
print(f"m_σ (CW) = {m_sigma_CW*1e9:.2e} eV")
print(f"Range of σ = {range_sigma:.1f} m")
print(f"\nσ₀ = g_σ × n_χ / m_σ² = {sigma_0:.2e} GeV")
print(f"ρ_σ = ½ m_σ² σ₀² = {rho_sigma:.2e} GeV⁴")
print(f"ρ_σ / ρ_Λ = {rho_sigma/rho_Lambda:.2e}")

# What m_σ would give ρ_Λ?
m_sigma_needed = np.sqrt(g_sigma**2 * n_chi**2 / (2 * rho_Lambda))
range_sigma_needed = 1 / m_sigma_needed * hbar_c_m

print(f"\n→ For ρ_σ = ρ_Λ, need m_σ = {m_sigma_needed:.2e} GeV = {m_sigma_needed*1e9:.2e} eV")
print(f"  Range = {range_sigma_needed:.2e} m = {range_sigma_needed/3.086e22:.2e} Mpc")

# ============= (3) Force perspective: F = ma ==============
print("\n" + "="*65)
print("  (3) F = ma — FORCE PERSPECTIVE")
print("="*65)

# The observed DE acceleration
a_DE = H_0**2 * rho_Lambda / rho_crit  # ~ H₀² Ω_Λ (in natural units, GeV)

# Force on a test DM particle from σ background gradient
# In a uniform medium: no gradient → no force
# At boundary of structure (halo edge): gradient exists
# Estimate: ∇σ ~ σ₀ / R_structure

R_halo = 100e3 * 3.086e22 * (1/hbar_c_m)  # 100 kpc in GeV⁻¹
# Hmm, let me use the Hubble radius for the cosmological gradient
R_H = 1 / H_0  # Hubble radius in GeV⁻¹

# From the σ equation, at cosmological scales, the "force" per unit mass:
# a_σ = (g_σ / m_χ) × ∇σ
# For a density gradient δρ over scale R:
# a_σ ~ g_σ² n_chi / (m_σ² m_χ R) × δρ/ρ .... complicated

# Simpler: the σ field contributes effective pressure
# P_σ = ½ σ̇² - V(σ). For slow-roll (static): P_σ = -V(σ) → w = -1!
# But V(σ) = ½ m_σ² σ₀² = g_σ² n_χ² / (2 m_σ²) which we already computed

# The key comparison: what is the ratio of φ-mediated force to gravity?
# At distance r between two χ particles:
# F_gravity / F_Yukawa = (G m_χ²) / (α_s / r²) × e^{m_φ r} for r >> 1/m_φ

# But at cosmological scales (r >> 1/m_φ), the Yukawa force is exponentially dead
# Gravity wins by e^{m_φ r} factor

r_cosmo = 1 / H_0 * hbar_c_m  # Hubble radius in meters
suppression = np.exp(-m_phi * (1/H_0))  # e^{-m_φ R_H}

print(f"\nHubble radius R_H = {r_cosmo:.2e} m")
print(f"m_φ × R_H = {m_phi/H_0:.2e}")
print(f"Yukawa suppression e^{{-m_φ R_H}} = {suppression}")
print(f"  (this is literally 10^{{-{m_phi/H_0/np.log(10):.0e}}})")

# ============= (4) What WORKS: ultra-light σ ==============
print("\n" + "="*65)
print("  (4) WHAT WOULD MAKE IT WORK — required σ mass")
print("="*65)

# Need: range of σ ≥ Hubble scale → m_σ ≲ H₀
# Need: ρ_σ = ρ_Λ → g² n² / (2 m_σ²) = ρ_Λ

# If m_σ = H₀:
rho_sigma_H0 = g_sigma**2 * n_chi**2 / (2 * H_0**2)
print(f"\nIf m_σ = H₀ = {H_0*1e9:.2e} eV:")
print(f"  ρ_σ = {rho_sigma_H0:.2e} GeV⁴")
print(f"  ρ_σ / ρ_Λ = {rho_sigma_H0/rho_Lambda:.2e}")

# Need: g_σ such that ρ_σ = ρ_Λ with m_σ = H₀
g_needed = np.sqrt(2 * rho_Lambda * H_0**2) / n_chi
print(f"\nFor ρ_σ = ρ_Λ with m_σ = H₀:")
print(f"  Need g_σ = {g_needed:.2e}")
print(f"  Have g_σ = {g_sigma:.2e}")
print(f"  Ratio (need/have) = {g_needed/g_sigma:.2e}")

# What f gives the right g_σ?
f_needed = y**2 * m_chi / (16 * np.pi**2 * g_needed)
print(f"\n  This requires f = {f_needed:.2e} GeV = {f_needed/M_Pl:.1f} M_Pl")

# ============= (5) CRITICAL INSIGHT — The real problem ==============
print("\n" + "="*65)
print("  (5) THE REAL PROBLEM")
print("="*65)
print(f"""
The accumulated φ/σ force CAN in principle give ρ_Λ.
But it requires:

  φ scenario: m_φ = {m_phi_needed*1e9:.1e} eV (range {range_needed/1e3:.0f} km)
    vs actual m_φ = {m_phi*1e9:.1e} eV (range {range_phi_m*1e15:.0f} fm)
    → OFF BY {np.log10(m_phi/m_phi_needed):.0f} ORDERS

  σ scenario: m_σ ≲ H₀ = {H_0*1e9:.1e} eV (range = Hubble)
    vs CW gives m_σ ~ 10⁻¹² eV (range = 30 m)
    → OFF BY {np.log10(m_sigma_CW/H_0):.0f} ORDERS

The force MUST have cosmological range to accumulate.
Yukawa (massive) = exponential death at distance >> 1/m.
Only MASSLESS or m ~ H₀ mediators can accumulate.

SIDM requirements: m_φ ~ MeV (for resonances)
DE requirements: m_mediator ~ H₀ (for accumulation)
Gap: {np.log10(m_phi/H_0):.0f} orders of magnitude.

→ Same mediator CANNOT do both SIDM and DE.
→ Need TWO mediators: φ (heavy, SIDM) + σ (ultralight, DE)
→ CW gives m_σ too heavy. Need protection mechanism.
""")

# ============= (6) IS THERE A WAY? ==============
print("="*65)
print("  (6) CAN σ BE PROTECTED?")
print("="*65)
print("""
In EM, m_photon = 0 is EXACT (gauge symmetry).
In our model, σ has shift symmetry σ → σ + c (approximate).
CW breaks it → m_σ ~ 10⁻¹² eV (too heavy).

KNOWN PROTECTION MECHANISMS:
1. Gauge symmetry (m = 0 exact)     — σ is scalar, not gauge boson ✗
2. Goldstone theorem (m = 0 from SSB) — σ IS a pseudo-Goldstone ✓
   BUT: explicit breaking (from y_s ≠ y_p) gives mass via CW
3. SUSY (cancels loops)             — if SUSY in dark sector:
   bosonic + fermionic loops cancel → m_σ protected
   Dark SUSY: χ (fermion) + φ (boson) could be superpartners!
4. Clockwork / relaxion mechanism    — exponentially small coupling
5. Sequestering (extra dimension)    — σ lives on a brane

Most promising: options 2+3. If dark sector has approximate SUSY
(m_χ ≈ m_φ), CW cancellation could suppress m_σ by (m_χ²-m_φ²)/f².
For MAP: m_χ = 94 MeV, m_φ = 1.94 MeV → partial cancellation only.
""")
