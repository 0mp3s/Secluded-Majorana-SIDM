#!/usr/bin/env python3
"""
freeze_out_trapping.py — Does DM freeze-out trap σ at θ_relic?
===============================================================

We solve the coupled system:
  1. σ equation of motion in FRW background
  2. χ number density (Boltzmann equation) with θ-dependent ⟨σv⟩
  3. Hubble expansion H(T)

The key coupling: ⟨σv⟩ depends on θ = σ/f through:
  ⟨σv⟩₀ = 2π α_s(θ) α_p(θ) / m²_χ = (y⁴ sin²(2θ)) / (128π m²_χ)

And the backreaction on σ:
  σ̈ + 3Hσ̇ + dV_eff/dσ = 0

where V_eff(σ) includes:
  - V_CW(θ): Coleman-Weinberg from χ loop
  - V_thermal(θ, T): finite-temperature correction ~ n_χ × ∂M_eff/∂σ
  - V_bare(σ): UV contribution (test with and without)

The question: starting from θ_initial ≠ θ_relic, does σ evolve
toward θ_relic = 19.47° during freeze-out?
"""

import numpy as np
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Constants (natural units: GeV)
# =============================================================================
GeV = 1.0
MeV = 1e-3
eV  = 1e-9
M_Pl = 2.435e18  # reduced Planck mass
g_star = 10.75    # relativistic dof at T ~ few MeV (below QCD, above e± annihilation)

# Conversion: 1 GeV⁻¹ = 6.58e-25 s
GeV_inv_to_sec = 6.58e-25  # seconds per GeV⁻¹

# Benchmark: MAP (best candidate for ρ_Λ match)
m_chi = 94.07 * MeV
m_phi = 11.10 * MeV
alpha = 5.734e-3

# Derived
theta_relic = np.arctan(1.0 / np.sqrt(8))  # 19.47°
y_sq = 4 * np.pi * alpha / np.cos(theta_relic)**2
y = np.sqrt(y_sq)

# Freeze-out temperature estimate: T_fo ~ m_χ / x_fo, x_fo ~ 20
x_fo = 20.0
T_fo = m_chi / x_fo

print("="*78)
print("FREEZE-OUT TRAPPING OF σ — MAP BENCHMARK")
print("="*78)
print(f"  m_χ = {m_chi/MeV:.2f} MeV")
print(f"  m_φ = {m_phi/MeV:.2f} MeV")
print(f"  α = {alpha:.4e}")
print(f"  y = {y:.5f}")
print(f"  θ_relic = {np.degrees(theta_relic):.2f}°")
print(f"  T_fo ≈ m_χ/{x_fo:.0f} = {T_fo/MeV:.2f} MeV")
print()

# =============================================================================
# Physics functions
# =============================================================================

def Hubble(T):
    """Hubble rate in radiation domination: H = √(π²g*/90) T²/M_Pl"""
    return np.sqrt(np.pi**2 * g_star / 90) * T**2 / M_Pl

def n_eq(T, m):
    """Equilibrium number density for non-relativistic species"""
    if T <= 0 or m/T > 500:
        return 0.0
    return 2 * (m * T / (2 * np.pi))**1.5 * np.exp(-m / T)  # Majorana: g=2

def sigma_v_s_wave(theta):
    """s-wave annihilation cross section × velocity, θ-dependent"""
    # ⟨σv⟩₀ = 2π α_s α_p / m²_χ
    alpha_s = (y**2 / (4*np.pi)) * np.cos(theta)**2
    alpha_p = (y**2 / (4*np.pi)) * np.sin(theta)**2
    return 2 * np.pi * alpha_s * alpha_p / m_chi**2

def dVCW_dtheta(theta, v_phi):
    """
    Derivative of CW potential w.r.t. θ = σ/f
    
    V_CW(θ) = -(1/32π²) M⁴_eff [ln(M²_eff/m²_χ) - 3/2]
    M²_eff = m²_χ + m_χ y cosθ v + y² v²/4
    
    dV/dθ = -(1/32π²) × d/dθ { M⁴ [ln(M²/μ²) - 3/2] }
           = -(1/32π²) × dM²/dθ × { 2M² [ln(M²/μ²) - 3/2] + M² × 2 }
           = -(1/16π²) × dM²/dθ × M² × [ln(M²/μ²) - 1]
    
    dM²/dθ = -m_χ y sinθ v
    """
    v = v_phi
    M2 = m_chi**2 + m_chi * y * np.cos(theta) * v + y**2 * v**2 / 4
    if M2 <= 0:
        return 0.0
    dM2_dth = -m_chi * y * np.sin(theta) * v
    return -(1/(16*np.pi**2)) * dM2_dth * M2 * (np.log(M2/m_chi**2) - 1.0)

def thermal_force(theta, n_chi, v_phi):
    """
    Backreaction force from χ number density on σ.
    
    The finite-density effective potential includes:
    δV_density ≈ n_χ × M_eff(θ)  (leading non-relativistic correction)
    
    dδV/dθ = n_χ × dM_eff/dθ
    
    M_eff(θ) = √(m²_χ + m_χ y cosθ v + y² v²/4)
    dM_eff/dθ = -m_χ y sinθ v / (2 M_eff)
    """
    v = v_phi
    M2 = m_chi**2 + m_chi * y * np.cos(theta) * v + y**2 * v**2 / 4
    if M2 <= 0:
        return 0.0
    M = np.sqrt(M2)
    dM_dth = -m_chi * y * np.sin(theta) * v / (2 * M)
    return n_chi * dM_dth

# =============================================================================
# Coupled ODE system
# =============================================================================
# Variables: use x = m_χ/T as time variable (increases as T drops)
#
# Boltzmann: dY/dx = -(s⟨σv⟩/H x) (Y² - Y²_eq)
#   where Y = n_χ/s, s = (2π²/45) g_star T³
#
# σ equation (in terms of θ = σ/f):
#   θ̈ + 3H θ̇ + (1/f²)[dV_CW/dθ + n_χ dM/dθ] = 0
#
# Rewrite using x = m_χ/T:
#   dt = -dx/(xH)  (since T = m_χ/x, dT/dt = -HT → dx/dt = xH)
#
#   dθ/dx = θ̇/(xH) ≡ ω/(xH)
#   dω/dx = [-3Hω - (1/f²)(dV_CW/dθ + n_χ dM/dθ)] / (xH)
#   dY/dx = -(s⟨σv⟩(θ))/(Hx) × (Y² - Y²_eq(x))

def entropy_density(T):
    """s = (2π²/45) g_star T³"""
    return (2 * np.pi**2 / 45) * g_star * T**3

def system(x, state, f_decay):
    """
    state = [θ, ω, Y]
    θ = σ/f
    ω = dθ/dt (proper time derivative)
    Y = n_χ / s
    """
    theta, omega, Y = state
    
    # Keep θ in [0, π/2] — physical range
    theta = np.clip(theta, 1e-6, np.pi/2 - 1e-6)
    
    T = m_chi / x
    H = Hubble(T)
    s = entropy_density(T)
    
    # Equilibrium
    n_eq_val = n_eq(T, m_chi)
    Y_eq = n_eq_val / s if s > 0 else 0.0
    
    # Current number density
    n_chi = Y * s
    
    # Annihilation cross section at current θ
    sv = sigma_v_s_wave(theta)
    
    # φ VEV (estimate: small, ~ fraction of m_phi from cannibal)
    v_phi = 0.5 * m_phi  # conservative estimate
    
    # Forces on σ
    F_CW = dVCW_dtheta(theta, v_phi)
    F_thermal = thermal_force(theta, n_chi, v_phi)
    
    # Total force per f²
    total_force = (F_CW + F_thermal) / f_decay**2
    
    # dθ/dx = ω / (xH)
    dtheta_dx = omega / (x * H)
    
    # dω/dx = [-3Hω - total_force] / (xH)
    domega_dx = (-3 * H * omega - total_force) / (x * H)
    
    # dY/dx = -(s ⟨σv⟩)/(H x) (Y² - Y²_eq)
    dY_dx = -(s * sv) / (H * x) * (Y**2 - Y_eq**2)
    
    return [dtheta_dx, domega_dx, dY_dx]


# =============================================================================
# Run simulations
# =============================================================================
print("="*78)
print("SIMULATION: σ EVOLUTION DURING FREEZE-OUT")
print("="*78)
print()

# x range: from x=1 (T=m_χ, hot) to x=200 (T=m_χ/200, frozen)
x_span = (1.0, 200.0)
x_eval = np.clip(np.logspace(0, np.log10(200), 1000), 1.0, 200.0)

# Initial conditions: Y starts at equilibrium
T_init = m_chi / x_span[0]
s_init = entropy_density(T_init)
Y_init = n_eq(T_init, m_chi) / s_init

# Test multiple initial θ values
theta_initials = [0.05, 0.15, 0.34, 0.50, 0.80, 1.2]  # radians
# θ_relic ≈ 0.3398 rad

# Test multiple f values
f_values = [
    (0.2 * M_Pl, "0.2 M_Pl"),
    (1.0 * M_Pl, "1.0 M_Pl"),
    (15.0 * M_Pl, "15 M_Pl"),
]

for f_decay, f_label in f_values:
    print(f"\n{'='*78}")
    print(f"f = {f_label}  (β = M_Pl/f = {M_Pl/f_decay:.2f})")
    print(f"{'='*78}")
    
    # Timescales
    H_fo = Hubble(T_fo)
    t_H_fo = 1.0 / H_fo  # Hubble time at freeze-out
    
    # σ mass from CW: m²_σ = (1/f²) d²V_CW/dθ²
    # Numerical: compute d²V/dθ² at the CW minimum (θ ≈ π/2)
    v_phi = 0.5 * m_phi
    def V_CW_func(theta):
        M2 = m_chi**2 + m_chi * y * np.cos(theta) * v_phi + y**2 * v_phi**2 / 4
        return -(1/(32*np.pi**2)) * M2**2 * (np.log(M2/m_chi**2) - 1.5) if M2 > 0 else 1e100
    dth = 1e-6
    th_min = np.pi/2  # CW minimum is at θ=π/2 (see verify_cw_bugs.py)
    d2V_dth2 = (V_CW_func(th_min + dth) - 2*V_CW_func(th_min) + V_CW_func(th_min - dth)) / dth**2
    m2_sigma = d2V_dth2 / f_decay**2
    m_sigma = np.sqrt(abs(m2_sigma))
    m_sigma_sign = '+' if d2V_dth2 > 0 else '-' # sign tells if minimum or maximum
    
    print(f"  H(T_fo) = {H_fo:.3e} GeV")
    print(f"  m_σ (CW at θ=π/2) = {m_sigma:.3e} GeV  (m²_σ sign: {m_sigma_sign})")
    print(f"  |m_σ| / H(T_fo) = {m_sigma/H_fo:.2e}")
    
    if d2V_dth2 < 0:
        print(f"  ⚠ m²_σ < 0 at θ=π/2: CW minimum is actually a MAXIMUM here!")
        # Compute dV/dθ magnitude to estimate rolling timescale
        F_CW_30 = abs(dVCW_dtheta(np.radians(30), v_phi))
        accel = F_CW_30 / f_decay**2
        t_roll = np.sqrt(2 * theta_relic / accel) if accel > 0 else np.inf
        print(f"    Rolling time from θ_relic: {t_roll:.3e} GeV⁻¹ vs 1/H = {1/H_fo:.3e} GeV⁻¹")
        print(f"    Rolls in {t_roll * H_fo:.2e} Hubble times")
    elif m_sigma > H_fo:
        print(f"  → m_σ > H: σ oscillates → settles to CW minimum (θ=π/2)")
        print(f"    Oscillation period / Hubble time = {H_fo/m_sigma:.2e}")
    else:
        print(f"  → m_σ < H: σ is FROZEN by Hubble friction → stays at initial θ")
    print()
    
    print(f"  {'θ_init (°)':<12} {'θ_final (°)':<14} {'→ θ_relic?':<12} "
          f"{'Y_final':<12} {'Y(θ_relic)':<12} {'ΔΩ/Ω':<10}")
    print(f"  {'-'*72}")
    
    # Reference: Y at θ_relic
    Y_ref = None
    
    results = []
    for theta_i in theta_initials:
        state0 = [theta_i, 0.0, Y_init]  # start at rest
        
        try:
            sol = solve_ivp(
                lambda x, s: system(x, s, f_decay),
                x_span, state0, 
                t_eval=x_eval,
                method='RK45',
                rtol=1e-8, atol=1e-12,
                max_step=0.5
            )
            
            if sol.success:
                theta_f = sol.y[0, -1]
                Y_f = sol.y[2, -1]
                
                if abs(theta_i - theta_relic) < 0.01:
                    Y_ref = Y_f
                
                results.append((theta_i, theta_f, Y_f))
                
                at_relic = "YES" if abs(theta_f - theta_relic) < 0.05 else "no"
                print(f"  {np.degrees(theta_i):<12.1f} {np.degrees(theta_f):<14.2f} "
                      f"{at_relic:<12} {Y_f:<12.4e}", end="")
                
                if Y_ref is not None and Y_ref > 0:
                    delta_omega = abs(Y_f - Y_ref) / Y_ref
                    print(f" {Y_ref:<12.4e} {delta_omega:<10.3e}")
                else:
                    print()
            else:
                print(f"  {np.degrees(theta_i):<12.1f} FAILED: {sol.message}")
        except Exception as e:
            print(f"  {np.degrees(theta_i):<12.1f} ERROR: {e}")
    
    # Track evolution of θ for the θ_init = 0.05 case
    if results:
        print()
        print(f"  Evolution snapshots for θ_init = {np.degrees(theta_initials[0]):.1f}°:")
        
        sol = solve_ivp(
            lambda x, s: system(x, s, f_decay),
            x_span, [theta_initials[0], 0.0, Y_init],
            t_eval=x_eval, method='RK45',
            rtol=1e-8, atol=1e-12, max_step=0.5
        )
        
        if sol.success:
            checkpoints = [1, 5, 10, 20, 50, 100, 200]
            print(f"    {'x=m/T':<8} {'T (MeV)':<10} {'θ (°)':<10} {'Y':<12} {'⟨σv⟩/⟨σv⟩_relic':<16}")
            for xc in checkpoints:
                idx = np.argmin(np.abs(sol.t - xc))
                T_c = m_chi / sol.t[idx]
                theta_c = sol.y[0, idx]
                Y_c = sol.y[2, idx]
                sv_c = sigma_v_s_wave(theta_c)
                sv_relic = sigma_v_s_wave(theta_relic)
                print(f"    {sol.t[idx]:<8.1f} {T_c/MeV:<10.3f} "
                      f"{np.degrees(theta_c):<10.2f} {Y_c:<12.4e} "
                      f"{sv_c/sv_relic:<16.4f}")

# =============================================================================
# PART 2: The backreaction strength
# =============================================================================
print()
print("="*78)
print("PART 2: BACKREACTION ANALYSIS")
print("="*78)
print()
print("Is the thermal force strong enough to move σ?")
print("Note: dV_CW/dθ < 0 at θ_relic → CW pushes σ TOWARD larger θ (toward π/2)")
print()

T = T_fo
H = Hubble(T)
s = entropy_density(T)
n_chi_fo = n_eq(T, m_chi)  # at freeze-out, n ~ n_eq
v_phi = 0.5 * m_phi

for f_decay, f_label in f_values:
    print(f"f = {f_label}:")
    
    # CW force
    F_CW = abs(dVCW_dtheta(theta_relic, v_phi))
    
    # Thermal force at freeze-out
    F_th = abs(thermal_force(theta_relic, n_chi_fo, v_phi))
    
    # σ acceleration: θ̈ = F/(f²)
    # Displacement in one Hubble time: Δθ ~ θ̈ / (2H²) ~ F/(f² 2H²)
    delta_theta_CW = F_CW / (f_decay**2 * 2 * H**2)
    delta_theta_th = F_th / (f_decay**2 * 2 * H**2)
    
    print(f"  F_CW at θ_relic = {F_CW:.3e} GeV⁴")
    print(f"  F_thermal at T_fo = {F_th:.3e} GeV⁴")
    print(f"  Ratio F_th/F_CW = {F_th/F_CW:.3e}")
    print(f"  Δθ from CW in 1/H = {delta_theta_CW:.3e} rad = {np.degrees(delta_theta_CW):.3e}°")
    print(f"  Δθ from thermal in 1/H = {delta_theta_th:.3e} rad = {np.degrees(delta_theta_th):.3e}°")
    print(f"  θ_relic = {theta_relic:.4f} rad = {np.degrees(theta_relic):.2f}°")
    
    if delta_theta_CW > 0.1:
        print(f"  → CW MOVES σ significantly in one Hubble time!")
    elif delta_theta_CW > 1e-6:
        print(f"  → CW moves σ slowly — many Hubble times to settle")
    else:
        print(f"  → CW force NEGLIGIBLE — σ is frozen by Hubble friction")
    print()

# =============================================================================
# PART 3: The adiabatic regime — m_σ vs H(T)
# =============================================================================
print("="*78)
print("PART 3: ADIABATIC REGIME — WHEN DOES σ START OSCILLATING?")
print("="*78)
print()
print("σ oscillates when m_σ > H(T).")
print("In radiation domination: H ∝ T² → H decreases as universe cools.")
print("m_σ from CW is roughly constant (set by m_χ, y, v_phi, f).")
print()

# Use numerical d²V/dθ² at the correct CW minimum (θ=π/2)
def V_CW_theta(theta):
    M2 = m_chi**2 + m_chi * y * np.cos(theta) * v_phi + y**2 * v_phi**2 / 4
    return -(1/(32*np.pi**2)) * M2**2 * (np.log(M2/m_chi**2) - 1.5) if M2 > 0 else 1e100

dth = 1e-6
print("Note: CW minimum is at θ = π/2 (verified numerically), NOT θ=0.")
print("V_CW decreases monotonically from θ=0 to θ=π/2.")
print("So CW pushes σ AWAY from θ_relic toward θ=π/2 (pure pseudoscalar).")
print()

summary_rows = []

for f_decay, f_label in f_values:
    # CW curvature at θ=π/2 (minimum)
    th_min = np.pi/2
    d2V = (V_CW_theta(th_min + dth) - 2*V_CW_theta(th_min) + V_CW_theta(th_min - dth)) / dth**2
    m2_sig = d2V / f_decay**2
    m_sig = np.sqrt(abs(m2_sig))
    sign_str = '+' if d2V > 0 else '-'
    
    # CW force at θ_relic (this drives the rolling speed)
    F_at_relic = abs(dVCW_dtheta(theta_relic, v_phi))
    accel_relic = F_at_relic / f_decay**2  # angular acceleration
    # Rolling timescale from θ_relic to π/2: rough estimate Δθ ∼ a t²/2
    delta_th = np.pi/2 - theta_relic
    t_roll = np.sqrt(2 * delta_th / accel_relic) if accel_relic > 0 else np.inf
    H_fo_val = Hubble(T_fo)
    
    print(f"  f = {f_label}:")
    print(f"    d²V/dθ² at θ=π/2: {d2V:.3e} GeV⁴ (sign: {sign_str})")
    print(f"    |m_σ| = {m_sig:.3e} GeV = {m_sig/eV:.2e} eV")
    print(f"    |m_σ|/H(T_fo) = {m_sig/H_fo_val:.2e}")
    print(f"    CW force at θ_relic: {F_at_relic:.3e} GeV⁴")
    print(f"    Rolling time to π/2: {t_roll:.3e} GeV⁻¹ vs 1/H = {1/H_fo_val:.3e} GeV⁻¹")
    nH_roll = t_roll * H_fo_val
    print(f"    → Rolls in {nH_roll:.2e} Hubble times")
    
    if d2V < 0:
        behavior = "CW min is tachyonic!"
        theta_final = "θ→π/2"
    elif nH_roll < 1:
        behavior = f"Fast roll to θ=π/2 in {nH_roll:.1e} t_H"
        theta_final = "θ→π/2 ✗"
    elif nH_roll < 100:
        behavior = f"Slow roll to π/2 ({nH_roll:.0f} t_H)"
        theta_final = "θ→π/2 ✗"
    else:
        behavior = "Frozen by Hubble friction"
        theta_final = "θ stays"
    
    print(f"    → {behavior}")
    summary_rows.append((f_label, f"{m_sig/H_fo_val:.2e}", behavior, theta_final))
    print()

# =============================================================================
# SUMMARY
# =============================================================================
print("="*78)
print("SUMMARY")
print("="*78)
print()
print("CRITICAL CORRECTION: CW minimum is at θ=π/2, NOT θ=0!")
print("  V_CW ∝ cosθ (leading θ-dependent term)")
print("  θ=0 is CW MAXIMUM (largest M_eff), θ=π/2 is MINIMUM (smallest ΔM)")
print("  CW pushes σ AWAY from θ_relic toward θ=π/2 (α_s=0, no SIDM)")
print()
print("The fate of σ depends on f:")
print()
print(f"  {'f':<12}| {'|mσ|/H':<12}| {'Behavior':<40}| {'Fate'}")
print(f"  {'-'*12}|{'-'*13}|{'-'*41}|{'-'*15}")
for f_label, mH_ratio, behavior, theta_final in summary_rows:
    print(f"  {f_label:<12}| {mH_ratio:<12}| {behavior:<40}| {theta_final}")
print()
print("CONCLUSION:")
print("  1. CW drives σ toward θ=π/2 (pure pseudoscalar, α_s=0).")
print("  2. Thermal backreaction (F_th/F_CW ~ 10⁻¹⁰) cannot resist this.")
print("  3. At θ=π/2: no scalar coupling → no SIDM. Phenomenologically dead.")
print()
print("  Freeze-out CANNOT trap σ at θ_relic. The CW potential is the enemy.")
print()
print("  Options:")
print("  1. V_bare(σ) with minimum at θ_relic that OVERCOMES V_CW")
print("  2. Non-minimal kinetic term K(σ)(∂σ)² providing extra friction")
print("  3. σ is not a dynamical field — the angle is fixed by UV completion")
print("  4. Additional particles in the dark sector modify V_CW topology")
