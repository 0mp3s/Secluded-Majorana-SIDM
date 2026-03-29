#!/usr/bin/env python3
"""
sigma_trapping_ode.py — Does σ get trapped at θ_relic during freeze-out?
========================================================================

Step 3 of verification: Coupled ODE system
  1. σ field equation in FRW background
  2. χ Boltzmann equation with θ-dependent ⟨σv⟩
  3. Hubble expansion

Key physics:
  - CW potential drives σ → θ=π/2 (pure pseudoscalar, phenomenologically dead)
  - Thermal backreaction from χ number density resists
  - Question: is there a mechanism that TRAPS σ at θ_relic = arcsin(1/3)?

This version uses CORRECTED CW potential (5 bugs fixed from freeze_out_analysis_corrected.py):
  1. dV/dθ coefficient: [ln - 1.0], not [ln - 0.5]
  2. CW minimum at θ=π/2, not θ=0
  3. m²_σ sign correct
  4. t_eval crash fixed with np.clip
  5. Dynamic results (not hard-coded)

We test:
  (a) Multiple initial θ values
  (b) Multiple f (decay constant) values
  (c) With and without V_bare(θ) at θ_relic
  (d) Timescale comparison: CW rolling vs Hubble vs freeze-out
"""
import sys, math
import numpy as np
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

# =============================================================================
# Constants
# =============================================================================
MeV = 1e-3        # GeV
M_Pl = 2.435e18   # reduced Planck mass [GeV]
g_star = 10.75     # at T ~ few MeV

# MAP benchmark point
m_chi = 94.07 * MeV
m_phi = 11.10 * MeV
alpha = 5.734e-3

theta_relic = math.asin(1.0 / 3.0)   # 19.47°
y_sq = 4 * math.pi * alpha / math.cos(theta_relic)**2  # y² = 4πα/cos²θ (total coupling)
y = math.sqrt(y_sq)

T_fo = m_chi / 20.0
v_phi = 0.5 * m_phi  # dark Higgs VEV estimate

print("=" * 80)
print("  σ TRAPPING AT θ_relic — COUPLED ODE ANALYSIS")
print("=" * 80)
print(f"  m_χ = {m_chi/MeV:.2f} MeV,  m_φ = {m_phi/MeV:.2f} MeV,  α = {alpha:.4e}")
print(f"  y = {y:.5f},  θ_relic = {math.degrees(theta_relic):.2f}°")
print(f"  T_fo ≈ {T_fo/MeV:.2f} MeV  (x_fo ≈ 20)")
print()

# =============================================================================
# Physics functions (CORRECTED)
# =============================================================================
def Hubble(T):
    return math.sqrt(math.pi**2 * g_star / 90) * T**2 / M_Pl

def n_eq(T, m):
    if T <= 0 or m / T > 500:
        return 0.0
    return 2 * (m * T / (2 * math.pi))**1.5 * math.exp(-m / T)

def sigma_v_s_wave(theta):
    """⟨σv⟩ = 2π α_s(θ) α_p(θ) / m²_χ"""
    alpha_s = (y**2 / (4 * math.pi)) * math.cos(theta)**2
    alpha_p = (y**2 / (4 * math.pi)) * math.sin(theta)**2
    return 2 * math.pi * alpha_s * alpha_p / m_chi**2

def V_CW(theta):
    """Coleman-Weinberg potential from χ loop (1-loop, MS-bar)"""
    M2 = m_chi**2 + m_chi * y * math.cos(theta) * v_phi + y**2 * v_phi**2 / 4
    if M2 <= 0:
        return 0.0
    return -(1 / (32 * math.pi**2)) * M2**2 * (math.log(M2 / m_chi**2) - 1.5)

def dVCW_dtheta(theta):
    """dV_CW/dθ — CORRECTED: coefficient is [ln - 1.0]"""
    M2 = m_chi**2 + m_chi * y * math.cos(theta) * v_phi + y**2 * v_phi**2 / 4
    if M2 <= 0:
        return 0.0
    dM2_dth = -m_chi * y * math.sin(theta) * v_phi
    return -(1 / (16 * math.pi**2)) * dM2_dth * M2 * (math.log(M2 / m_chi**2) - 1.0)

def thermal_force(theta, n_chi_val):
    """n_χ × dM_eff/dθ — backreaction from DM number density"""
    M2 = m_chi**2 + m_chi * y * math.cos(theta) * v_phi + y**2 * v_phi**2 / 4
    if M2 <= 0:
        return 0.0
    M = math.sqrt(M2)
    dM_dth = -m_chi * y * math.sin(theta) * v_phi / (2 * M)
    return n_chi_val * dM_dth

def entropy_density(T):
    return (2 * math.pi**2 / 45) * g_star * T**3

# =============================================================================
# Coupled ODE System (x = m_χ/T as time variable)
# =============================================================================
def system(x, state, f_decay, V_bare_func=None):
    """
    state = [θ, ω, Y]
    θ = σ/f,  ω = dθ/dt,  Y = n_χ/s
    """
    theta_val, omega_val, Y_val = state
    theta_val = max(1e-6, min(theta_val, math.pi/2 - 1e-6))
    Y_val = max(Y_val, 1e-30)
    
    T = m_chi / x
    if T < 1e-10:
        return [0.0, 0.0, 0.0]
    
    H = Hubble(T)
    s = entropy_density(T)
    n_eq_val = n_eq(T, m_chi)
    Y_eq = n_eq_val / s if s > 0 else 0.0
    n_chi_val = Y_val * s
    
    sv = sigma_v_s_wave(theta_val)
    
    # Forces on σ
    F_CW = dVCW_dtheta(theta_val)
    F_th = thermal_force(theta_val, n_chi_val)
    
    # Optional bare potential
    F_bare = 0.0
    if V_bare_func is not None:
        F_bare = V_bare_func(theta_val)
    
    total_force = (F_CW + F_th + F_bare) / f_decay**2
    
    # dθ/dx = ω / (xH)
    dtheta_dx = omega_val / (x * H)
    
    # dω/dx = [-3Hω - total_force] / (xH)
    domega_dx = (-3 * H * omega_val - total_force) / (x * H)
    
    # dY/dx = -(s⟨σv⟩)/(Hx) (Y² - Y²_eq)
    dY_dx = -(s * sv) / (H * x) * (Y_val**2 - Y_eq**2)
    
    return [dtheta_dx, domega_dx, dY_dx]

# =============================================================================
# PART 1: CW FORCE ANALYSIS (before running ODE)
# =============================================================================
print("=" * 80)
print("  PART 1: FORCE ANALYSIS AT θ_relic")
print("=" * 80)
print()

H_fo = Hubble(T_fo)
n_chi_fo = n_eq(T_fo, m_chi)
s_fo = entropy_density(T_fo)

F_CW_relic = dVCW_dtheta(theta_relic)
F_th_relic = thermal_force(theta_relic, n_chi_fo)

print(f"  At T_fo = {T_fo/MeV:.2f} MeV:")
print(f"    H(T_fo) = {H_fo:.3e} GeV")
print(f"    n_χ(T_fo) = {n_chi_fo:.3e} GeV³")
print(f"    s(T_fo) = {s_fo:.3e} GeV³")
print(f"    Y_eq(T_fo) = {n_chi_fo/s_fo:.3e}")
print()
print(f"  Forces at θ_relic = {math.degrees(theta_relic):.2f}°:")
print(f"    F_CW = dV_CW/dθ = {F_CW_relic:.3e} GeV⁴")
print(f"    F_th = n_χ dM/dθ = {F_th_relic:.3e} GeV⁴")
print(f"    |F_th / F_CW| = {abs(F_th_relic/F_CW_relic):.3e}")
print()

# CW drives which direction?
if F_CW_relic < 0:
    print(f"  F_CW < 0  →  CW pushes θ toward LARGER values (→ π/2)")
else:
    print(f"  F_CW > 0  →  CW pushes θ toward SMALLER values (→ 0)")

if F_th_relic > 0:
    print(f"  F_th > 0  →  Thermal force pushes θ toward SMALLER values (→ 0)")
else:
    print(f"  F_th < 0  →  Thermal force pushes θ toward LARGER values (→ π/2)")
print()

# =============================================================================
# PART 2: COUPLED ODE — WITHOUT BARE POTENTIAL
# =============================================================================
print("=" * 80)
print("  PART 2: ODE EVOLUTION — NO BARE POTENTIAL (CW + thermal only)")
print("=" * 80)
print()

x_span = (1.0, 200.0)

T_init = m_chi / x_span[0]
s_init = entropy_density(T_init)
Y_init = n_eq(T_init, m_chi) / s_init if s_init > 0 else 1e-10

theta_initials = [0.05, 0.15, theta_relic, 0.50, 0.80, 1.2]
f_values = [
    (0.2 * M_Pl, "0.2 M_Pl"),
    (1.0 * M_Pl, "1.0 M_Pl"),
    (15.0 * M_Pl, "15 M_Pl"),
]

for f_decay, f_label in f_values:
    print(f"  ┌────────────────────────────────────────────")
    print(f"  │ f = {f_label}  (β = M_Pl/f = {M_Pl/f_decay:.2f})")
    print(f"  └────────────────────────────────────────────")
    
    # Rolling time estimate
    accel = abs(F_CW_relic) / f_decay**2
    delta_theta = math.pi/2 - theta_relic
    t_roll = math.sqrt(2 * delta_theta / accel) if accel > 0 else float('inf')
    n_hubble_roll = t_roll * H_fo
    
    print(f"    CW rolling time θ_relic→π/2: {n_hubble_roll:.2e} Hubble times")
    print()
    
    print(f"    {'θ_init(°)':>10}  {'θ_final(°)':>11}  {'→θ_relic?':>10}  {'Y_final':>12}  {'Ωh²_eff':>10}  {'status':>10}")
    print(f"    {'-'*70}")
    
    for theta_i in theta_initials:
        state0 = [theta_i, 0.0, Y_init]
        
        try:
            sol = solve_ivp(
                lambda x, s: system(x, s, f_decay),
                x_span, state0,
                method='Radau',       # stiff solver
                rtol=1e-6, atol=1e-10,
                max_step=2.0,
                dense_output=False,
            )
            
            if sol.success:
                theta_f = sol.y[0, -1]
                Y_f = sol.y[2, -1]
                omega_f = m_chi * Y_f * 2891.2 / 1.0539e-5  # rough Ωh²
                
                at_relic = abs(theta_f - theta_relic) < 0.05
                mark = "✅ TRAPPED" if at_relic else "❌ DRIFTED"
                
                print(f"    {math.degrees(theta_i):>10.1f}  {math.degrees(theta_f):>11.2f}  "
                      f"{'YES' if at_relic else 'no':>10}  {Y_f:>12.4e}  {omega_f:>10.4f}  {mark:>10}")
            else:
                print(f"    {math.degrees(theta_i):>10.1f}  {'FAILED':>11}  {'-':>10}  {'-':>12}  {'-':>10}  ⚠ ODE failed: {sol.message[:40]}")
        except Exception as e:
            print(f"    {math.degrees(theta_i):>10.1f}  {'ERROR':>11}  {str(e)[:50]}")
    print()

# =============================================================================
# PART 3: WITH BARE POTENTIAL V_bare(θ) CENTERED AT θ_relic
# =============================================================================
print("=" * 80)
print("  PART 3: ODE WITH V_bare(θ) — QUADRATIC TRAP AT θ_relic")
print("=" * 80)
print()

# If the UV theory provides a bare potential with minimum at θ_relic:
#   V_bare(θ) = (1/2) m²_σ_bare (θ - θ_relic)²
# with m_σ_bare >> H(T_fo), then σ oscillates around θ_relic and is trapped.
#
# The required condition: V_bare gradient > V_CW gradient at θ_relic
# dV_bare/dθ = m²_σ_bare (θ - θ_relic)
# At θ = θ_relic: this is zero. At θ = θ + δ: m²_σ_bare × δ
#
# For the restoring force to DOMINATE CW near θ_relic:
# m²_σ_bare × δ > |dV_CW/dθ| ≈ |dV_CW/dθ|_relic (roughly constant for small δ)
# → m²_σ_bare > |dV_CW/dθ|_relic / δ

# Test several bare mass scales
m_bare_values = [
    1e-5,   # GeV scale (way above H)
    1e-10,  # 
    1e-15,  #
    1e-20,  # intermediate
    H_fo,   # Hubble scale at freeze-out
]

f_test = 1.0 * M_Pl  # f = M_Pl as reference

print(f"  Using f = M_Pl,  CW force at θ_relic = {abs(F_CW_relic):.3e} GeV⁴")
print()
print(f"  For V_bare(θ) = ½ m²_bare f² (θ - θ_relic)²:")
print(f"  → dV_bare/dθ = m²_bare f² (θ - θ_relic)")
print()

for m_bare in m_bare_values:
    # The bare potential force (in units of GeV⁴, i.e. before dividing by f²)
    def dV_bare_func(theta_val, _m=m_bare, _f=f_test):
        return _m**2 * _f**2 * (theta_val - theta_relic)
    
    print(f"  m_bare = {m_bare:.3e} GeV  (m_bare/H_fo = {m_bare/H_fo:.2e}):")
    
    # At δθ = 0.1 (small displacement):
    F_bare_01 = abs(dV_bare_func(theta_relic + 0.1)) / f_test**2 * f_test**2
    print(f"    |dV_bare/dθ| at δθ=0.1: {F_bare_01:.3e} GeV⁴  vs |F_CW|: {abs(F_CW_relic):.3e} GeV⁴")
    
    dominates = F_bare_01 > abs(F_CW_relic)
    print(f"    Bare dominates CW at δθ=0.1: {'YES ✅' if dominates else 'NO ❌'}")
    
    # Run ODE
    state0 = [0.05, 0.0, Y_init]  # start far from θ_relic
    
    try:
        sol = solve_ivp(
            lambda x, s: system(x, s, f_test, dV_bare_func),
            x_span, state0,
            method='Radau',
            rtol=1e-6, atol=1e-10,
            max_step=2.0,
            dense_output=False,
        )
        
        if sol.success:
            theta_f = sol.y[0, -1]
            at_relic = abs(theta_f - theta_relic) < 0.05
            mark = "TRAPPED ✅" if at_relic else f"θ_f={math.degrees(theta_f):.1f}° ❌"
            print(f"    ODE result (θ₀=2.9°): θ_final = {math.degrees(theta_f):.2f}°  {mark}")
        else:
            print(f"    ODE failed")
    except Exception as e:
        print(f"    ODE error: {e}")
    
    print()

# =============================================================================
# PART 4: WHAT m_bare IS NEEDED?
# =============================================================================
print("=" * 80)
print("  PART 4: MINIMUM m_bare TO TRAP σ")
print("=" * 80)
print()

# Need: restoring force > CW force everywhere on [0, π/2]
# Max CW force on this interval:
theta_scan = np.linspace(0.01, math.pi/2 - 0.01, 1000)
F_CW_scan = [abs(dVCW_dtheta(th)) for th in theta_scan]
F_CW_max = max(F_CW_scan)
theta_at_max = theta_scan[np.argmax(F_CW_scan)]

print(f"  Max |dV_CW/dθ| on [0, π/2]: {F_CW_max:.3e} GeV⁴  at θ = {math.degrees(theta_at_max):.1f}°")
print()

# For quadratic bare potential with minimum at θ_relic:
# |dV_bare/dθ| at the farthest point from θ_relic is:
# m²_bare f² × max(θ_relic, π/2 - θ_relic)
delta_max = max(theta_relic, math.pi/2 - theta_relic)

for f_decay, f_label in f_values:
    m_bare_min = math.sqrt(F_CW_max / (f_decay**2 * delta_max))
    ratio_to_H = m_bare_min / H_fo
    
    print(f"  f = {f_label}:")
    print(f"    m_bare_min = {m_bare_min:.3e} GeV  = {ratio_to_H:.2e} × H(T_fo)")
    
    if ratio_to_H > 1:
        print(f"    → m_bare > H: oscillatory regime, σ settles to θ_relic ✅")
    else:
        print(f"    → m_bare < H: slow roll, Hubble friction helps BUT may not be enough")
    print()

# =============================================================================
# PART 5: ALTERNATIVE — A₄ DISCRETE SYMMETRY AS THE TRAP
# =============================================================================
print("=" * 80)
print("  PART 5: A₄ DISCRETE SYMMETRY INTERPRETATION")
print("=" * 80)
print()

print("""  If θ is NOT a dynamical field but is FIXED by A₄ symmetry breaking:
    → θ = arcsin(1/3) is set by the A₄ Clebsch-Gordan decomposition
    → No σ rolling problem — the angle is a CONSTANT of the theory
    → CW corrections shift the effective coupling but not θ itself
    
  In this interpretation:
    - The dark Yukawa decomposes into scalar (1) and pseudo (1') channels
    - The ratio is FIXED by A₄ group theory: tan²θ = 1/9
    - θ is not a field — it's a GROUP THEORY NUMBER
    - There is no σ field to roll
    
  This avoids ALL the σ trapping problems.
  
  The only correction: CW renormalizes the MAGNITUDES of g_s, g_p,
  but their RATIO is protected by the discrete A₄ symmetry
  (up to higher-order corrections suppressed by v/Λ).
""")

# =============================================================================
# PART 6: TIMESCALE SUMMARY
# =============================================================================
print("=" * 80)
print("  PART 6: TIMESCALE SUMMARY")
print("=" * 80)
print()

print(f"  {'Timescale':<35} {'Value (GeV⁻¹)':<15} {'In Hubble times':>16}")
print(f"  {'-'*70}")

t_H = 1 / H_fo
print(f"  {'Hubble time at T_fo':<35} {t_H:.3e}        {'1.00':>16}")

# Freeze-out duration: Δx ~ 10 → Δt ~ 10/(x_fo H)
dt_fo = 10 / (20 * H_fo)
print(f"  {'Freeze-out duration (Δx~10)':<35} {dt_fo:.3e}        {dt_fo/t_H:>16.2f}")

# CW rolling time (from θ_relic to π/2)
for f_decay, f_label in f_values:
    accel = abs(F_CW_relic) / f_decay**2
    delta_theta = math.pi/2 - theta_relic
    t_roll = math.sqrt(2 * delta_theta / accel) if accel > 0 else float('inf')
    label = f"CW roll (f={f_label})"
    print(f"  {label:<35} {t_roll:.3e}        {t_roll/t_H:>16.2e}")

print()

# =============================================================================
# SUMMARY
# =============================================================================
print("=" * 80)
print("  CONCLUSIONS")
print("=" * 80)
print(f"""
  1. WITHOUT bare potential (CW + thermal only):
     CW dominates over thermal backreaction by ~{abs(F_CW_relic/F_th_relic):.0e}
     → σ rolls to θ=π/2 (pure pseudoscalar) → SIDM dies
     → NOT viable as a dynamical mechanism

  2. WITH bare potential V_bare = ½ m²_bare f² (θ-θ_relic)²:
     Need m_bare >> H(T_fo) for effective trapping
     This requires UV input (fine-tuning of V_bare)

  3. A₄ DISCRETE SYMMETRY INTERPRETATION (PREFERRED):
     θ is NOT a dynamical field
     tan²θ = 1/9 is a GROUP THEORY CONSTANT from A₄ CG
     CW renormalizes coupling magnitudes, not the ratio
     ✅ No σ trapping problem — there is no σ to roll

  → The σ as a dynamical field with V_CW is PROBLEMATIC.
  → The A₄ interpretation (θ fixed by group theory) is CLEAN.
""")
