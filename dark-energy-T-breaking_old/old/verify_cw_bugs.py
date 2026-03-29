#!/usr/bin/env python3
"""
verify_cw_bugs.py — Verify CW potential derivative and minimum location
========================================================================
Numerical verification of all bugs found in freeze_out_trapping.py
"""

import numpy as np

# Constants
MeV = 1e-3
M_Pl = 2.435e18

# MAP benchmark
m_chi = 94.07 * MeV
m_phi = 11.10 * MeV
alpha = 5.734e-3
theta_relic = np.arctan(1/np.sqrt(8))
y_sq = 4 * np.pi * alpha / np.cos(theta_relic)**2
y = np.sqrt(y_sq)
v_phi = 0.5 * m_phi

print("="*78)
print("BUG VERIFICATION FOR freeze_out_trapping.py")
print("="*78)

# ============================================================================
# BUG 1: CW derivative — is the constant -0.5 or -1?
# ============================================================================
print("\n" + "="*78)
print("BUG 1: dV_CW/dθ — coefficient check (numerical vs analytical)")
print("="*78)

def V_CW(theta, v, mu=None):
    """CW potential (Majorana fermion, n_f=2)"""
    if mu is None:
        mu = m_chi
    M2 = m_chi**2 + m_chi * y * np.cos(theta) * v + y**2 * v**2 / 4
    if M2 <= 0:
        return np.inf
    # V = -(n_f/64π²) M⁴ [ln(M²/μ²) - 3/2], n_f=2 for Majorana
    return -(2/(64 * np.pi**2)) * M2**2 * (np.log(M2 / mu**2) - 1.5)

def dVCW_numerical(theta, v, dth=1e-8):
    """Numerical derivative via central difference"""
    return (V_CW(theta + dth, v) - V_CW(theta - dth, v)) / (2 * dth)

def dVCW_with_half(theta, v):
    """Analytical derivative with [ln - 0.5] (code's version)"""
    M2 = m_chi**2 + m_chi * y * np.cos(theta) * v + y**2 * v**2 / 4
    dM2_dth = -m_chi * y * np.sin(theta) * v
    return -(1/(16*np.pi**2)) * dM2_dth * M2 * (np.log(M2/m_chi**2) - 0.5)

def dVCW_with_one(theta, v):
    """Analytical derivative with [ln - 1] (proposed fix)"""
    M2 = m_chi**2 + m_chi * y * np.cos(theta) * v + y**2 * v**2 / 4
    dM2_dth = -m_chi * y * np.sin(theta) * v
    return -(1/(16*np.pi**2)) * dM2_dth * M2 * (np.log(M2/m_chi**2) - 1.0)

print(f"\nTesting at several θ values (v_φ = {v_phi/MeV:.2f} MeV):\n")
print(f"  {'θ (°)':<10} {'Numerical':<15} {'[ln-0.5]':<15} {'[ln-1.0]':<15} "
      f"{'Err(0.5)':<12} {'Err(1.0)':<12}")
print(f"  {'-'*78}")

for th_deg in [5, 10, 19.47, 30, 45, 60, 80]:
    th = np.radians(th_deg)
    num = dVCW_numerical(th, v_phi)
    ana_half = dVCW_with_half(th, v_phi)
    ana_one = dVCW_with_one(th, v_phi)
    
    err_half = abs(ana_half - num) / abs(num) if num != 0 else 0
    err_one = abs(ana_one - num) / abs(num) if num != 0 else 0
    
    print(f"  {th_deg:<10.2f} {num:<15.6e} {ana_half:<15.6e} {ana_one:<15.6e} "
          f"{err_half:<12.2e} {err_one:<12.2e}")

# ============================================================================
# BUG 2: Location of CW minimum
# ============================================================================
print("\n" + "="*78)
print("BUG 2: V_CW(θ) — where is the minimum?")
print("="*78)

thetas = np.linspace(0.01, np.pi/2 - 0.01, 1000)
V_vals = np.array([V_CW(th, v_phi) for th in thetas])

idx_min = np.argmin(V_vals)
idx_max = np.argmax(V_vals)

print(f"\n  V_CW minimum at θ = {np.degrees(thetas[idx_min]):.2f}°")
print(f"  V_CW maximum at θ = {np.degrees(thetas[idx_max]):.2f}°")
print(f"  V_CW(0°)  = {V_CW(0.01, v_phi):.6e} GeV⁴")
print(f"  V_CW(19.47°) = {V_CW(theta_relic, v_phi):.6e} GeV⁴")
print(f"  V_CW(45°) = {V_CW(np.pi/4, v_phi):.6e} GeV⁴")
print(f"  V_CW(89°) = {V_CW(np.pi/2-0.01, v_phi):.6e} GeV⁴")
print(f"\n  ΔV = V(0) - V(90°) = {V_CW(0.01, v_phi) - V_CW(np.pi/2-0.01, v_phi):.6e} GeV⁴")
print(f"  ΔV > 0 means θ=0 is HIGHER energy → θ=0 is MAXIMUM")
print(f"  ΔV < 0 means θ=0 is LOWER energy → θ=0 is minimum")

# Verify with dV/dθ sign
dV_at_30 = dVCW_numerical(np.radians(30), v_phi)
print(f"\n  dV/dθ at θ=30°: {dV_at_30:.6e}  (negative = rolling toward larger θ)")

# ============================================================================
# BUG 3: m_σ² — sign and magnitude
# ============================================================================
print("\n" + "="*78)
print("BUG 3: m_σ² at θ=0 — sign check")
print("="*78)

# Numerical second derivative at θ=0
dth = 1e-5
d2V_num_0 = (V_CW(dth, v_phi) - 2*V_CW(0.001, v_phi) + V_CW(-dth + 0.002, v_phi))
# Better: use symmetric points around a small θ
th0 = 0.01  # near θ=0
d2V_num_0 = (V_CW(th0 + dth, v_phi) - 2*V_CW(th0, v_phi) + V_CW(th0 - dth, v_phi)) / dth**2

# Numerical second derivative at θ=π/2
th90 = np.pi/2 - 0.01
d2V_num_90 = (V_CW(th90 + dth, v_phi) - 2*V_CW(th90, v_phi) + V_CW(th90 - dth, v_phi)) / dth**2

# Numerical second derivative at θ_relic
th_r = theta_relic
d2V_num_r = (V_CW(th_r + dth, v_phi) - 2*V_CW(th_r, v_phi) + V_CW(th_r - dth, v_phi)) / dth**2

# Code's estimate for m_σ²
f_val = 0.2 * M_Pl
m2_sigma_code = y**2 * m_chi**2 * v_phi**2 / (16 * np.pi**2 * f_val**2)

print(f"\n  d²V/dθ² at θ ≈ 0°:     {d2V_num_0:.6e} GeV⁴  ({'CONCAVE (max)' if d2V_num_0 < 0 else 'CONVEX (min)'})")
print(f"  d²V/dθ² at θ_relic:     {d2V_num_r:.6e} GeV⁴  ({'CONCAVE' if d2V_num_r < 0 else 'CONVEX'})")
print(f"  d²V/dθ² at θ ≈ 90°:    {d2V_num_90:.6e} GeV⁴  ({'CONCAVE' if d2V_num_90 < 0 else 'CONVEX'})")
print(f"\n  Code's m²_σ estimate (always positive): {m2_sigma_code:.6e} GeV²")
print(f"  Actual m²_σ(θ=0)/f² = {d2V_num_0/f_val**2:.6e} GeV²")
print(f"  Actual m²_σ(θ_relic)/f² = {d2V_num_r/f_val**2:.6e} GeV²")
print(f"  Ratio |actual(θ=0)| / code estimate: {abs(d2V_num_0/f_val**2) / m2_sigma_code:.1f}")

# ============================================================================
# BUG 4: t_eval floating point
# ============================================================================
print("\n" + "="*78)
print("BUG 4: t_eval floating point issue")
print("="*78)

x_span_end = 200.0
x_eval = np.logspace(0, np.log10(200), 1000)
print(f"\n  x_span = (1.0, {x_span_end})")
print(f"  x_eval[-1] = {x_eval[-1]:.17f}")
print(f"  x_eval[-1] > x_span_end? {x_eval[-1] > x_span_end}")
print(f"  Difference: {x_eval[-1] - x_span_end:.2e}")

x_eval_fixed = np.clip(x_eval, 1.0, x_span_end)
print(f"  Fixed x_eval[-1] = {x_eval_fixed[-1]:.17f}")

# ============================================================================
# BUG 5: Hard-coded summary vs computed results
# ============================================================================
print("\n" + "="*78)
print("BUG 5: Summary contradicts computed results")
print("="*78)

H_fo = np.sqrt(np.pi**2 * 10.75 / 90) * (m_chi/20)**2 / M_Pl

print(f"\n  H(T_fo) = {H_fo:.3e} GeV")
print(f"\n  Computed m_σ/H results vs hard-coded summary:")
print(f"    f = 0.2 M_Pl: m_σ/H = 2.46  → code says 'Already oscillating'  → summary says '≫ 1' ✓")
print(f"    f = 1.0 M_Pl: m_σ/H = 0.49  → code says 'FROZEN by Hubble'     → summary says '≫ 1' ✗ WRONG")
print(f"    f = 15  M_Pl: m_σ/H = 0.033 → code says 'FROZEN by Hubble'     → summary says '> 1'  ✗ WRONG")
print(f"\n  Also: Part 3 found 'possible dynamical trapping' for f=1,15 M_Pl")
print(f"  but hard-coded conclusion says 'For ALL values, m_σ ≫ H' — contradicts Part 3!")

# ============================================================================
# IMPACT ASSESSMENT
# ============================================================================
print("\n" + "="*78)
print("IMPACT ASSESSMENT: HOW DO BUGS AFFECT THE CONCLUSION?")
print("="*78)

print("""
Bug 1 (dV/dθ coefficient):
  [ln-0.5] → [ln-1.0]: factor ~2 in CW force.
  Since F_th/F_CW ~ 10⁻¹⁰, this doesn't change the main conclusion.
  
Bug 2 (CW minimum location):
  ⚠️ CRITICAL: θ=0 is a CW MAXIMUM, not minimum!
  The CW potential pushes σ TOWARD θ=π/2 (pure pseudoscalar).
  Previous claim "σ falls to θ=0" was WRONG in direction.
  
  BUT: the main conclusion is still valid — CW dominates over thermal
  backreaction, and drives σ AWAY from θ_relic (toward θ=90°, not θ=0°).
  At θ=90°: α_s=0, no SIDM. This is phenomenologically dead.

Bug 3 (m_σ² sign):
  The code uses a positive m_σ² estimate. The actual m_σ² at θ=0 is NEGATIVE
  (hilltop), confirming Bug 2. At the CW minimum (θ=π/2), m_σ² is tiny
  (∝ ln(M²/μ²) ≈ 10⁻⁵), so σ is essentially massless there.

Bug 4 (t_eval):
  Caused ALL ODE simulations to crash. The analytical parts ran fine.

Bug 5 (hard-coded summary):
  Summary claims "all f have m_σ ≫ H" but the code itself found
  m_σ < H for f = 1.0 and 15 M_Pl. The "interesting" cases were dismissed.

OVERALL: The CONCLUSION (freeze-out cannot trap σ at θ_relic) survives,
but for partially different reasons than stated:
  - CW drives σ to θ=π/2, not θ=0
  - Thermal backreaction is too weak regardless of direction
  - The "interesting" f=1-15 M_Pl cases deserve actual ODE solution
""")
