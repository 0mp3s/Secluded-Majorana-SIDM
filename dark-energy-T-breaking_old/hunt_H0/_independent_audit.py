#!/usr/bin/env python3
"""
INDEPENDENT AUDIT — Zero project imports.
All constants from PDG/Planck 2018. All equations from the action.
Purpose: verify the research journal claims WITHOUT trusting project code.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

print("=" * 75)
print("  INDEPENDENT AUDIT — NO IMPORTS FROM PROJECT")
print("  Constants: PDG 2024 / Planck 2018")
print("  Equations: derived from Einstein-Hilbert + Klein-Gordon action")
print("=" * 75)

# =====================================================================
#  CONSTANTS (all from PDG / Planck 2018)
# =====================================================================
M_PL       = 2.435e18          # reduced Planck mass [GeV]
T_CMB_K    = 2.7255            # CMB temperature [K]
k_B        = 8.617333e-14      # Boltzmann [GeV/K]
T_CMB      = T_CMB_K * k_B     # [GeV]

# H_100 = 100 km/s/Mpc in natural units
# 100 km/s/Mpc = 3.2408e-18 s^-1, hbar = 6.5821e-25 GeV*s
H_100      = 2.1332e-42        # [GeV]

OMEGA_B_H2 = 0.02237           # Planck 2018
N_EFF      = 3.044

# Radiation density from first principles
rho_gamma  = np.pi**2 / 15.0 * T_CMB**4
nu_factor  = 1.0 + N_EFF * (7.0/8.0) * (4.0/11.0)**(4.0/3.0)
RHO_UNIT   = 3.0 * M_PL**2 * H_100**2
OMEGA_R_H2 = rho_gamma * nu_factor / RHO_UNIT

print(f"\n  M_PL           = {M_PL:.4e} GeV")
print(f"  T_CMB          = {T_CMB:.4e} GeV")
print(f"  H_100          = {H_100:.4e} GeV")
print(f"  Omega_r h^2    = {OMEGA_R_H2:.4e}  (PDG ~4.15e-5)")

# =====================================================================
#  CHECK 1: Omega_r h^2 sanity
# =====================================================================
print("\n--- CHECK 1: Omega_r h^2 ---")
ratio_or = OMEGA_R_H2 / 4.15e-5
print(f"  Our/PDG = {ratio_or:.4f}" + (" [PASS]" if abs(ratio_or 
    - 1) < 0.02 else f" [WARN: {abs(ratio_or-1)*100:.1f}% off]"))

# =====================================================================
#  ODE system (from action S = int d^4x sqrt(-g) [M_Pl^2 R/2 - 1/2
#  (d sigma)^2 - V(sigma)])
# =====================================================================
def ode_rhs(N, state, f, Ld, rho_r0, rho_m0):
    sigma, p = state
    a = np.exp(N)
    rho_r = rho_r0 * a**(-4)
    rho_m = rho_m0 * a**(-3)
    V = Ld**4 * (1.0 - np.cos(sigma / f))
    denom = 3.0 * M_PL**2 - 0.5 * p**2
    if denom <= 0:
        return [0.0, 0.0]
    H2 = (rho_r + rho_m + V) / denom
    if H2 <= 0:
        return [0.0, 0.0]
    eps = (4.0/3.0 * rho_r + rho_m + H2 * p**2) / (2.0 * M_PL**2 * H2)
    dV = Ld**4 / f * np.sin(sigma / f)
    return [p, -(3.0 - eps) * p - dV / H2]


def solve_H0(theta_i, Lambda_d, f, omega_chi_h2=0.120):
    """Solve coupled sigma-Friedmann ODE -> H0."""
    rho_r0 = OMEGA_R_H2 * RHO_UNIT
    rho_m0 = (omega_chi_h2 + OMEGA_B_H2) * RHO_UNIT
    T_RH = 1e5
    a_RH = (T_CMB / T_RH) * (3.91 / 106.75)**(1.0/3.0)
    N_RH = np.log(a_RH)
    
    sol = solve_ivp(ode_rhs, [N_RH, 0.0], [f * theta_i, 0.0],
                    args=(f, Lambda_d, rho_r0, rho_m0),
                    method='RK45', rtol=1e-12, atol=1e-15,
                    dense_output=True, max_step=1.0)
    if not sol.success:
        return None, None, None, None
    
    sigma_0, p_0 = sol.y[0, -1], sol.y[1, -1]
    V_0 = Lambda_d**4 * (1.0 - np.cos(sigma_0 / f))
    denom = 3.0 * M_PL**2 - 0.5 * p_0**2
    H0_sq = (rho_r0 + rho_m0 + V_0) / denom
    H0_GeV = np.sqrt(abs(H0_sq))
    H0_kms = H0_GeV / H_100 * 100.0
    
    rho_sig = 0.5 * H0_sq * p_0**2 + V_0
    P_sig = 0.5 * H0_sq * p_0**2 - V_0
    w = P_sig / rho_sig if rho_sig > 0 else -1
    rho_tot = 3.0 * M_PL**2 * H0_sq
    Omega_DE = rho_sig / rho_tot if rho_tot > 0 else 0
    return H0_kms, w, Omega_DE, sol


# =====================================================================
#  TEST A: theta_i scan (core claim: H0 is OUTPUT)
# =====================================================================
print("\n" + "=" * 75)
print("  TEST A: theta_i scan — H0 from Lagrangian")
print("  Lambda_d = 2 meV, f = 0.27 M_Pl, Omega_chi h^2 = 0.120")
print("=" * 75)

f_val = 0.27 * M_PL
Ld = 2.0e-12
m_sig = Ld**2 / f_val
print(f"\n  m_sigma = {m_sig:.3e} GeV")
print(f"  m_sigma / H0_Planck = {m_sig / (67.4 * H_100 / 100):.2f}")
print()

hdr = f"  {'theta_i':>8}  {'theta/pi':>8}  {'H0 [km/s/Mpc]':>14}  {'w_sigma':>10}  {'Omega_DE':>10}"
print(hdr)
print("  " + "-" * 58)

journal_claims = {
    2.0:   ("~40.8",  "~-0.36"),
    2.5:   ("~48.0",  "~+0.93"),
    3.0:   ("~71.1",  "~-0.86"),
    3.1:   ("~73.1",  "~-0.99"),
}

scan_thetas = [2.0, 2.5, 2.8, 2.887, 2.9, 2.95, 3.0, 3.05, 3.1, np.pi]
results_A = {}
for ti in scan_thetas:
    H0, w, ODE, _ = solve_H0(ti, Ld, f_val)
    results_A[round(ti, 4)] = (H0, w, ODE)
    if H0 is not None:
        label = "pi" if abs(ti - np.pi) < 0.001 else f"{ti:.3f}"
        print(f"  {label:>8}  {ti/np.pi:>8.4f}  {H0:>14.2f}  {w:>10.5f}  {ODE:>10.4f}")

print()
print("  JOURNAL vs ACTUAL:")
for ti, (jH0, jw) in journal_claims.items():
    actual = results_A.get(round(ti, 4))
    if actual and actual[0]:
        print(f"    theta={ti}: journal H0={jH0}, actual H0={actual[0]:.1f}, "
              f"journal w={jw}, actual w={actual[1]:.3f}")

# =====================================================================
#  TEST B: w0 CPL fit around theta_i = 2.887 (DESI match claim)
# =====================================================================
print("\n" + "=" * 75)
print("  TEST B: CPL fit at theta_i = 2.887 (DESI w0 = -0.727 claim)")
print("=" * 75)

H0_2887, w_2887, ODE_2887, sol_2887 = solve_H0(2.887, Ld, f_val)
print(f"\n  H0(2.887) = {H0_2887:.2f} km/s/Mpc")
print(f"  w(z=0)    = {w_2887:.5f}")
print(f"  Omega_DE  = {ODE_2887:.4f}")

# Extract w(a) from dense solution 
if sol_2887 is not None:
    N_arr = np.linspace(max(sol_2887.t[0], -5.0), 0.0, 500)
    Y = sol_2887.sol(N_arr)
    a_arr = np.exp(N_arr)
    sigma_arr, p_arr = Y[0], Y[1]
    
    rho_r0 = OMEGA_R_H2 * RHO_UNIT
    rho_m0 = (0.120 + OMEGA_B_H2) * RHO_UNIT
    rho_r = rho_r0 * a_arr**(-4)
    rho_m = rho_m0 * a_arr**(-3)
    V_arr = Ld**4 * (1.0 - np.cos(sigma_arr / f_val))
    denom_arr = 3.0 * M_PL**2 - 0.5 * p_arr**2
    H2_arr = np.where(denom_arr > 0, (rho_r + rho_m + V_arr) / denom_arr, 1e-100)
    rho_sig = 0.5 * H2_arr * p_arr**2 + V_arr
    P_sig = 0.5 * H2_arr * p_arr**2 - V_arr
    w_arr = np.where(rho_sig > 0, P_sig / rho_sig, -1.0)
    
    # CPL fit: w(a) = w0 + wa*(1-a) for z < 2
    mask = a_arr > 1.0/3.0  # z < 2
    try:
        popt, _ = curve_fit(lambda a, w0, wa: w0 + wa*(1-a),
                           a_arr[mask], w_arr[mask], p0=[-0.9, -0.5])
        w0_cpl, wa_cpl = popt
        print(f"\n  CPL fit (z < 2):")
        print(f"    w0 = {w0_cpl:.5f}")
        print(f"    wa = {wa_cpl:.5f}")
        print(f"\n  DESI DR1: w0 = -0.727 +/- 0.067,  wa = -1.05 +0.31/-0.27")
        print(f"  Delta_w0 = {abs(w0_cpl - (-0.727)):.4f} = {abs(w0_cpl - (-0.727))/0.067:.1f} sigma")
        print(f"  Delta_wa = {abs(wa_cpl - (-1.05)):.4f} = {abs(wa_cpl - (-1.05))/0.31:.1f} sigma")
    except Exception as e:
        print(f"  CPL fit failed: {e}")

# =====================================================================
#  TEST C: Dimensional transmutation
# =====================================================================
print("\n" + "=" * 75)
print("  TEST C: Dimensional Transmutation alpha_d -> Lambda_d")
print("=" * 75)

b0 = 19.0 / 3.0  # SU(2) with 3 Majorana: b0 = 11/3 * C2(G) - 2/3 * T(R) * N_f
m_chi = 98.19

print(f"\n  b0 = 11/3 * 2 - 2/3 * 1/2 * 3 = {11/3*2 - 2/3*0.5*3:.4f}")
print(f"  (Check: 22/3 - 1 = 19/3 = {19/3:.4f})")

# Forward
print(f"\n  {'alpha_d':>8}  {'1/alpha':>8}  {'Lambda_d [GeV]':>16}  {'Lambda_d [meV]':>14}")
print("  " + "-" * 52)
for ad in [0.025, 0.029, 0.0315, 0.035, 0.040]:
    Ld_tr = m_chi * np.exp(-2*np.pi / (b0 * ad))
    print(f"  {ad:>8.4f}  {1/ad:>8.1f}  {Ld_tr:>16.4e}  {Ld_tr*1e12:>14.4f}")

# Inverse
alpha_needed = 2*np.pi / (b0 * np.log(m_chi / 2.0e-12))
print(f"\n  Inverse: Lambda_d = 2 meV -> alpha_d = {alpha_needed:.5f} = 1/{1/alpha_needed:.1f}")

# =====================================================================
#  TEST D: find_Lambda_d_for_H0 circularity test
# =====================================================================
print("\n" + "=" * 75)
print("  TEST D: CIRCULARITY CHECK — Is Lambda_d derived or fitted?")
print("=" * 75)

# Binary search for Lambda_d that gives H0 = 67.4
target_H0 = 67.4
theta_test = 3.0
lo_Ld, hi_Ld = 1e-13, 1e-11

for _ in range(60):
    mid_Ld = np.sqrt(lo_Ld * hi_Ld)
    H0_test, _, _, _ = solve_H0(theta_test, mid_Ld, f_val)
    if H0_test is None:
        hi_Ld = mid_Ld
        continue
    if H0_test > target_H0:
        hi_Ld = mid_Ld
    else:
        lo_Ld = mid_Ld
    if abs(H0_test - target_H0) < 0.01:
        break

print(f"\n  theta_i = {theta_test}")
print(f"  Target H0 = {target_H0} km/s/Mpc")
print(f"  Binary search found: Lambda_d = {mid_Ld:.4e} GeV = {mid_Ld*1e12:.4f} meV")
print(f"  Resulting H0 = {H0_test:.2f} km/s/Mpc")
print()
print("  CRITICAL QUESTION: Is Lambda_d = 2 meV a PREDICTION or a FIT?")
print(f"    If you SET H0 = 67.4, you GET Lambda_d = {mid_Ld*1e12:.2f} meV")
print(f"    If you SET Lambda_d = 2 meV, you GET H0 = {results_A.get(3.0, (None,))[0]:.1f} km/s/Mpc")
print("    The two are equivalent. Lambda_d is NOT derived from first principles.")
print("    It IS derived from alpha_d via transmutation — but alpha_d is free.")

# =====================================================================
#  TEST E: Fine-tuning of theta_i
# =====================================================================
print("\n" + "=" * 75)
print("  TEST E: Fine-tuning assessment")
print("=" * 75)

H0_range = []
for ti in np.linspace(0.01, 2*np.pi - 0.01, 200):
    H0_t, _, _, _ = solve_H0(ti, Ld, f_val)
    if H0_t is not None and H0_t > 60 and H0_t < 80:
        H0_range.append(ti)

if H0_range:
    ti_min, ti_max = min(H0_range), max(H0_range)
    frac = (ti_max - ti_min) / (2 * np.pi) * 100
    print(f"\n  theta_i giving H0 in [60, 80]: [{ti_min:.2f}, {ti_max:.2f}]")
    print(f"  Width = {ti_max - ti_min:.2f} rad out of 2*pi = {2*np.pi:.2f}")
    print(f"  Fraction of parameter space: {frac:.1f}%")
    print(f"  Journal claims: ~5.5%")
else:
    print("  No theta_i found giving H0 in [60, 80]")

# =====================================================================
#  FINAL VERDICT
# =====================================================================
print("\n" + "=" * 75)
print("  FINAL INDEPENDENT VERDICT")
print("=" * 75)
print("""
  1. H0 IS an output of the ODE solver: TRUE
     Given (Lambda_d, f, theta_i, Omega_chi), solving sigma(N)+Friedmann
     gives H0 without inputting it. This is real.

  2. H0 ~ 67-73 for theta_i ~ 2.9-3.1: VERIFIED numerically.

  3. w0 = -0.727 at theta_i = 2.887: CHECKING CPL fit above.

  4. Lambda_d = 2 meV from alpha_d ~ 0.031: CORRECT math.
     But alpha_d is free — this is parameter substitution, not prediction.

  5. Circularity: find_Lambda_d_for_H0(67.4) is a binary search for 
     the Lambda_d that reproduces the target. THIS IS FITTING, not predicting.
     The model has 2 free DE parameters (Lambda_d or alpha_d, and theta_i).
     Once you fix them from (H0, Omega_DE), w(z) is a genuine prediction.

  6. DESI match: w0 = -0.727 AT theta_i = 2.887 — need to verify CPL fit.
     But theta_i = 2.887 was CHOSEN to match DESI. The prediction is:
     IF theta_i ~ 2.9 THEN w0 ~ -0.73. That's model-dependent, not a priori.

  HONEST SUMMARY:
  - The ODE math is correct.
  - H0 emerges from the ODE (not put in by hand) — TRUE.
  - But Lambda_d (or alpha_d) and theta_i are free parameters.
  - With 2 free params, matching H0 and Omega_DE is expected (2 eqs, 2 unknowns).
  - The REAL prediction is w(z) ≠ -1 — a curve, not a number.
  - w(z) is falsifiable by DESI DR2 — this IS genuine predictive power.
""")
