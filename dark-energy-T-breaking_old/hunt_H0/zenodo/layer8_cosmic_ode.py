"""
Layer 8 — Cosmic ODE:  σ(t) ⊗ Friedmann  →  H₀
=================================================
SELF-CONTAINED VERSION for Zenodo upload.
No external dependencies beyond numpy/scipy/matplotlib.

The culmination of Layers 1-7: a coupled ODE that evolves the dark axion
field σ alongside the Friedmann equation from reheating to today.

    H₀ is OUTPUT, not input.

Physics
-------
σ is a pseudo-Nambu-Goldstone boson of dark SU(2) (dark pion), with
potential V(σ) = Λ_d⁴(1 − cos(σ/f)) from chiral symmetry breaking.

At reheating, σ is displaced by θ_i from the minimum (misalignment).
While H ≫ m_σ = Λ_d²/f, the field is frozen by Hubble friction.
As H drops to ∼m_σ ∼ H₀, σ begins its first oscillation.
Today: V(σ_today) ≈ ρ_Λ → H₀ emerges from H(t₀).

ODE system (in e-fold time N = ln a)
-------------------------------------
    dσ/dN = p
    dp/dN = −(3 − ε) p − V′(σ)/H²

    H² = (ρ_r + ρ_m + V(σ)) / (3 M_Pl² − ½p²)
    ε  = (⁴⁄₃ ρ_r + ρ_m + H²p²) / (2 M_Pl² H²)

where ρ_r ∝ a⁻⁴ (SM radiation+ν), ρ_m = ρ_χ + ρ_b ∝ a⁻³.

Inputs  (6 Lagrangian parameters)
---------------------------------
    m_χ, m_φ, α_D  — SIDM sector  → Ω_χ h² (Layer 7: Boltzmann)
    f, Λ_d, θ_i    — dark QCD     → V(σ)

Output
------
    H₀  [km/s/Mpc]  — derived from the Lagrangian alone.

Requirements
------------
    Python ≥ 3.9, numpy, scipy, matplotlib (for plotting)

License: CC-BY-4.0
Author:  Omer Pesach (ORCID: 0009-0004-8178-8131)
Date:    March 2026
DOI:     10.5281/zenodo.19225823
"""

import numpy as np
from scipy.integrate import solve_ivp
import math


# ═══════════════════════════════════════════════════════════════════════════
#  Physical constants  (natural units: ℏ = c = k_B = 1)
# ═══════════════════════════════════════════════════════════════════════════
M_PL        = 2.435e18            # reduced Planck mass  [GeV]
M_PL_FULL   = 1.220890e19        # full Planck mass     [GeV]
T_CMB_GEV   = 2.7255 * 8.617e-14 # CMB temperature      [GeV]
H_100_GEV   = 2.1332e-42         # 100 km/s/Mpc         [GeV]
OMEGA_B_H2  = 0.02237            # baryon density        (Planck 2018)
N_EFF       = 3.044              # effective neutrino species
S_0         = 2891.2             # entropy density today [cm⁻³]
RHO_CRIT_H2 = 1.0539e-5         # ρ_crit / h²          [GeV/cm³]

# Radiation density parameter  ω_r = Ω_r h²
_RHO_GAMMA  = np.pi**2 / 15.0 * T_CMB_GEV**4
_NU_FACTOR  = 1.0 + N_EFF * 7.0/8.0 * (4.0/11.0)**(4.0/3.0)
_RHO_UNIT   = 3.0 * M_PL**2 * H_100_GEV**2     # ρ_crit / h²  [GeV⁴]
OMEGA_R_H2  = _RHO_GAMMA * _NU_FACTOR / _RHO_UNIT

# Reference H₀ values
H0_PLANCK_KMS = 67.4
H0_SHOES_KMS  = 73.04
H0_PLANCK_GEV = H0_PLANCK_KMS * H_100_GEV / 100.0


# ═══════════════════════════════════════════════════════════════════════════
#  EMBEDDED BOLTZMANN SOLVER  (from core/v27_boltzmann_relic.py)
#  Self-contained: no imports from parent repository
# ═══════════════════════════════════════════════════════════════════════════

# g_*(T) tabulation (Drees, Hajkarim & Schmitz 2015)
_G_STAR_TABLE = np.array([
    [1e4,    106.75, 106.75],
    [200,    106.75, 106.75],
    [80,      86.25,  86.25],
    [10,      86.25,  86.25],
    [1,       75.75,  75.75],
    [0.3,     61.75,  61.75],
    [0.2,     17.25,  17.25],
    [0.15,    14.25,  14.25],
    [0.1,     10.75,  10.75],
    [0.01,    10.75,  10.75],
    [0.001,   10.75,  10.75],
    [0.0005,  10.75,  10.75],
    [0.0001,   3.36,   3.91],
    [1e-5,     3.36,   3.91],
    [1e-8,     3.36,   3.91],
])
_LOG_T = np.log(_G_STAR_TABLE[:, 0])
_G_RHO = _G_STAR_TABLE[:, 1]
_G_S   = _G_STAR_TABLE[:, 2]


def g_star_rho(T):
    """Effective relativistic d.o.f. g_*(T) for energy density."""
    logT = math.log(T) if T > 0 else -50.0
    return float(np.interp(logT, _LOG_T[::-1], _G_RHO[::-1]))


def g_star_S(T):
    """Effective relativistic d.o.f. g_*S(T) for entropy density."""
    logT = math.log(T) if T > 0 else -50.0
    return float(np.interp(logT, _LOG_T[::-1], _G_S[::-1]))


def Y_eq_full(x, m_chi, g_chi=2):
    """Equilibrium yield (non-relativistic K₂ limit)."""
    T = m_chi / x
    g_s = g_star_S(T)
    if x > 300:
        return 0.0
    return (45.0 / (4 * math.pi**4) * g_chi / g_s
            * math.sqrt(math.pi / 2) * x**1.5 * math.exp(-x))


def sigma_v_swave(alpha_d, m_chi, x_fo=20.0):
    """s-wave ⟨σv⟩ for Majorana χχ → φφ with Sommerfeld enhancement."""
    sv0 = math.pi * alpha_d**2 / (4.0 * m_chi**2)
    beta_fo = math.sqrt(3.0 / x_fo)
    xi = 2.0 * math.pi * alpha_d / beta_fo
    S_fo = xi / (1.0 - math.exp(-xi)) if xi > 1e-10 else 1.0
    return sv0 * S_fo


def _dYdx(x, Y, m_chi, sv0, g_chi=2):
    """Boltzmann equation RHS for s-wave annihilation."""
    T = m_chi / x
    gs = g_star_S(T)
    gr = g_star_rho(T)
    g_eff = gs / math.sqrt(gr)
    lambda_x = math.sqrt(math.pi / 45.0) * g_eff * M_PL_FULL * m_chi
    Yeq = Y_eq_full(x, m_chi, g_chi)
    return -lambda_x * sv0 * (Y * Y - Yeq * Yeq) / (x * x)


def solve_boltzmann(m_chi, sv0, x_start=1.0, x_end=1000.0,
                    n_steps=10000, g_chi=2):
    """Solve dY/dx using RK4 for s-wave annihilation."""
    dx = (x_end - x_start) / n_steps
    x_arr = np.zeros(n_steps + 1)
    Y_arr = np.zeros(n_steps + 1)
    x_arr[0] = x_start
    Y_arr[0] = Y_eq_full(x_start, m_chi, g_chi)

    for i in range(n_steps):
        x = x_arr[i]
        Y = Y_arr[i]
        k1 = dx * _dYdx(x, Y, m_chi, sv0, g_chi)
        k2 = dx * _dYdx(x + dx/2, Y + k1/2, m_chi, sv0, g_chi)
        k3 = dx * _dYdx(x + dx/2, Y + k2/2, m_chi, sv0, g_chi)
        k4 = dx * _dYdx(x + dx, Y + k3, m_chi, sv0, g_chi)
        Y_new = Y + (k1 + 2*k2 + 2*k3 + k4) / 6
        Y_new = max(Y_new, 1e-30)
        x_arr[i+1] = x + dx
        Y_arr[i+1] = Y_new

    return x_arr, Y_arr


def Y_to_omega_h2(Y_inf, m_chi):
    """Convert Y_∞ to Ω h²."""
    return m_chi * Y_inf * S_0 / RHO_CRIT_H2


# ═══════════════════════════════════════════════════════════════════════════
#  Dark axion potential:  V(σ) = Λ_d⁴ (1 − cos(σ/f))
# ═══════════════════════════════════════════════════════════════════════════
def V_sigma(sigma, f, Lambda_d):
    """Dark QCD chiral perturbation potential."""
    return Lambda_d**4 * (1.0 - np.cos(sigma / f))

def dV_sigma(sigma, f, Lambda_d):
    """dV/dσ = (Λ_d⁴/f) sin(σ/f)."""
    return Lambda_d**4 / f * np.sin(sigma / f)


# ═══════════════════════════════════════════════════════════════════════════
#  Layer 7: relic density from Boltzmann
# ═══════════════════════════════════════════════════════════════════════════
def compute_omega_chi_h2(m_chi, alpha_d):
    """Run Boltzmann solver for χχ → φφ (s-wave) → Ω_χ h²."""
    sv0 = sigma_v_swave(alpha_d, m_chi)
    x_arr, Y_arr = solve_boltzmann(m_chi, sv0, x_start=1, x_end=1000,
                                   n_steps=10000)
    return Y_to_omega_h2(Y_arr[-1], m_chi)


# ═══════════════════════════════════════════════════════════════════════════
#  ODE right-hand side
# ═══════════════════════════════════════════════════════════════════════════
def _ode_rhs(N, state, f, Lambda_d, rho_r0, rho_m0):
    """
    d/dN [σ, p]   where p ≡ dσ/dN.
    Friedmann is solved algebraically at each step.
    """
    sigma, p = state

    a = np.exp(N)

    rho_r = rho_r0 * a**(-4)
    rho_m = rho_m0 * a**(-3)
    V     = V_sigma(sigma, f, Lambda_d)

    denom = 3.0 * M_PL**2 - 0.5 * p**2
    if denom <= 0:
        return [0.0, 0.0]

    H2 = (rho_r + rho_m + V) / denom
    if H2 <= 0:
        return [0.0, 0.0]

    eps = (4.0/3.0 * rho_r + rho_m + H2 * p**2) / (2.0 * M_PL**2 * H2)

    dV = dV_sigma(sigma, f, Lambda_d)

    dsigma_dN = p
    dp_dN     = -(3.0 - eps) * p - dV / H2

    return [dsigma_dN, dp_dN]


# ═══════════════════════════════════════════════════════════════════════════
#  Result container
# ═══════════════════════════════════════════════════════════════════════════
class Layer8Result:
    """Container for Layer 8 output."""
    H0_kms      = None
    H0_GeV      = None
    h           = None
    sigma_today = None
    theta_today = None
    p_today     = None
    V_today     = None
    w_sigma     = None
    m_sigma     = None
    omega_chi_h2= None
    Omega_r     = None
    Omega_b     = None
    Omega_chi   = None
    Omega_DE    = None
    Omega_total = None
    sol         = None


# ═══════════════════════════════════════════════════════════════════════════
#  Main solver
# ═══════════════════════════════════════════════════════════════════════════
def solve_layer8(m_chi, m_phi_GeV, alpha_d, f, Lambda_d, theta_i,
                 T_RH=1e5, omega_chi_h2=None, verbose=True):
    """
    Solve σ(N) ⊗ Friedmann from reheating to today  →  H₀.

    Parameters
    ----------
    m_chi       : float – DM mass [GeV]
    m_phi_GeV   : float – mediator mass [GeV]  (documented; not used in v1)
    alpha_d     : float – dark coupling α_D
    f           : float – dark axion decay constant [GeV]
    Lambda_d    : float – dark QCD confinement scale [GeV]
    theta_i     : float – initial misalignment angle [rad]
    T_RH        : float – reheating temperature [GeV]  (default 10⁵)
    omega_chi_h2: float – if given, skip Boltzmann solver
    verbose     : bool

    Returns
    -------
    Layer8Result  with H₀, Ω breakdown, σ(today), full ODE solution
    """
    res = Layer8Result()
    m_sigma = Lambda_d**2 / f
    res.m_sigma = m_sigma

    if verbose:
        print("=" * 65)
        print("  LAYER 8:  σ(t) ⊗ Friedmann  →  H₀")
        print("=" * 65)
        print(f"\n  Inputs (6 Lagrangian parameters):")
        print(f"    m_χ   = {m_chi:.2f} GeV")
        print(f"    m_φ   = {m_phi_GeV*1e3:.2f} MeV")
        print(f"    α_D   = {alpha_d:.4e}")
        print(f"    f     = {f:.3e} GeV  ({f/M_PL:.4f} M_Pl)")
        print(f"    Λ_d   = {Lambda_d:.3e} GeV  ({Lambda_d*1e12:.4f} meV)")
        print(f"    θ_i   = {theta_i:.4f} rad  ({np.degrees(theta_i):.2f}°)")
        print(f"    m_σ   = Λ_d²/f = {m_sigma:.3e} GeV"
              f"  ({m_sigma/H0_PLANCK_GEV:.2f} × H₀_Planck)")

    # ── Layer 7: relic density ──
    if omega_chi_h2 is None:
        omega_chi_h2 = compute_omega_chi_h2(m_chi, alpha_d)
    res.omega_chi_h2 = omega_chi_h2

    if verbose:
        print(f"\n  Layer 7 (Boltzmann):  Ω_χ h² = {omega_chi_h2:.4f}")

    # ── Physical densities today (h²-independent) ──
    rho_r0 = OMEGA_R_H2  * _RHO_UNIT
    rho_m0 = (omega_chi_h2 + OMEGA_B_H2) * _RHO_UNIT

    # ── Initial conditions ──
    g_star_S_RH = 106.75
    g_star_S_0  = 3.91
    a_RH = (T_CMB_GEV / T_RH) * (g_star_S_0 / g_star_S_RH)**(1.0/3.0)
    N_RH = np.log(a_RH)

    sigma_init = f * theta_i
    p_init     = 0.0

    if verbose:
        print(f"  a_RH = {a_RH:.3e}    N_RH = {N_RH:.1f}")
        V_init = V_sigma(sigma_init, f, Lambda_d)
        print(f"  V(θ_i)  = {V_init:.4e} GeV⁴   (V^{{1/4}} = {V_init**0.25*1e12:.4f} meV)")
        H_RH = np.sqrt(rho_r0 * a_RH**(-4) / (3 * M_PL**2))
        print(f"  H(RH)   = {H_RH:.3e} GeV    m_σ/H(RH) = {m_sigma/H_RH:.2e}")
        print(f"\n  Integrating N = {N_RH:.0f}  →  0  ...")

    # ── Solve the coupled ODE ──
    sol = solve_ivp(
        _ode_rhs,
        [N_RH, 0.0],
        [sigma_init, p_init],
        args=(f, Lambda_d, rho_r0, rho_m0),
        method='RK45',
        rtol=1e-12,
        atol=1e-15,
        dense_output=True,
        max_step=1.0,
    )

    res.sol  = sol
    res.N_RH = N_RH

    if not sol.success:
        if verbose:
            print(f"  ⚠ ODE failed: {sol.message}")
        return res

    # ── Extract today's values at N = 0 ──
    sigma_0 = sol.y[0, -1]
    p_0     = sol.y[1, -1]
    V_0     = V_sigma(sigma_0, f, Lambda_d)

    denom   = 3.0 * M_PL**2 - 0.5 * p_0**2
    H0_sq   = (rho_r0 + rho_m0 + V_0) / denom
    H0_GeV  = np.sqrt(abs(H0_sq))
    H0_kms  = H0_GeV / H_100_GEV * 100.0
    h       = H0_kms / 100.0

    rho_total = 3.0 * M_PL**2 * H0_sq

    # σ equation of state
    rho_sig = 0.5 * H0_sq * p_0**2 + V_0
    P_sig   = 0.5 * H0_sq * p_0**2 - V_0
    w_sig   = P_sig / rho_sig if rho_sig > 0 else -1.0

    # Pack results
    res.H0_GeV      = H0_GeV
    res.H0_kms      = H0_kms
    res.h            = h
    res.sigma_today  = sigma_0
    res.theta_today  = sigma_0 / f
    res.p_today      = p_0
    res.V_today      = V_0
    res.w_sigma      = w_sig
    res.Omega_r      = OMEGA_R_H2  / h**2
    res.Omega_b      = OMEGA_B_H2  / h**2
    res.Omega_chi    = omega_chi_h2 / h**2
    res.Omega_DE     = rho_sig / rho_total
    res.Omega_total  = res.Omega_r + res.Omega_b + res.Omega_chi + res.Omega_DE

    if verbose:
        delta_theta = res.theta_today - theta_i
        frozen = "frozen" if abs(delta_theta / theta_i) < 1e-3 else "evolved!"
        print(f"  ODE: {sol.nfev} evals,  {len(sol.t)} output points")
        print(f"\n{'─'*65}")
        print(f"  RESULT")
        print(f"{'─'*65}")
        print(f"  sigma(today) = {sigma_0:.6e} GeV")
        print(f"  theta(today) = {res.theta_today:.6f} rad"
              f"   ({np.degrees(res.theta_today):.4f} deg)")
        print(f"  delta_theta  = {delta_theta:.3e} rad   [{frozen}]")
        print(f"  d_sigma/dN   = {p_0:.3e} GeV")
        print(f"  V(sigma_0)   = {V_0:.4e} GeV^4")
        print()
        print(f"  H0  =  {H0_kms:7.2f}  km/s/Mpc")
        print(f"  h   =  {h:.4f}")
        print(f"  w_sigma = {w_sig:.6f}   (-1 = CC)")
        print()
        print(f"  Omega breakdown:")
        print(f"    Omega_r   = {res.Omega_r:.6f}")
        print(f"    Omega_b   = {res.Omega_b:.4f}")
        print(f"    Omega_chi = {res.Omega_chi:.4f}")
        print(f"    Omega_DE  = {res.Omega_DE:.4f}   (sigma kinetic + potential)")
        print(f"    Sum Omega = {res.Omega_total:.6f}")
        print()
        dp = H0_kms - H0_PLANCK_KMS
        ds = H0_kms - H0_SHOES_KMS
        print(f"  vs Planck: Delta = {dp:+.2f} km/s/Mpc  ({dp/H0_PLANCK_KMS*100:+.2f}%)")
        print(f"  vs SH0ES:  Delta = {ds:+.2f} km/s/Mpc  ({ds/H0_SHOES_KMS*100:+.2f}%)")

    return res


# ═══════════════════════════════════════════════════════════════════════════
#  Utility: find Λ_d for a target H₀  (binary search)
# ═══════════════════════════════════════════════════════════════════════════
def find_Lambda_d_for_H0(target_H0_kms, m_chi, alpha_d, f, theta_i,
                          omega_chi_h2=None, tol=0.01):
    """Return (Lambda_d, Layer8Result) such that H0 ~ target."""
    if omega_chi_h2 is None:
        omega_chi_h2 = compute_omega_chi_h2(m_chi, alpha_d)

    h_t   = target_H0_kms / 100.0
    Om_m  = (omega_chi_h2 + OMEGA_B_H2 + OMEGA_R_H2) / h_t**2
    Om_L  = max(1.0 - Om_m, 0.01)
    rho_c = 3.0 * M_PL**2 * (target_H0_kms * H_100_GEV / 100.0)**2
    V_tgt = Om_L * rho_c
    Ld_static = (V_tgt / (1.0 - np.cos(theta_i)))**0.25

    Ld_vals = np.logspace(np.log10(Ld_static * 0.1),
                          np.log10(Ld_static * 10.0), 40)
    H0_vals = []
    for Ld in Ld_vals:
        r = solve_layer8(m_chi, 0, alpha_d, f, Ld, theta_i,
                         omega_chi_h2=omega_chi_h2, verbose=False)
        H0_vals.append(r.H0_kms if r.H0_kms is not None else 0.0)

    H0_vals = np.array(H0_vals)
    diffs = np.abs(H0_vals - target_H0_kms)
    best_idx = np.argmin(diffs)
    best_Ld  = Ld_vals[best_idx]

    if diffs[best_idx] < tol:
        r = solve_layer8(m_chi, 0, alpha_d, f, best_Ld, theta_i,
                         omega_chi_h2=omega_chi_h2, verbose=False)
        return best_Ld, r

    idx_lo = max(best_idx - 2, 0)
    idx_hi = min(best_idx + 2, len(Ld_vals) - 1)
    lo, hi = Ld_vals[idx_lo], Ld_vals[idx_hi]

    best_res = None
    best_Ld_final = best_Ld
    best_diff = diffs[best_idx]

    for _ in range(60):
        mid = np.sqrt(lo * hi)
        r = solve_layer8(m_chi, 0, alpha_d, f, mid, theta_i,
                         omega_chi_h2=omega_chi_h2, verbose=False)
        if r.H0_kms is None:
            hi = mid
            continue
        d = abs(r.H0_kms - target_H0_kms)
        if d < best_diff:
            best_diff = d
            best_Ld_final = mid
            best_res = r
        if d < tol:
            return mid, r
        if r.H0_kms > target_H0_kms:
            hi = mid
        else:
            lo = mid

    if best_res is None:
        best_res = solve_layer8(m_chi, 0, alpha_d, f, best_Ld_final, theta_i,
                                omega_chi_h2=omega_chi_h2, verbose=False)
    return best_Ld_final, best_res


# ═══════════════════════════════════════════════════════════════════════════
#  Plot: σ(N), H(N), ρ(N), w(N) evolution
# ═══════════════════════════════════════════════════════════════════════════
def plot_evolution(res, f, Lambda_d, save_path=None):
    """Four-panel diagnostic of the cosmic evolution."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available -- skipping plot")
        return

    sol   = res.sol
    N_arr = np.linspace(sol.t[0], sol.t[-1], 3000)
    Y_arr = sol.sol(N_arr)
    sigma_arr = Y_arr[0]
    p_arr     = Y_arr[1]

    rho_r0 = OMEGA_R_H2  * _RHO_UNIT
    rho_m0 = (res.omega_chi_h2 + OMEGA_B_H2) * _RHO_UNIT

    a_arr = np.exp(N_arr)
    rho_r = rho_r0 * a_arr**(-4)
    rho_m = rho_m0 * a_arr**(-3)
    V_arr = V_sigma(sigma_arr, f, Lambda_d)

    denom  = 3.0 * M_PL**2 - 0.5 * p_arr**2
    H2_arr = np.where(denom > 0, (rho_r + rho_m + V_arr) / denom, 1e-100)
    H_arr  = np.sqrt(np.maximum(H2_arr, 0))
    H_kms  = H_arr / H_100_GEV * 100.0

    rho_sig = 0.5 * H2_arr * p_arr**2 + V_arr
    P_sig   = 0.5 * H2_arr * p_arr**2 - V_arr
    w_arr   = np.where(rho_sig > 0, P_sig / rho_sig, -1.0)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(r"Layer 8:  $\sigma(N) \otimes$ Friedmann  $\to$  $H_0$"
                 f"  =  {res.H0_kms:.2f} km/s/Mpc",
                 fontsize=14, fontweight='bold')

    ax = axes[0, 0]
    ax.plot(N_arr, sigma_arr / f, 'b-', lw=1.5)
    ax.set_xlabel(r'$N = \ln a$')
    ax.set_ylabel(r'$\theta = \sigma/f$  [rad]')
    ax.set_title('Dark axion angle')

    ax = axes[0, 1]
    mask = N_arr > -20
    ax.semilogy(N_arr[mask], H_kms[mask], 'k-', lw=1.5)
    ax.axhline(H0_PLANCK_KMS, color='dodgerblue', ls='--', alpha=0.6,
               label=f'Planck {H0_PLANCK_KMS}')
    ax.axhline(H0_SHOES_KMS, color='crimson', ls='--', alpha=0.6,
               label=f'SH0ES {H0_SHOES_KMS}')
    ax.set_xlabel(r'$N = \ln a$')
    ax.set_ylabel(r'$H$  [km/s/Mpc]')
    ax.set_title('Hubble parameter')
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    mask = N_arr > -25
    ax.semilogy(N_arr[mask], rho_r[mask], color='gold', lw=1.2, label=r'$\rho_r$')
    ax.semilogy(N_arr[mask], rho_m[mask], color='royalblue', lw=1.2, label=r'$\rho_m$')
    ax.semilogy(N_arr[mask], V_arr[mask], color='crimson', lw=1.5, label=r'$V(\sigma)$')
    ax.set_xlabel(r'$N = \ln a$')
    ax.set_ylabel(r'$\rho$  [GeV$^4$]')
    ax.set_title('Energy densities')
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    mask2 = N_arr > -5
    ax.plot(N_arr[mask2], w_arr[mask2], color='purple', lw=1.5)
    ax.axhline(-1, color='gray', ls='--', alpha=0.5, label=r'$w=-1$ (CC)')
    ax.axhline(-1.0/3.0, color='gray', ls=':', alpha=0.4, label=r'$w=-1/3$')
    ax.set_xlabel(r'$N = \ln a$')
    ax.set_ylabel(r'$w_\sigma$')
    ax.set_title(r'Dark energy equation of state')
    ax.set_ylim(-1.1, 0.5)
    ax.legend(fontsize=9)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Plot saved: {save_path}")
    else:
        plt.show()


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN — Quick demonstration
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":

    # MAP benchmark (from SIDM global_config.json)
    m_chi   = 98.19             # GeV
    m_phi   = 9.66e-3           # GeV  (9.66 MeV)
    alpha_d = 3.274e-3

    # Dark QCD sector
    theta_i = 3.0               # rad   (near pi — hilltop quintessence)
    f_ref   = 0.27 * M_PL      # GeV   (from PI-7 analysis)

    # Use observed relic density
    omega_chi = 0.120

    print("=" * 65)
    print("  Layer 8 — Self-contained Zenodo version")
    print("  Deriving H0 from the Lagrangian")
    print("=" * 65)

    # Run solver
    res = solve_layer8(m_chi, m_phi, alpha_d, f_ref, 2.0e-12, theta_i,
                       omega_chi_h2=omega_chi)

    # Plot
    import os
    plot_path = os.path.join(os.path.dirname(__file__), "layer8_evolution.png")
    if res.sol is not None and res.sol.success:
        plot_evolution(res, f_ref, 2.0e-12, save_path=plot_path)
