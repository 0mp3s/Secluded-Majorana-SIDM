"""
boltzmann_relic.py — Standalone Boltzmann solver for relic density (s-wave)
============================================================================
Extracted from core/v27_boltzmann_relic.py for Paper 2 (dark-energy-T-breaking).
No external dependencies (no core/, no global_config).

Functions exported:
    sigma_v_swave(alpha_d, m_chi, x_fo=20.0) -> float
    solve_boltzmann(m_chi, sv0, ...) -> (x_arr, Y_arr)
    Y_to_omega_h2(Y_inf, m_chi) -> float
"""
import math
import numpy as np

# ── Physical constants (hardcoded; sourced from global_config.json) ───────
M_PL         = 1.22089e19    # Planck mass [GeV]  (non-reduced, for freeze-out)
RHO_CRIT_H2  = 1.0539e-5    # ρ_crit h² [GeV cm⁻³]
S_0          = 2891.2        # entropy density today [cm⁻³]

# ── g_*(T) tabulation (Drees, Hajkarim & Schmitz 2015) ────────────────────
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
    logT = math.log(T) if T > 0 else -50
    return float(np.interp(logT, _LOG_T[::-1], _G_RHO[::-1]))


def g_star_S(T):
    """Effective relativistic d.o.f. g_*S(T) for entropy density."""
    logT = math.log(T) if T > 0 else -50
    return float(np.interp(logT, _LOG_T[::-1], _G_S[::-1]))


def Y_eq_full(x, m_chi, g_chi=2):
    """Equilibrium yield Y_eq (non-relativistic K₂ limit)."""
    T = m_chi / x
    g_s = g_star_S(T)
    if x > 300:
        return 0.0
    return 45.0 / (4 * math.pi**4) * g_chi / g_s * math.sqrt(math.pi / 2) * x**1.5 * math.exp(-x)


def _dYdx(x, Y, m_chi, sv0, g_chi=2):
    """dY/dx for s-wave: ⟨σv⟩ = sv0 = constant."""
    T = m_chi / x
    gs = g_star_S(T)
    gr = g_star_rho(T)
    g_eff = gs / math.sqrt(gr)
    lambda_x = math.sqrt(math.pi / 45.0) * g_eff * M_PL * m_chi
    Yeq = Y_eq_full(x, m_chi, g_chi)
    return -lambda_x * sv0 * (Y * Y - Yeq * Yeq) / (x * x)


# ── Public API ────────────────────────────────────────────────────────────

def sigma_v_swave(alpha_d, m_chi, x_fo=20.0):
    """s-wave ⟨σv⟩ for Majorana χχ → φφ (scalar mediator, t/u-channel),
    with Sommerfeld enhancement at freeze-out velocity.

        ⟨σv⟩ = ⟨σv⟩₀ × S(β_fo)
        S(β)  = (2πα/β) / (1 − e^{−2πα/β})   [Coulomb limit]
        ⟨σv⟩₀ = πα²/(4 m_χ²)
    """
    sv0 = math.pi * alpha_d**2 / (4.0 * m_chi**2)
    beta_fo = math.sqrt(3.0 / x_fo)
    xi = 2.0 * math.pi * alpha_d / beta_fo
    S_fo = xi / (1.0 - math.exp(-xi)) if xi > 1e-10 else 1.0
    return sv0 * S_fo


def solve_boltzmann(m_chi, sv0, x_start=1.0, x_end=1000.0, n_steps=10000, g_chi=2):
    """Solve dY/dx using RK4 for s-wave annihilation.

    Returns:
        x_arr, Y_arr  (numpy arrays)
    """
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
