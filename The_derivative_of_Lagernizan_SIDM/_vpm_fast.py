"""
Numba-accelerated VPM cross section with configurable partial-wave weights.

Reuses the Numba-JIT phase shift solver from core/v22_raw_scan.py,
giving ~100-1000x speedup over the pure-Python reimplementation in T9-T12.

Usage:
    from _vpm_fast import sigma_T_weighted, get_phase_shifts
    sigma = sigma_T_weighted(m_chi, m_phi, alpha, v_km_s, 1.0, 3.0)  # Majorana
    sigma = sigma_T_weighted(m_chi, m_phi, alpha, v_km_s, 1.0, 1.0)  # Dirac
"""
import math
import sys
from pathlib import Path

from numba import jit

# ── path setup (ensure core/ is importable) ──────────────────────────
_HERE = Path(__file__).resolve().parent
_SIDM_ROOT = _HERE.parent
_CORE = _SIDM_ROOT / "core"
if str(_SIDM_ROOT) not in sys.path:
    sys.path.insert(0, str(_SIDM_ROOT))
if str(_CORE) not in sys.path:
    sys.path.insert(0, str(_CORE))

from core.v22_raw_scan import (
    vpm_phase_shift,
    C_KM_S,
    GEV2_TO_CM2,
    GEV_IN_G,
)


@jit(nopython=True, cache=True)
def sigma_T_weighted(m_chi, m_phi, alpha, v_km_s, w_even, w_odd):
    r"""σ_T/m [cm²/g] with configurable partial-wave weights.

    Majorana: w_even=1.0, w_odd=3.0.  Dirac: w_even=1.0, w_odd=1.0.
    Uses Numba-JIT phase shifts from v22_raw_scan.
    """
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi

    if kappa < 1e-15:
        return 0.0

    if kappa < 5:
        x_max, n_steps = 50.0, 4000
    elif kappa < 50:
        x_max, n_steps = 80.0, 8000
    else:
        x_max, n_steps = 100.0, 12000

    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)

    sigma_sum = 0.0
    peak_contrib = 0.0
    n_small = 0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, n_steps)
        weight = w_even if l % 2 == 0 else w_odd
        contrib = weight * (2 * l + 1) * math.sin(delta) ** 2
        sigma_sum += contrib
        if contrib > peak_contrib:
            peak_contrib = contrib
        if peak_contrib > 0.0 and contrib / peak_contrib < 1e-4:
            n_small += 1
            if n_small >= 5:
                break
        else:
            n_small = 0

    sigma_gev2 = 2.0 * math.pi * sigma_sum / (k * k)
    sigma_cm2 = sigma_gev2 * GEV2_TO_CM2
    return sigma_cm2 / (m_chi * GEV_IN_G)


def get_phase_shifts(m_chi, m_phi, alpha, v_km_s, l_max_cap=60):
    """Compute individual phase shifts δ_ℓ using Numba-accelerated solver."""
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi

    if kappa < 1e-15:
        return [], kappa, lam

    if kappa < 5:
        x_max, n_steps = 50.0, 4000
    elif kappa < 50:
        x_max, n_steps = 80.0, 8000
    else:
        x_max, n_steps = 100.0, 12000

    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)),
                l_max_cap)

    deltas = []
    for l in range(l_max + 1):
        d = vpm_phase_shift(l, kappa, lam, x_max, n_steps)
        deltas.append(d)

    return deltas, kappa, lam
