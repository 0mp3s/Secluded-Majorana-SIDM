#!/usr/bin/env python3
"""
vev_scan.py — Effective Potential & VEV from CP Violation  (V2 — CORRECTED)
============================================================================
BUGFIX V2: The V1 parameterization forced y_s² + y_p² = 4πα (a circle),
which is WRONG.  The physics requires only α_s α_p = α²/8 (a hyperbola).
α_s and α_p are independent; y_s = √(4πα_s), y_p = √(4πα_p).

The natural point α_s = α (SIDM unchanged), α_p = α/8 lies OFF the old
circle.  This V2 scans along the relic-matched hyperbola.

Derivation chain:
─────────────────
  Lagrangian:   L ⊃ -½ χ̄(y_s + iy_p γ⁵)χ φ + ½(∂φ)² - V(φ)

  Tree-level:   V₀(φ) = ½ m²_φ φ² + (μ₃/6) φ³ + (λ₄/24) φ⁴

  Coleman-Weinberg 1-loop (Majorana χ, n_f = 2 real dof, μ_R = m_χ):
    V_CW(φ) = -(1/32π²) M⁴_eff(φ) [ln(M²_eff(φ)/m²_χ) - 3/2]

    M²_eff(φ) = (m_χ + y_s φ/2)² + (y_p φ/2)²

  Renormalization:
    δV(φ) ≡ V_CW(φ) - V_CW(0) - V'_CW(0)·φ - ½V''_CW(0)·φ²

  Effective:    V_eff(φ) = V₀(φ) + δV(φ)

  Relic constraint (mixed Majorana s-wave):
    ⟨σv⟩₀ = 2π α_s α_p / m_χ²
    Matching to pure-scalar formula: α_s α_p = α²/8

  Parameterization (CORRECTED — V2):
    Scan variable:  r = α_s / α   (r = 1 means SIDM unchanged)
    Then:           α_s = r α,   α_p = α/(8r)
    Couplings:      y_s = √(4π r α),   y_p = √(4π α/(8r)) = √(π α/(2r))
    VPM parameter:  λ = α_s m_χ/m_φ = r × λ_orig

  Note: r = 1 → α_s = α, α_p = α/8   (preprint natural point)
        r = 1/√8 → α_s = α_p = α/√8   (CP-symmetric point)

Scan dimensions:
  - 5 benchmark points (BP1, BP9, BP16, MAP, MAP_relic)
  - r = α_s/α ∈ [0.01, 10]  (log-spaced, 30 points)
  - μ₃/m_φ ∈ [1.7, 80]      (cannibal lower bound = 1.7)
  - λ₄ ∈ [10⁻⁴, 4π]        (perturbativity upper bound)

Outputs:
  CSV:      dark_energy_exploration/vev_scan_results.csv
  Terminal: per-BP summary + global summary
"""

import sys
import os
import time
import csv
import math
import numpy as np
from pathlib import Path
from scipy.optimize import minimize_scalar

# ── path setup ────────────────────────────────────────────────
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "core"))
from global_config import GC
from run_logger import RunLogger

# ── output ────────────────────────────────────────────────────
OUT_DIR = Path(__file__).resolve().parent
OUT_CSV = OUT_DIR / "vev_scan_results.csv"

# ── observed vacuum energy density ────────────────────────────
RHO_LAMBDA_GEV4 = 2.58e-47  # GeV⁴ (Planck 2018: Ω_Λ h² = 0.3153)

# ── scan grid ─────────────────────────────────────────────────
N_RATIO = 30        # α_s/α ratio points
N_MU3   = 12
N_LAM4  = 12

# r = α_s/α.  Range covers:
#   r = 0.01 → α_s = 0.01α (very weak SIDM)
#   r = 1/√8 ≈ 0.354 → CP-symmetric (α_s = α_p)
#   r = 1 → natural point (α_s = α, α_p = α/8)
#   r = 10 → strong SIDM (α_p very small)
RATIOS        = np.geomspace(0.01, 10.0, N_RATIO)
MU3_RATIOS    = np.geomspace(1.7, 80.0, N_MU3)    # μ₃ / m_φ
LAM4_VALUES   = np.geomspace(1e-4, 4 * np.pi, N_LAM4)

R_CP_SYMMETRIC = 1.0 / math.sqrt(8)  # ≈ 0.354
R_NATURAL      = 1.0                  # α_s = α


# ══════════════════════════════════════════════════════════════
#  Physics: Effective Potential
# ══════════════════════════════════════════════════════════════

def M2_eff(phi, m_chi, y_s, y_p):
    """
    Effective Majorana mass² in classical background field φ.

    From L_int = -½ χ̄(y_s + iy_p γ⁵)χ φ, in background ⟨φ⟩ = φ₀:
      scalar mass:       m_s(φ) = m_χ + y_s φ/2
      pseudoscalar mass: m_p(φ) = y_p φ/2
      M²_eff = m_s² + m_p²
    """
    return (m_chi + y_s * phi / 2.0)**2 + (y_p * phi / 2.0)**2


def V_tree(phi, m_phi, mu3, lam4):
    """
    Tree-level scalar potential in GeV⁴.
    V₀(φ) = ½ m²_φ φ² + (μ₃/3!) φ³ + (λ₄/4!) φ⁴
    """
    return (0.5 * m_phi**2 * phi**2
            + (mu3 / 6.0) * phi**3
            + (lam4 / 24.0) * phi**4)


def _V_CW_raw(phi, m_chi, y_s, y_p):
    """
    Raw Coleman-Weinberg one-loop correction from Majorana χ.

    n_f = 2 (Majorana real dof),  renorm scale μ_R = m_χ.

    V_CW(φ) = -(n_f / 64π²) M⁴_eff(φ) [ln(M²_eff / μ²_R) - 3/2]
            = -(1 / 32π²)   M⁴_eff(φ) [ln(M²_eff / m²_χ) - 3/2]
    """
    M2 = M2_eff(phi, m_chi, y_s, y_p)
    if M2 < 1e-300:
        return 0.0
    return -(1.0 / (32.0 * np.pi**2)) * M2**2 * (np.log(M2 / m_chi**2) - 1.5)


def build_V_eff(m_chi, m_phi, y_s, y_p, mu3, lam4):
    """
    Build a fast V_eff(φ) function with pre-computed CW subtraction.

    Returns (V_vec, V_scalar):
      V_vec(phi_array) — vectorized, for coarse scans
      V_scalar(phi)    — scalar, for minimization
    """
    # ── Pre-compute CW renormalization constants at φ = 0 ──
    #    δV(φ) = V_CW(φ) - c₀ - c₁φ - ½c₂φ²
    c0 = _V_CW_raw(0.0, m_chi, y_s, y_p)
    eps = 1e-7 * m_chi
    cwp = _V_CW_raw(+eps, m_chi, y_s, y_p)
    cwm = _V_CW_raw(-eps, m_chi, y_s, y_p)
    c1 = (cwp - cwm) / (2.0 * eps)           # V'_CW(0)
    c2 = (cwp - 2.0 * c0 + cwm) / eps**2     # V''_CW(0)

    mu_R2 = m_chi**2
    inv_32pi2 = 1.0 / (32.0 * np.pi**2)

    def V_vec(phi_arr):
        """Vectorized effective potential for numpy arrays."""
        # Tree
        vt = (0.5 * m_phi**2 * phi_arr**2
              + (mu3 / 6.0) * phi_arr**3
              + (lam4 / 24.0) * phi_arr**4)
        # CW
        M2 = (m_chi + y_s * phi_arr / 2.0)**2 + (y_p * phi_arr / 2.0)**2
        M2 = np.maximum(M2, 1e-300)
        vcw = -inv_32pi2 * M2**2 * (np.log(M2 / mu_R2) - 1.5)
        # Renormalized CW
        dv = vcw - c0 - c1 * phi_arr - 0.5 * c2 * phi_arr**2
        return vt + dv

    def V_scalar(phi):
        """Scalar version for scipy minimize."""
        M2 = (m_chi + y_s * phi / 2.0)**2 + (y_p * phi / 2.0)**2
        if M2 < 1e-300:
            vcw = 0.0
        else:
            vcw = -inv_32pi2 * M2**2 * (np.log(M2 / mu_R2) - 1.5)
        vt = (0.5 * m_phi**2 * phi**2
              + (mu3 / 6.0) * phi**3
              + (lam4 / 24.0) * phi**4)
        dv = vcw - c0 - c1 * phi - 0.5 * c2 * phi**2
        return vt + dv

    return V_vec, V_scalar


# ══════════════════════════════════════════════════════════════
#  Minimization
# ══════════════════════════════════════════════════════════════

def find_minimum(V_vec, V_scalar, phi_range):
    """
    Find global minimum of V_eff in [-phi_range, +phi_range].
    Returns (phi_min, V_min, is_stable).
    """
    # Coarse vectorized scan
    phis = np.linspace(-phi_range, phi_range, 2000)
    Vs = V_vec(phis)

    # Collect local-minimum candidates
    candidates = [(0.0, V_scalar(0.0))]
    for i in range(1, len(Vs) - 1):
        if Vs[i] < Vs[i - 1] and Vs[i] < Vs[i + 1]:
            lo = phis[max(0, i - 3)]
            hi = phis[min(len(phis) - 1, i + 3)]
            res = minimize_scalar(V_scalar, bounds=(lo, hi),
                                  method='bounded',
                                  options={'xatol': 1e-14})
            candidates.append((res.x, res.fun))

    # Deepest minimum
    candidates.sort(key=lambda c: c[1])
    phi_min, V_min = candidates[0]

    # Stability: V''(φ_min) > 0
    eps = max(1e-8 * abs(phi_min), 1e-12)
    d2V = (V_scalar(phi_min + eps) - 2.0 * V_min
           + V_scalar(phi_min - eps)) / eps**2

    return phi_min, V_min, d2V > 0


# ══════════════════════════════════════════════════════════════
#  Main scan
# ══════════════════════════════════════════════════════════════

CSV_FIELDS = [
    "BP", "m_chi_GeV", "m_phi_GeV", "alpha",
    "r_alpha_s_over_alpha", "alpha_s", "alpha_p",
    "y_s", "y_p", "alpha_s_x_alpha_p",
    "lambda_VPM",
    "sidm_regime",
    "mu3_over_mphi", "lambda4",
    "phi_min_GeV", "V_eff_min_GeV4", "DeltaV_GeV4",
    "log10_abs_DeltaV", "log10_ratio_to_Lambda",
    "nontrivial", "stable",
    "delta_mchi_frac",
    "sigma_v_s_wave_GeV2",
]


def sidm_regime(lam_vpm):
    if lam_vpm < 0.3:
        return "Born"
    elif lam_vpm < 3.0:
        return "single-res"
    elif lam_vpm < 50.0:
        return "multi-res"
    else:
        return "classical"


def main():
    t_start = time.time()

    bps = GC.all_benchmarks()
    n_total = len(bps) * N_RATIO * N_MU3 * N_LAM4

    sep = "=" * 78
    print(sep)
    print("  V_eff SCAN V2 — CORRECTED: Hyperbolic Parameterization")
    print(sep)
    print()
    print("  BUGFIX: V1 forced y_s² + y_p² = 4πα (circle) — WRONG.")
    print("  V2 uses the correct constraint: α_s α_p = α²/8 (hyperbola).")
    print("  Scan variable: r = α_s/α,  then α_p = α/(8r)")
    print()
    print("  Derivation:")
    print("    L_int = -½ χ̄(y_s + iy_p γ⁵)χ φ")
    print("    ⟨σv⟩₀ = 2π α_s α_p / m_χ²  (mixed Majorana s-wave)")
    print("    Match to V10: 2π α_s α_p / m² = π α² / (4 m²)  →  α_s α_p = α²/8")
    print("    V₀(φ) = ½m²_φ φ² + (μ₃/6)φ³ + (λ₄/24)φ⁴")
    print("    V_CW(φ) = -(1/32π²) M⁴_eff [ln(M²_eff/m²_χ) - 3/2]")
    print("    M²_eff(φ) = (m_χ + y_s φ/2)² + (y_p φ/2)²")
    print()
    print(f"  Scan: {len(bps)} BPs × {N_RATIO} r × {N_MU3} μ₃ × {N_LAM4} λ₄"
          f" = {n_total:,} points")
    print(f"  r range: [{RATIOS[0]:.3f}, {RATIOS[-1]:.1f}]  "
          f"(r=1 = natural point, r={R_CP_SYMMETRIC:.3f} = CP-symmetric)")
    print(f"  ρ_Λ(obs) = {RHO_LAMBDA_GEV4:.2e} GeV⁴")
    print()

    fcsv = open(OUT_CSV, "w", newline="", encoding="utf-8")
    writer = csv.DictWriter(fcsv, fieldnames=CSV_FIELDS)
    writer.writeheader()

    n_done = 0
    n_nontrivial = 0
    closest_row = None
    closest_log_ratio = 999.0

    for bp in bps:
        label = bp["label"]
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3   # → GeV
        alpha = bp["alpha"]
        lam_orig = alpha * m_chi / m_phi

        print(f"  ── {label}:  m_χ = {m_chi} GeV,  m_φ = {m_phi*1e3:.2f} MeV,"
              f"  α = {alpha:.3e},  λ_orig = {lam_orig:.2f} ──")

        bp_nontrivial = 0
        bp_best_ratio = 999.0
        bp_best_info = ""

        for r in RATIOS:
            # Correct parameterization: hyperbola α_s α_p = α²/8
            a_s = r * alpha
            a_p = alpha**2 / (8.0 * a_s)     # = α/(8r)
            y_s = math.sqrt(4.0 * math.pi * a_s)
            y_p = math.sqrt(4.0 * math.pi * a_p)

            # VPM coupling
            lam_vpm = a_s * m_chi / m_phi
            regime = sidm_regime(lam_vpm)

            # s-wave (constant along hyperbola = α²/8 × 2π/m²)
            sv0 = 2.0 * math.pi * a_s * a_p / m_chi**2

            for mu3r in MU3_RATIOS:
                mu3 = mu3r * m_phi

                for lam4 in LAM4_VALUES:
                    V_vec, V_scl = build_V_eff(m_chi, m_phi, y_s, y_p,
                                               mu3, lam4)

                    phi_range = max(100.0, 5.0 * m_chi)
                    phi_min, V_min, stable = find_minimum(V_vec, V_scl,
                                                          phi_range)

                    V_at_0 = V_scl(0.0)
                    DeltaV = V_min - V_at_0
                    nontrivial = abs(phi_min) > 1e-10

                    if nontrivial:
                        bp_nontrivial += 1
                        n_nontrivial += 1

                    abs_dv = abs(DeltaV) if DeltaV != 0.0 else 1e-300
                    log_dv = math.log10(abs_dv)
                    ratio_val = abs_dv / RHO_LAMBDA_GEV4
                    log_ratio = math.log10(ratio_val) if ratio_val > 0 else -999.0

                    dm_frac = abs(y_s * phi_min / 2.0) / m_chi if nontrivial else 0.0

                    row = {
                        "BP":                    label,
                        "m_chi_GeV":             m_chi,
                        "m_phi_GeV":             f"{m_phi:.6e}",
                        "alpha":                 f"{alpha:.4e}",
                        "r_alpha_s_over_alpha":  f"{r:.6f}",
                        "alpha_s":               f"{a_s:.6e}",
                        "alpha_p":               f"{a_p:.6e}",
                        "y_s":                   f"{y_s:.6e}",
                        "y_p":                   f"{y_p:.6e}",
                        "alpha_s_x_alpha_p":     f"{a_s * a_p:.6e}",
                        "lambda_VPM":            f"{lam_vpm:.4f}",
                        "sidm_regime":           regime,
                        "mu3_over_mphi":         f"{mu3r:.2f}",
                        "lambda4":               f"{lam4:.6e}",
                        "phi_min_GeV":           f"{phi_min:.6e}",
                        "V_eff_min_GeV4":        f"{V_min:.6e}",
                        "DeltaV_GeV4":           f"{DeltaV:.6e}",
                        "log10_abs_DeltaV":      f"{log_dv:.2f}",
                        "log10_ratio_to_Lambda": f"{log_ratio:.2f}",
                        "nontrivial":            nontrivial,
                        "stable":                stable,
                        "delta_mchi_frac":       f"{dm_frac:.6e}",
                        "sigma_v_s_wave_GeV2":   f"{sv0:.6e}",
                    }
                    writer.writerow(row)
                    n_done += 1

                    if nontrivial and abs(log_ratio) < abs(closest_log_ratio):
                        closest_log_ratio = log_ratio
                        closest_row = dict(row)

                    if nontrivial and abs(log_ratio) < abs(bp_best_ratio):
                        bp_best_ratio = log_ratio
                        bp_best_info = (
                            f"r={r:.3f}  λ_VPM={lam_vpm:.2f} ({regime})  "
                            f"μ₃/m_φ={mu3r:.1f}  λ₄={lam4:.2e}  "
                            f"⟨φ⟩={phi_min:.3e} GeV  "
                            f"Δm_χ/m_χ={dm_frac:.2e}"
                        )

        bp_total = N_RATIO * N_MU3 * N_LAM4
        print(f"     Non-trivial VEV: {bp_nontrivial}/{bp_total}"
              f"  ({100.0 * bp_nontrivial / bp_total:.1f}%)")
        if bp_nontrivial > 0:
            print(f"     Closest to Λ:  log₁₀(|ΔV|/ρ_Λ) = {bp_best_ratio:.1f}")
            print(f"       {bp_best_info}")
        print()

    fcsv.close()
    elapsed = time.time() - t_start

    # ══════════════════════════════════════════════════════════
    #  Global Summary
    # ══════════════════════════════════════════════════════════
    print(sep)
    print("  SUMMARY")
    print(sep)
    print(f"  Total scanned:    {n_done:,}")
    print(f"  Non-trivial VEV:  {n_nontrivial:,}"
          f"  ({100.0 * n_nontrivial / n_done:.1f}%)")
    print(f"  CSV written:      {OUT_CSV}")
    print(f"  Elapsed:          {elapsed:.1f}s")
    print()

    if closest_row:
        print("  ── Closest to observed ρ_Λ ──")
        print(f"     BP            = {closest_row['BP']}")
        print(f"     r = α_s/α     = {closest_row['r_alpha_s_over_alpha']}")
        print(f"     α_s           = {closest_row['alpha_s']}")
        print(f"     α_p           = {closest_row['alpha_p']}")
        print(f"     λ_VPM         = {closest_row['lambda_VPM']}  ({closest_row['sidm_regime']})")
        print(f"     μ₃/m_φ        = {closest_row['mu3_over_mphi']}")
        print(f"     λ₄            = {closest_row['lambda4']}")
        print(f"     ⟨φ⟩           = {closest_row['phi_min_GeV']} GeV")
        print(f"     ΔV            = {closest_row['DeltaV_GeV4']} GeV⁴")
        print(f"     log₁₀(|ΔV|/ρ_Λ) = {closest_log_ratio:.2f}")
        print(f"     Δm_χ/m_χ      = {closest_row['delta_mchi_frac']}")
        print()

        orders = abs(closest_log_ratio)
        if orders < 5:
            print(f"     >>> Within 5 orders of ρ_Λ — promising!")
        elif orders < 20:
            print(f"     >>> {orders:.0f} orders from ρ_Λ — interesting gap reduction")
        else:
            print(f"     >>> {orders:.0f} orders from ρ_Λ — CC problem persists")

    print()
    print("  ── Physics Notes (V2 — CORRECTED) ──")
    print()
    print("  1. PARAMETERIZATION FIX:")
    print("     V1 (WRONG): y_s = y cos δ, y_p = y sin δ  →  α_s + α_p = α")
    print("     V2 (CORRECT): α_s α_p = α²/8 (hyperbola, independent couplings)")
    print("     r = α_s/α is the free parameter; α_p = α/(8r)")
    print()
    print("  2. SPECIAL POINTS along the relic hyperbola:")
    print(f"     r = 1       : natural point (α_s = α, α_p = α/8)")
    print(f"                   SIDM unchanged from pipeline, y_p suppressed")
    print(f"     r = {R_CP_SYMMETRIC:.3f}  : CP-symmetric (α_s = α_p = α/√8)")
    print(f"                   maximal CP violation at fixed relic density")
    print(f"     r → 0       : SIDM dies (Born regime)")
    print(f"     r → ∞       : α_p → 0 → s-wave vanishes (p-wave only)")
    print()
    print("  3. CW LOOP: V_CW depends on y_s² + y_p² = 4π(α_s + α_p).")
    print("     On the hyperbola α_s α_p = const, the sum α_s + α_p has")
    print("     a MINIMUM at α_s = α_p (AM-GM), and grows for r ≫ 1 or r ≪ 1.")
    print("     → ΔV depends on r!  (Unlike V1 which wrongly showed flat ΔV.)")
    print()
    print("  4. KEY CHECK: Δm_χ/m_χ ≪ 1 ensures the existing SIDM pipeline")
    print("     remains valid when φ develops a VEV.")
    print()
    print("  5. CONSISTENCY: ⟨σv⟩₀ = 2π α_s α_p / m_χ² = π α²/(4 m_χ²)")
    print("     is CONSTANT along the entire hyperbola.  Relic density is")
    print("     identical for all r — only SIDM and vacuum structure change.")
    print(sep)

    return 0


if __name__ == "__main__":
    with RunLogger(
        "dark_energy_exploration/vev_scan.py",
        stage="9 - Dark Energy Exploration",
        params={
            "N_ratio": N_RATIO, "N_mu3": N_MU3, "N_lam4": N_LAM4,
            "r_range": f"[{RATIOS[0]:.3f}, {RATIOS[-1]:.1f}]",
            "mu3_range": f"[1.7, 80] × m_phi",
            "lam4_range": f"[1e-4, 4π]",
            "parameterization": "hyperbola: alpha_s * alpha_p = alpha^2/8",
        },
        data_source="global_config.json (benchmark points)",
    ) as rl:
        status = main()
        rl.add_output(str(OUT_CSV))
        rl.set_notes("VEV scan V2: CORRECTED parameterization (hyperbola)")
    sys.exit(status)
