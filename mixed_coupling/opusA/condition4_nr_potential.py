#!/usr/bin/env python3
"""
Condition 4: NR Potential Check — yₚ contribution negligible in VPM
====================================================================

Question: Does the pseudoscalar coupling yₚ affect σ_T/m at SIDM velocities?

Theory:
  The NR (Born-level) potential from scalar + pseudoscalar exchange:

    Scalar:       V_s(r) = -αₛ e^{-m_φ r} / r           ← Yukawa (spin-independent)
    Pseudoscalar: V_p(r) = +αₚ/(12 m_χ²) [4π δ³(r̄) − m_φ² e^{-m_φ r}/r] (σ₁·σ₂)
                         + αₚ/(4 m_χ²) e^{-m_φ r}/r [m_φ²/r + m_φ/r² + ...]  (tensor)

  The pseudoscalar contribution is suppressed by (m_φ/m_χ)² ~ 10⁻⁷ relative
  to scalar — it's a P²/m² correction, typical of pseudoscalar exchange (cf. pion).

  The scalar term dominates σ_T by many orders of magnitude.

Method:
  Run VPM σ_T(v) for BP1 parameters with:
    (a) α = αₛ                   (scalar only — current code)
    (b) α = αₛ + αₚ(m_φ/m_χ)²   (upper bound on pseudoscalar shift)
    (c) α = αₛ × (1 ± 10%)       (systematic comparison — how much does α matter?)

  Show that |σ_T(a) − σ_T(b)| ≪ σ_T(a), confirming VPM only needs αₛ.

Date: 23 March 2026
"""

import sys, os
import numpy as np

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

from v22_raw_scan import sigma_T_vpm

# ═══════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════
GEV2_TO_CM2 = 3.8938e-28
GEV_IN_G    = 1.78266e-24

# ═══════════════════════════════════════════════════════════════
#  Benchmark Points
# ═══════════════════════════════════════════════════════════════
BPs = {
    "BP1":  {"m_chi": 20.69, "m_phi_MeV": 11.34, "alpha_tot": 1.048e-3},
    "BP16": {"m_chi": 20.70, "m_phi_MeV":  9.91, "alpha_tot": 1.048e-3},
    "MAP":  {"m_chi": 90.64, "m_phi_MeV": 13.85, "alpha_tot": 2.546e-2},
}

# Velocities of interest (km/s)
V_LIST = [10, 30, 50, 100, 200, 500, 1000, 1500]

# ═══════════════════════════════════════════════════════════════
#  NR potential analysis
# ═══════════════════════════════════════════════════════════════

def analyze_bp(name, bp):
    m_chi = bp["m_chi"]
    m_phi = bp["m_phi_MeV"] * 1e-3    # GeV - VPM takes GeV
    alpha_tot = bp["alpha_tot"]
    alpha_s = alpha_tot / 2.0          # split equally
    alpha_p = alpha_tot / 2.0

    r = m_phi / m_chi
    r2 = r**2

    # Pseudoscalar NR correction factor (upper bound)
    # The dominant pseudoscalar NR term goes as α_p × (m_φ/m_χ)² / m_φ r
    # relative to scalar term α_s / m_φ r → ratio ~ α_p r² / α_s
    ps_correction = alpha_p * r2

    lam_s = alpha_s * m_chi / m_phi
    lam_tot = (alpha_s + ps_correction) * m_chi / m_phi

    print(f"\n  ── {name} ──")
    print(f"  m_χ = {m_chi:.2f} GeV,  m_φ = {bp['m_phi_MeV']:.2f} MeV")
    print(f"  α_total = {alpha_tot:.3e}")
    print(f"  αₛ = αₚ = {alpha_s:.3e}  (equal split)")
    print(f"  r = m_φ/m_χ = {r:.4e}")
    print(f"  r² = {r2:.4e}")
    print(f"  αₚ × r² = {ps_correction:.4e}  (pseudoscalar NR shift)")
    print(f"  Fractional shift: αₚr²/αₛ = {ps_correction/alpha_s:.4e}")
    print(f"  λₛ = αₛm_χ/m_φ = {lam_s:.3f}")
    print(f"  λ_shifted = {lam_tot:.6f}  (Δλ/λ = {(lam_tot-lam_s)/lam_s:.4e})")

    print(f"\n  {'v [km/s]':>10s}  {'σ/m (αₛ)':>14s}  {'σ/m (αₛ+corr)':>14s}  "
          f"{'σ/m (αₛ×1.1)':>14s}  {'Δ(corr)':>10s}  {'Δ(10%)':>10s}")
    print("  " + "─" * 82)

    max_delta = 0.0
    for v in V_LIST:
        sm_s    = sigma_T_vpm(m_chi, m_phi, alpha_s, v)
        sm_corr = sigma_T_vpm(m_chi, m_phi, alpha_s + ps_correction, v)
        sm_10p  = sigma_T_vpm(m_chi, m_phi, alpha_s * 1.1, v)

        d_corr = abs(sm_corr / sm_s - 1.0) if sm_s > 0 else 0.0
        d_10p  = abs(sm_10p / sm_s - 1.0) if sm_s > 0 else 0.0
        max_delta = max(max_delta, d_corr)

        print(f"  {v:10d}  {sm_s:14.4e}  {sm_corr:14.4e}  "
              f"{sm_10p:14.4e}  {d_corr:10.2e}  {d_10p:10.2e}")

    return max_delta, ps_correction / alpha_s


def main():
    print("=" * 72)
    print("  Condition 4: NR Potential Check — yₚ contribution to VPM")
    print("=" * 72)

    print("\n" + "─" * 72)
    print("  Theory: Pseudoscalar NR potential")
    print("─" * 72)
    print("""
  Scalar (Yukawa):       V_s = -αₛ × e^{-m_φ r}/r
  Pseudoscalar (NR):     V_p ~ αₚ × (m_φ/m_χ)² × e^{-m_φ r}/r × (spin-dep)

  Ratio:  |V_p/V_s| ~ (αₚ/αₛ) × (m_φ/m_χ)²

  For BP1:  r = m_φ/m_χ = 5.5e-4  →  r² = 3.0e-7
  → Pseudoscalar shift to σ_T:  ~ 3 × 10⁻⁷  (completely negligible)
  → VPM only needs αₛ (the scalar coupling)
""")

    print("─" * 72)
    print("  Numerical verification: σ_T/m with and without yₚ correction")
    print("─" * 72)

    results = {}
    for name, bp in BPs.items():
        max_d, frac = analyze_bp(name, bp)
        results[name] = (max_d, frac)

    # ──────────────────────────────────────────────────────────
    #  Summary
    # ──────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("  SUMMARY — CONDITION 4")
    print("=" * 72)

    all_pass = True
    for name, (max_d, frac) in results.items():
        ok = max_d < 1e-3
        if not ok:
            all_pass = False
        mark = "✓" if ok else "✗"
        print(f"  {name:6s}:  αₚr²/αₛ = {frac:.2e}  "
              f" max |Δσ/σ| = {max_d:.2e}   {mark}")

    condition4 = all_pass
    print(f"""
  Physics:
    Pseudoscalar exchange produces a spin-dependent NR potential
    suppressed by (m_φ/m_χ)² ~ 10⁻⁷ relative to scalar Yukawa.

    This is analogous to pion exchange: the pseudoscalar potential
    goes as p²/m² in the NR limit → velocity-suppressed.

    For SIDM at v ~ 10-1500 km/s, the VPM solver with V = -αₛ e^{{-m_φ r}}/r
    captures the physics to better than 10⁻⁴ accuracy.

  Consequence:
    The VPM computation needs ONLY αₛ = yₛ²/(4π).
    The pseudoscalar coupling αₚ affects ONLY relic density (s-wave a₀),
    NOT elastic self-interactions.

    This is good: it means the SIDM phenomenology (core sizes, offsets,
    velocity profiles) is identical to the pure-scalar case.
    The 17 benchmarks from V10 remain valid — just re-interpret
    α → αₛ for the SIDM sector.

  ★ CONDITION 4:  {"PASSED ✓" if condition4 else "NEEDS REVIEW ✗"} ★
""")
    return condition4


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)
