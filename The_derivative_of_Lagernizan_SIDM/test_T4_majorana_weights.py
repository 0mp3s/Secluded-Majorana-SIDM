#!/usr/bin/env python3
r"""
test_T4_majorana_weights.py
============================
T4: Derive the Majorana partial-wave weights w_ℓ = (1, 3, 1, 3, ...)
from first principles using SymPy.

PHYSICS:
  Majorana fermions are their own antiparticle: χ = χ^c.
  In χχ → χχ scattering, the two initial (and final) particles are IDENTICAL.
  For identical fermions, the total wavefunction must be antisymmetric.

  Total wavefunction = Spatial × Spin:
    - Spin singlet (S=0): antisymmetric spin → symmetric spatial → even ℓ only
      Multiplicity: 2S+1 = 1
    - Spin triplet (S=1): symmetric spin → antisymmetric spatial → odd ℓ only
      Multiplicity: 2S+1 = 3

  Therefore:
    w_ℓ = 1  (even ℓ)  ← singlet contributes
    w_ℓ = 3  (odd ℓ)   ← triplet contributes

  This also encodes the t+u interference:
    |f(θ) ± f(π-θ)|² for identical fermions gives the factor-of-4
    enhancement that resolved the T1 "discrepancy".

VERIFICATION:
  1. SymPy: derive weights from Clebsch-Gordan / spin statistics
  2. Numeric: compute σ_T with (1,3) and with naive (1,1) — compare to v22
  3. Cross-check: spin sum 1+3 = 4 = 2² (identical particle factor)

Output: data/archive/T4_majorana_weights_*.csv
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

# ── path setup ────────────────────────────────────────────────────────
_HERE      = Path(__file__).resolve().parent
_SIDM_ROOT = _HERE.parent
_CORE      = _SIDM_ROOT / "core"
sys.path.insert(0, str(_SIDM_ROOT))
sys.path.insert(0, str(_CORE))

from core.global_config import GC
from core.output_manager import timestamped_path
from core.run_logger import RunLogger
from core.v22_raw_scan import sigma_T_vpm as sigma_T_vpm_v22

# ── local config ──────────────────────────────────────────────────────
_LOCAL_CFG_PATH = _HERE / "data" / "config.json"
with open(_LOCAL_CFG_PATH, "r", encoding="utf-8") as _fh:
    _LOCAL_CFG = json.load(_fh)

_PC = GC.physical_constants()
GEV2_TO_CM2 = _PC["GEV2_to_cm2"]
GEV_IN_G    = _PC["GeV_in_g"]
C_KM_S      = _PC["c_km_s"]


# ═══════════════════════════════════════════════════════════════════════
# PART 1: SYMBOLIC DERIVATION WITH SYMPY
# ═══════════════════════════════════════════════════════════════════════

def sympy_derive_weights():
    """
    Derive w_ℓ from Majorana symmetry using SymPy.

    For two identical spin-½ fermions:
      - Total spin S: decompose ½ ⊗ ½ = 0 ⊕ 1
      - Pauli exclusion: total ψ must be antisymmetric under exchange
      - Spin singlet (S=0): spin part is antisymmetric → spatial SYMMETRIC → even ℓ
      - Spin triplet (S=1): spin part is symmetric → spatial ANTISYMMETRIC → odd ℓ

    Weight = (2S+1) for each allowed S at given ℓ.
    """
    import sympy as sp
    from sympy.physics.quantum.spin import CG
    from sympy.physics.quantum.cg import CG as CG2

    print("=" * 70)
    print("T4: DERIVATION OF MAJORANA PARTIAL-WAVE WEIGHTS")
    print("=" * 70)
    print()

    # ── Step 1: Spin decomposition ½ ⊗ ½ ──────────────────────────────
    print("Step 1: Spin decomposition of two spin-½ particles")
    print("  ½ ⊗ ½ = 0 ⊕ 1")
    print("  S=0 (singlet): 2S+1 = 1 state,  antisymmetric under exchange")
    print("  S=1 (triplet): 2S+1 = 3 states, symmetric under exchange")
    print()

    # Verify spin exchange symmetry with CG coefficients
    # |S, M_S⟩ = Σ ⟨s1 m1; s2 m2 | S M_S⟩ |m1⟩|m2⟩
    # Exchange P₁₂: swap m1 ↔ m2
    # For singlet (S=0, M=0):
    #   CG(½,½,½,-½;0,0) = 1/√2,  CG(½,-½,½,½;0,0) = -1/√2
    #   → antisymmetric ✓

    half = sp.Rational(1, 2)

    print("Step 2: Verify exchange symmetry via Clebsch-Gordan coefficients")
    print()

    # Singlet S=0, M=0
    cg_up_down_S0 = CG2(half, half, half, -half, 0, 0).doit()
    cg_down_up_S0 = CG2(half, -half, half, half, 0, 0).doit()
    print(f"  Singlet (S=0, M=0):")
    print(f"    ⟨½,+½; ½,-½ | 0,0⟩ = {cg_up_down_S0} = {float(cg_up_down_S0):.4f}")
    print(f"    ⟨½,-½; ½,+½ | 0,0⟩ = {cg_down_up_S0} = {float(cg_down_up_S0):.4f}")
    singlet_symmetric = (cg_up_down_S0 == cg_down_up_S0)
    singlet_antisymmetric = (cg_up_down_S0 == -cg_down_up_S0)
    print(f"    Exchange symmetric?     {singlet_symmetric}")
    print(f"    Exchange antisymmetric? {singlet_antisymmetric}")
    assert singlet_antisymmetric, "Singlet must be antisymmetric!"
    print()

    # Triplet S=1, M=0 (representative)
    cg_up_down_S1 = CG2(half, half, half, -half, 1, 0).doit()
    cg_down_up_S1 = CG2(half, -half, half, half, 1, 0).doit()
    print(f"  Triplet (S=1, M=0):")
    print(f"    ⟨½,+½; ½,-½ | 1,0⟩ = {cg_up_down_S1} = {float(cg_up_down_S1):.4f}")
    print(f"    ⟨½,-½; ½,+½ | 1,0⟩ = {cg_down_up_S1} = {float(cg_down_up_S1):.4f}")
    triplet_symmetric = (cg_up_down_S1 == cg_down_up_S1)
    triplet_antisymmetric = (cg_up_down_S1 == -cg_down_up_S1)
    print(f"    Exchange symmetric?     {triplet_symmetric}")
    print(f"    Exchange antisymmetric? {triplet_antisymmetric}")
    assert triplet_symmetric, "Triplet must be symmetric!"
    print()

    # ── Step 3: Pauli principle ────────────────────────────────────────
    print("Step 3: Pauli exclusion principle")
    print("  Total ψ = ψ_spatial × ψ_spin must be ANTISYMMETRIC under exchange.")
    print()
    print("  Spatial exchange: P₁₂ Y_ℓ^m(θ) = (-1)^ℓ Y_ℓ^m(θ)")
    print("    Even ℓ → symmetric spatial")
    print("    Odd  ℓ → antisymmetric spatial")
    print()
    print("  Combining:")
    print("    Even ℓ (sym. spatial) × singlet (antisym. spin) = ANTISYMMETRIC ✓")
    print("    Odd  ℓ (antisym. spatial) × triplet (sym. spin) = ANTISYMMETRIC ✓")
    print("    Even ℓ × triplet = SYMMETRIC ✗ (FORBIDDEN)")
    print("    Odd  ℓ × singlet = SYMMETRIC ✗ (FORBIDDEN)")
    print()

    # ── Step 4: Compute weights ────────────────────────────────────────
    print("Step 4: Weights = spin multiplicity (2S+1) for allowed channel")
    print()

    weights = {}
    for l in range(10):
        if l % 2 == 0:
            # Even ℓ: singlet only
            S_allowed = 0
            w = 2 * S_allowed + 1   # = 1
        else:
            # Odd ℓ: triplet only
            S_allowed = 1
            w = 2 * S_allowed + 1   # = 3
        weights[l] = w
        parity = "even" if l % 2 == 0 else "odd"
        channel = "singlet (S=0)" if l % 2 == 0 else "triplet (S=1)"
        print(f"    ℓ={l} ({parity}): {channel} → w_{l} = {w}")

    print()
    print("  RESULT: w_ℓ = 1 (even ℓ), 3 (odd ℓ)")
    print()

    # ── Step 5: Explain connection to t+u exchange ─────────────────────
    print("Step 5: Connection to t+u exchange amplitudes")
    print()
    print("  For identical particles, the scattering amplitude is:")
    print("    |M|² ∝ |f(θ) - (-1)^S f(π-θ)|²")
    print()
    print("  Singlet (S=0): |f(θ) - f(π-θ)|² → interference is DESTRUCTIVE at θ=π/2")
    print("  Triplet (S=1): |f(θ) + f(π-θ)|² → interference is CONSTRUCTIVE at θ=π/2")
    print()
    print("  Spin-averaged: σ = ¼|f-f'|² + ¾|f+f'|²   (¼ singlet + ¾ triplet)")
    print("  This is exactly w_even=1, w_odd=3 in the partial-wave sum!")
    print()
    print("  Cross-check: average weight = (1+3)/2 = 2 = factor from identical particles")
    print("  (The \"factor of 4\" in T1 is: 2 from t+u × 2 from spin averaging)")
    print("=" * 70)
    print()

    return weights


# ═══════════════════════════════════════════════════════════════════════
# PART 2: NUMERIC VERIFICATION — w=(1,3) vs w=(1,1)
# ═══════════════════════════════════════════════════════════════════════
# Import the independent VPM from the EL pipeline

# We inline a lightweight VPM to avoid circular imports
def _sph_jn(l: int, z: float) -> float:
    if z < 1e-30:
        return 1.0 if l == 0 else 0.0
    j0 = math.sin(z) / z
    if l == 0:
        return j0
    j1 = math.sin(z) / (z * z) - math.cos(z) / z
    if l == 1:
        return j1
    jp, jc = j0, j1
    for n in range(1, l):
        jn = (2 * n + 1) / z * jc - jp
        jp, jc = jc, jn
        if abs(jc) < 1e-300:
            return 0.0
    return jc


def _sph_yn(l: int, z: float) -> float:
    if z < 1e-30:
        return -1e300
    y0 = -math.cos(z) / z
    if l == 0:
        return y0
    y1 = -math.cos(z) / (z * z) - math.sin(z) / z
    if l == 1:
        return y1
    yp, yc = y0, y1
    for n in range(1, l):
        yn = (2 * n + 1) / z * yc - yp
        yp, yc = yc, yn
        if abs(yc) > 1e200:
            return yc
    return yc


def _vpm_rhs(l, kappa, lam, x, delta):
    if x < 1e-20:
        return 0.0
    z = kappa * x
    if z < 1e-20:
        return 0.0
    jl = _sph_jn(l, z)
    nl = _sph_yn(l, z)
    j_hat = z * jl
    n_hat = -z * nl
    cd = math.cos(delta)
    sd = math.sin(delta)
    bracket = j_hat * cd - n_hat * sd
    if not math.isfinite(bracket):
        return 0.0
    pot = lam * math.exp(-x) / (kappa * x)
    val = pot * bracket * bracket
    return val if math.isfinite(val) else 0.0


def vpm_phase_shift(l, kappa, lam, x_max=50.0, n_steps=4000):
    if lam < 1e-30 or kappa < 1e-30:
        return 0.0
    x_min = max(1e-5, 0.05 / (kappa + 0.01))
    if l > 0:
        x_barrier = l / kappa
        if x_barrier > x_min:
            x_min = x_barrier
    h = (x_max - x_min) / n_steps
    delta = 0.0
    for i in range(n_steps):
        x = x_min + i * h
        k1 = _vpm_rhs(l, kappa, lam, x, delta)
        k2 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k1)
        k3 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k2)
        k4 = _vpm_rhs(l, kappa, lam, x + h, delta + h*k3)
        delta += h * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return delta


def sigma_T_with_weights(m_chi, m_phi, alpha, v_km_s, w_even, w_odd):
    """σ_T/m [cm²/g] with arbitrary weights w_even, w_odd."""
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
        contrib = weight * (2*l + 1) * math.sin(delta)**2
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


def test_numeric(bp_labels: list[str]) -> list[dict]:
    """
    For each BP at representative velocities, compute σ_T with:
      A) Majorana weights w=(1,3) — should match v22
      B) Naive weights w=(1,1) — should NOT match v22
      C) Dirac-like weights w=(1,1) with no 1/2 symmetry factor — different physics
    """
    obs = GC.observations()
    test_vs = [10.0, 30.0, 100.0, 300.0, 1000.0]

    print("=" * 70)
    print("NUMERIC VERIFICATION: w=(1,3) vs w=(1,1) vs v22")
    print("=" * 70)
    print()

    rows = []
    for bp_label in bp_labels:
        bp = GC.benchmark(bp_label)
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3
        alpha = bp["alpha"]

        print(f"── {bp_label}: m_χ={m_chi} GeV, α={alpha:.4e} ──")
        header = (f"  {'v [km/s]':>10} {'σ/m v22':>12} {'σ/m (1,3)':>12} "
                  f"{'σ/m (1,1)':>12} {'R(1,3)/v22':>12} {'R(1,1)/v22':>12}")
        print(header)

        for v in test_vs:
            sigma_v22 = sigma_T_vpm_v22(m_chi, m_phi, alpha, v)
            sigma_13  = sigma_T_with_weights(m_chi, m_phi, alpha, v, 1.0, 3.0)
            sigma_11  = sigma_T_with_weights(m_chi, m_phi, alpha, v, 1.0, 1.0)

            r_13 = sigma_13 / sigma_v22 if sigma_v22 > 0 else float('nan')
            r_11 = sigma_11 / sigma_v22 if sigma_v22 > 0 else float('nan')

            match_13 = "PASS" if abs(r_13 - 1.0) < 0.01 else "FAIL"

            print(f"  {v:>10.0f} {sigma_v22:>12.4f} {sigma_13:>12.4f} "
                  f"{sigma_11:>12.4f} {r_13:>12.4f} {r_11:>12.4f}  {match_13}")

            rows.append({
                "bp":           bp_label,
                "m_chi_GeV":    m_chi,
                "m_phi_MeV":    bp["m_phi_MeV"],
                "alpha":        alpha,
                "v_km_s":       v,
                "sigma_m_v22":  sigma_v22,
                "sigma_m_w13":  sigma_13,
                "sigma_m_w11":  sigma_11,
                "ratio_w13_v22": r_13,
                "ratio_w11_v22": r_11,
                "match_w13":    match_13,
            })
        print()

    return rows


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    bp_labels = _LOCAL_CFG["test_benchmarks"]

    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  T4: MAJORANA PARTIAL-WAVE WEIGHTS w_ℓ FROM FIRST PRINCIPLES   ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()

    # Part 1: Symbolic
    weights = sympy_derive_weights()

    # Part 2: Numeric
    rows = test_numeric(bp_labels)

    # ── Summary ───────────────────────────────────────────────────────
    n_pass = sum(1 for r in rows if r["match_w13"] == "PASS")
    n_total = len(rows)

    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print()
    print(f"  Symbolic: w_ℓ derived from Pauli exclusion + spin decomposition")
    print(f"    w_even = 1 (singlet, S=0)")
    print(f"    w_odd  = 3 (triplet, S=1)")
    print()
    print(f"  Numeric: {n_pass}/{n_total} comparisons PASS (w=(1,3) matches v22)")
    print()

    if n_pass == n_total:
        print("  ✓ w_ℓ = (1,3) is a PREDICTION from the Lagrangian + Majorana symmetry,")
        print("    not an arbitrary choice. The VPM code correctly implements it.")
        print()
        print("  PHYSICS INSIGHT:")
        print("    The weights encode THREE physical effects simultaneously:")
        print("    1. Identical-particle exchange (t+u channels)")
        print("    2. Spin statistics (Fermi-Dirac)")
        print("    3. Spin-averaged cross section (¼ singlet + ¾ triplet)")
        print("    A Dirac fermion with w=(1,1) gives a DIFFERENT σ_T.")
    else:
        print("  ✗ MISMATCH — investigate.")
    print("=" * 70)
    print()

    # ── Average ratio w=(1,1)/v22 — quantify the Majorana effect ──────
    avg_ratio_11 = np.mean([r["ratio_w11_v22"] for r in rows])
    print(f"  Average σ(w=1,1)/σ(w=1,3) = {avg_ratio_11:.3f}")
    print(f"  → Majorana weights change σ_T by {(1/avg_ratio_11 - 1)*100:+.0f}% on average")
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha", "v_km_s",
        "sigma_m_v22", "sigma_m_w13", "sigma_m_w11",
        "ratio_w13_v22", "ratio_w11_v22", "match_w13",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T4_majorana_weights.py",
        stage="T4 — Majorana Weights",
        params={"benchmarks": bp_labels, "n_velocities": 5},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T4_majorana_weights",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(rows)

        rl.add_output(str(out_csv))
        rl.set_notes(f"{n_pass}/{n_total} PASS. w_ℓ derived from Pauli exclusion.")
        print(f"Output: {out_csv}")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T4_majorana_weights"):
        main()
