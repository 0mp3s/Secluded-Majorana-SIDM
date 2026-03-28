#!/usr/bin/env python3
r"""
test_T9_majorana_vs_dirac.py
==============================
T9: Compute R(v) = σ_T^Maj / σ_T^Dir as a function of velocity.

THE KEY QUESTION:
  Given the SAME Lagrangian parameters (m_χ, m_φ, α), does the
  Majorana vs Dirac nature of the dark matter leave an observable
  imprint on the velocity-dependent cross section?

PHYSICS:
  Majorana (χ = χ^c):
    σ_T = (2π/k²) Σ_ℓ w_ℓ^M (2ℓ+1) sin²δ_ℓ
    w_ℓ^M = 1 (even), 3 (odd) — from identical-particle + spin statistics

  Dirac (χ ≠ χ̄):
    σ_T = (2π/k²) Σ_ℓ (2ℓ+1) sin²δ_ℓ
    All w_ℓ^D = 1 — distinguishable particles, no exchange

  The ratio R(v) = σ_T^Maj / σ_T^Dir encodes the pure Majorana effect.

  PREDICTION:
    - In Born limit (low v, few partial waves): R ≈ 2 (average weight = 2)
    - In classical limit (high v, many ℓ): R → 2 (statistical average)
    - In resonance regime: R can fluctuate — specific ℓ values dominate
    - If R(v) is NON-MONOTONIC, specific velocity windows discriminate better

  OBSERVATIONAL SIGNATURE:
    If we observe σ/m at multiple velocities (dwarfs ~30, groups ~300,
    clusters ~1000 km/s), the velocity dependence σ(v) has different
    SHAPE for Majorana vs Dirac, even with the same underlying potential.

Output: data/archive/T9_majorana_vs_dirac_*.csv
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

# ── local config ──────────────────────────────────────────────────────
_LOCAL_CFG_PATH = _HERE / "data" / "config.json"
with open(_LOCAL_CFG_PATH, "r", encoding="utf-8") as _fh:
    _LOCAL_CFG = json.load(_fh)

_PC = GC.physical_constants()
GEV2_TO_CM2 = _PC["GEV2_to_cm2"]
GEV_IN_G    = _PC["GeV_in_g"]
C_KM_S      = _PC["c_km_s"]


# ═══════════════════════════════════════════════════════════════════════
# VPM PHASE SHIFTS (same independent implementation)
# ═══════════════════════════════════════════════════════════════════════

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


# ═══════════════════════════════════════════════════════════════════════
# σ_T FOR ARBITRARY WEIGHTS
# ═══════════════════════════════════════════════════════════════════════

def sigma_T_weighted(m_chi, m_phi, alpha, v_km_s, w_even, w_odd):
    """σ_T/m [cm²/g] with given weights."""
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


def compute_ratio(m_chi, m_phi, alpha, v_km_s):
    """Compute R(v) = σ_T^Maj / σ_T^Dir and both cross sections."""
    sigma_maj = sigma_T_weighted(m_chi, m_phi, alpha, v_km_s, 1.0, 3.0)
    sigma_dir = sigma_T_weighted(m_chi, m_phi, alpha, v_km_s, 1.0, 1.0)

    if sigma_dir > 0:
        R = sigma_maj / sigma_dir
    else:
        R = float('nan')

    return sigma_maj, sigma_dir, R


# ═══════════════════════════════════════════════════════════════════════
# BORN LIMIT ANALYTIC CHECK
# ═══════════════════════════════════════════════════════════════════════

def born_ratio_analytic():
    """
    In the Born limit (weak coupling), σ ∝ α² and partial waves are small.
    The σ_T is dominated by ℓ=0 (s-wave).

    For ℓ=0 only: σ^Maj/σ^Dir = w_0/1 = 1/1 = 1.

    But once ℓ=1 contributes: the triplet weight 3 kicks in.
    In the full Born limit with all ℓ:
      σ^Maj = (2π/k²) [1·1·sin²δ₀ + 3·3·sin²δ₁ + 1·5·sin²δ₂ + ...]
      σ^Dir = (2π/k²) [1·1·sin²δ₀ + 1·3·sin²δ₁ + 1·5·sin²δ₂ + ...]

    With Born δ_ℓ ∝ small:
      R → [Σ w_ℓ(2ℓ+1)δ_ℓ²] / [Σ (2ℓ+1)δ_ℓ²]

    This is velocity-dependent through which ℓ values contribute!
    """
    print("  Born limit analysis:")
    print("    R = [Σ w_ℓ(2ℓ+1)sin²δ_ℓ] / [Σ (2ℓ+1)sin²δ_ℓ]")
    print("    If only ℓ=0 contributes: R = 1")
    print("    If all ℓ contribute equally: R → Σw_ℓ(2ℓ+1)/Σ(2ℓ+1) → 2")
    print("    Intermediate: 1 < R < 2+ε (resonances can push higher)")
    print()


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    bp_labels = _LOCAL_CFG["test_benchmarks"]

    # Dense velocity grid to capture resonance structure
    velocities = np.concatenate([
        np.arange(5, 50, 5),        # 5, 10, ..., 45
        np.arange(50, 200, 10),      # 50, 60, ..., 190
        np.arange(200, 500, 25),     # 200, 225, ..., 475
        np.arange(500, 1100, 50),    # 500, 550, ..., 1050
    ]).tolist()

    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  T9: MAJORANA vs DIRAC — R(v) = σ_T^Maj / σ_T^Dir            ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()

    born_ratio_analytic()

    all_rows = []

    for bp_label in bp_labels:
        bp = GC.benchmark(bp_label)
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3
        alpha = bp["alpha"]
        lam = alpha * m_chi / m_phi

        print(f"{'━' * 70}")
        print(f"  {bp_label}: m_χ={m_chi} GeV, m_φ={bp['m_phi_MeV']} MeV, "
              f"α={alpha:.4e}, λ={lam:.1f}")
        print(f"{'━' * 70}")
        print()

        header = (f"  {'v [km/s]':>10} {'σ/m Maj':>10} {'σ/m Dir':>10} "
                  f"{'R(v)':>8} {'regime':>12}")
        print(header)
        print("  " + "─" * 55)

        Rs = []
        for v in velocities:
            sigma_maj, sigma_dir, R = compute_ratio(m_chi, m_phi, alpha, v)

            # Classify regime
            kappa = (m_chi / 2 * v / C_KM_S) / m_phi
            if kappa < 0.5:
                regime = "s-wave"
            elif kappa < 5:
                regime = "few-ℓ"
            elif kappa < 50:
                regime = "semi-class"
            else:
                regime = "classical"

            Rs.append(R)

            # Only print a subset but store all
            if v in [10, 30, 50, 100, 200, 300, 500, 1000]:
                print(f"  {v:>10.0f} {sigma_maj:>10.4f} {sigma_dir:>10.4f} "
                      f"{R:>8.4f} {regime:>12}")

            all_rows.append({
                "bp":           bp_label,
                "m_chi_GeV":    m_chi,
                "m_phi_MeV":    bp["m_phi_MeV"],
                "alpha":        alpha,
                "lambda":       lam,
                "v_km_s":       v,
                "sigma_m_maj":  sigma_maj,
                "sigma_m_dir":  sigma_dir,
                "R_v":          R,
                "kappa":        (m_chi / 2 * v / C_KM_S) / m_phi,
                "regime":       regime,
            })

        # Statistics for this BP
        R_arr = np.array([r for r in Rs if math.isfinite(r)])
        if len(R_arr) > 0:
            print(f"\n  R(v) statistics:")
            print(f"    min  = {R_arr.min():.4f}")
            print(f"    max  = {R_arr.max():.4f}")
            print(f"    mean = {R_arr.mean():.4f}")
            print(f"    std  = {R_arr.std():.4f}")
            print(f"    monotonic? {_is_monotonic(Rs)}")
        print()

    # ── Global analysis ───────────────────────────────────────────────
    print("=" * 70)
    print("GLOBAL ANALYSIS: R(v) = σ_T^Maj / σ_T^Dir")
    print("=" * 70)
    print()

    # Group by velocity
    print("  Average R(v) across all BPs:")
    print(f"  {'v [km/s]':>10} {'⟨R⟩':>8} {'min':>8} {'max':>8} {'spread':>8}")
    print("  " + "─" * 50)

    for v in [10, 30, 50, 100, 200, 300, 500, 1000]:
        v_rows = [r for r in all_rows if r["v_km_s"] == v]
        if v_rows:
            Rs = [r["R_v"] for r in v_rows if math.isfinite(r["R_v"])]
            if Rs:
                R_arr = np.array(Rs)
                print(f"  {v:>10.0f} {R_arr.mean():>8.4f} {R_arr.min():>8.4f} "
                      f"{R_arr.max():>8.4f} {R_arr.max()-R_arr.min():>8.4f}")

    print()

    # Observational window analysis
    print("  OBSERVATIONAL WINDOWS:")
    print()

    for bp_label in bp_labels:
        bp_rows = [r for r in all_rows if r["bp"] == bp_label]
        R_at_30 = next((r["R_v"] for r in bp_rows if r["v_km_s"] == 30), None)
        R_at_300 = next((r["R_v"] for r in bp_rows if r["v_km_s"] == 300), None)
        R_at_1000 = next((r["R_v"] for r in bp_rows if r["v_km_s"] == 1000), None)

        if all(x is not None for x in [R_at_30, R_at_300, R_at_1000]):
            print(f"  {bp_label}:")
            print(f"    Dwarf (30 km/s):    R = {R_at_30:.4f}")
            print(f"    Group (300 km/s):   R = {R_at_300:.4f}")
            print(f"    Cluster (1000 km/s): R = {R_at_1000:.4f}")
            spread = max(R_at_30, R_at_300, R_at_1000) - min(R_at_30, R_at_300, R_at_1000)
            print(f"    Max spread: {spread:.4f}")
            if spread > 0.1:
                print(f"    → VELOCITY-DEPENDENT Majorana imprint detected!")
            print()

    # ── Conclusion ────────────────────────────────────────────────────
    print("CONCLUSION:")
    print()

    all_R = np.array([r["R_v"] for r in all_rows if math.isfinite(r["R_v"])])
    if len(all_R) > 0:
        overall_spread = all_R.max() - all_R.min()
        overall_mean = all_R.mean()

        print(f"  Overall R(v) range: [{all_R.min():.4f}, {all_R.max():.4f}]")
        print(f"  Overall mean: {overall_mean:.4f}")
        print()

        if overall_spread > 0.3:
            print("  ✓ R(v) is STRONGLY velocity-dependent.")
            print("    The Majorana nature of DM changes σ_T(v) shape significantly.")
            print("    Multi-velocity observations CAN distinguish Majorana from Dirac.")
        elif overall_spread > 0.1:
            print("  ~ R(v) is MODERATELY velocity-dependent.")
            print("    The Majorana signal is present but requires precision measurements.")
        else:
            print("  ✗ R(v) is approximately CONSTANT ≈ {:.2f}.".format(overall_mean))
            print("    The Majorana vs Dirac distinction is a simple rescaling of σ,")
            print("    absorbed into a redefinition of α. Not independently measurable.")

        if abs(overall_mean - 2.0) < 0.3:
            print()
            print(f"  R ≈ 2: consistent with average weight ⟨w⟩ = (1+3)/2 = 2.")
            print(f"  The identical-particle exchange effectively doubles σ_T.")

    print("=" * 70)
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha", "lambda", "v_km_s",
        "sigma_m_maj", "sigma_m_dir", "R_v", "kappa", "regime",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T9_majorana_vs_dirac.py",
        stage="T9 — Majorana vs Dirac R(v)",
        params={"benchmarks": bp_labels, "n_velocities": len(velocities)},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T9_majorana_vs_dirac",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_rows)

        rl.add_output(str(out_csv))

        if len(all_R) > 0:
            rl.set_notes(f"R(v) range=[{all_R.min():.3f}, {all_R.max():.3f}], "
                         f"mean={all_R.mean():.3f}. "
                         f"{len(bp_labels)} BPs × {len(velocities)} velocities.")
        print(f"Output: {out_csv}")


def _is_monotonic(values):
    """Check if a list is monotonically increasing or decreasing."""
    clean = [v for v in values if math.isfinite(v)]
    if len(clean) < 3:
        return "N/A"
    diffs = np.diff(clean)
    if np.all(diffs >= -1e-6):
        return "increasing"
    elif np.all(diffs <= 1e-6):
        return "decreasing"
    else:
        return "NON-MONOTONIC ← interesting!"


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T9_majorana_vs_dirac"):
        main()
