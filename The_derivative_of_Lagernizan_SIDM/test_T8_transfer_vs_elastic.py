#!/usr/bin/env python3
r"""
test_T8_transfer_vs_elastic.py
===============================
T8: Compare transfer and elastic (viscosity) cross sections for Majorana vs Dirac.

PHYSICS:
  The total elastic cross section used in SIDM is:
    σ_T = (2π/k²) Σ_ℓ w_ℓ (2ℓ+1) sin²δ_ℓ

  The transfer (momentum-transfer) cross section is:
    σ_tr = (4π/k²) Σ_ℓ w_ℓ (ℓ+1) sin²(δ_ℓ - δ_{ℓ+1})

  The viscosity cross section is:
    σ_vis = (4π/k²) Σ_ℓ w_ℓ (ℓ+1)(ℓ+2)/(2ℓ+3) sin²(δ_ℓ - δ_{ℓ+2})

  For Majorana: w_ℓ = (1, 3, 1, 3, ...) — from T4
  For Dirac:    all w_ℓ = 1 (distinguishable particles)

  KEY QUESTION: Does the ratio σ_tr/σ_T depend on fermion type?
  If σ_tr/σ_T is different for Majorana vs Dirac at the same velocities,
  this provides an OBSERVATIONAL SIGNATURE distinguishing the two.

  SIDM observations typically constrain σ_T/m, but N-body simulations
  sometimes use σ_tr/m. The mapping between them depends on the
  particle type — this test quantifies that dependence.

PREDICTION:
  If the ratio R_tr(v) = σ_tr/σ_T is velocity-dependent AND differs
  between Majorana and Dirac, any SIDM observation interpreted with
  the wrong fermion type gets the wrong σ/m ↔ v mapping.

Output: data/archive/T8_transfer_vs_elastic_*.csv
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
# VPM PHASE SHIFTS (shared code from T4)
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
# CROSS SECTION FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def compute_all_cross_sections(m_chi, m_phi, alpha, v_km_s,
                                w_even, w_odd):
    """
    Compute σ_T, σ_tr, σ_vis simultaneously (share phase shift calculation).

    Returns dict with σ/m [cm²/g] for each type.
    """
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi

    if kappa < 1e-15:
        return {"sigma_T": 0.0, "sigma_tr": 0.0, "sigma_vis": 0.0}

    if kappa < 5:
        x_max, n_steps = 50.0, 4000
    elif kappa < 50:
        x_max, n_steps = 80.0, 8000
    else:
        x_max, n_steps = 100.0, 12000

    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)

    # Compute ALL phase shifts first
    deltas = []
    for l in range(l_max + 2):  # +2 for transfer/viscosity formulas
        delta = vpm_phase_shift(l, kappa, lam, x_max, n_steps)
        deltas.append(delta)

    # σ_T = (2π/k²) Σ_ℓ w_ℓ (2ℓ+1) sin²δ_ℓ
    sum_T = 0.0
    for l in range(l_max + 1):
        w = w_even if l % 2 == 0 else w_odd
        sum_T += w * (2*l + 1) * math.sin(deltas[l])**2

    # σ_tr = (4π/k²) Σ_ℓ w_ℓ (ℓ+1) sin²(δ_ℓ - δ_{ℓ+1})
    sum_tr = 0.0
    for l in range(l_max):
        w = w_even if l % 2 == 0 else w_odd
        sum_tr += w * (l + 1) * math.sin(deltas[l] - deltas[l+1])**2

    # σ_vis = (4π/k²) Σ_ℓ w_ℓ (ℓ+1)(ℓ+2)/(2ℓ+3) sin²(δ_ℓ - δ_{ℓ+2})
    sum_vis = 0.0
    for l in range(l_max - 1):
        w = w_even if l % 2 == 0 else w_odd
        sum_vis += w * (l+1)*(l+2) / (2*l+3) * math.sin(deltas[l] - deltas[l+2])**2

    # Convert to cm²/g
    factor = GEV2_TO_CM2 / (m_chi * GEV_IN_G)

    sigma_T   = 2.0 * math.pi * sum_T   / (k * k) * factor
    sigma_tr  = 4.0 * math.pi * sum_tr  / (k * k) * factor
    sigma_vis = 4.0 * math.pi * sum_vis / (k * k) * factor

    return {"sigma_T": sigma_T, "sigma_tr": sigma_tr, "sigma_vis": sigma_vis}


# ═══════════════════════════════════════════════════════════════════════
# MAIN TEST
# ═══════════════════════════════════════════════════════════════════════

def main():
    bp_labels = _LOCAL_CFG["test_benchmarks"]
    velocities = [10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0, 500.0, 1000.0]

    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  T8: TRANSFER vs ELASTIC CROSS SECTIONS — Majorana vs Dirac    ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()
    print("Definitions:")
    print("  σ_T   = (2π/k²) Σ w_ℓ(2ℓ+1) sin²δ_ℓ           [elastic/total]")
    print("  σ_tr  = (4π/k²) Σ w_ℓ(ℓ+1) sin²(δ_ℓ-δ_{ℓ+1})  [momentum transfer]")
    print("  σ_vis = (4π/k²) Σ w_ℓ(ℓ+1)(ℓ+2)/(2ℓ+3) sin²(δ_ℓ-δ_{ℓ+2}) [viscosity]")
    print()
    print("  Majorana: w = (1,3)    Dirac: w = (1,1)")
    print()

    all_rows = []

    for bp_label in bp_labels:
        bp = GC.benchmark(bp_label)
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3
        alpha = bp["alpha"]

        print(f"{'━' * 70}")
        print(f"  {bp_label}: m_χ={m_chi} GeV, m_φ={bp['m_phi_MeV']} MeV, α={alpha:.4e}")
        print(f"{'━' * 70}")
        print()

        header = (f"  {'v':>6} │ {'σ_T^M':>8} {'σ_tr^M':>8} {'σ_vis^M':>8} "
                  f"{'R_tr^M':>7} │ {'σ_T^D':>8} {'σ_tr^D':>8} {'σ_vis^D':>8} "
                  f"{'R_tr^D':>7} │ {'ΔR_tr%':>7}")
        print(header)
        print("  " + "─" * (len(header) - 2))

        for v in velocities:
            # Majorana w=(1,3)
            cs_maj = compute_all_cross_sections(m_chi, m_phi, alpha, v, 1.0, 3.0)
            # Dirac w=(1,1)
            cs_dir = compute_all_cross_sections(m_chi, m_phi, alpha, v, 1.0, 1.0)

            # Ratios σ_tr/σ_T
            R_tr_maj = cs_maj["sigma_tr"] / cs_maj["sigma_T"] if cs_maj["sigma_T"] > 0 else float('nan')
            R_tr_dir = cs_dir["sigma_tr"] / cs_dir["sigma_T"] if cs_dir["sigma_T"] > 0 else float('nan')
            R_vis_maj = cs_maj["sigma_vis"] / cs_maj["sigma_T"] if cs_maj["sigma_T"] > 0 else float('nan')
            R_vis_dir = cs_dir["sigma_vis"] / cs_dir["sigma_T"] if cs_dir["sigma_T"] > 0 else float('nan')

            # Difference in ratio
            delta_R_tr = (R_tr_maj - R_tr_dir) / R_tr_dir * 100 if (
                math.isfinite(R_tr_dir) and R_tr_dir != 0) else float('nan')

            print(f"  {v:>6.0f} │ {cs_maj['sigma_T']:>8.3f} {cs_maj['sigma_tr']:>8.3f} "
                  f"{cs_maj['sigma_vis']:>8.3f} {R_tr_maj:>7.3f} │ "
                  f"{cs_dir['sigma_T']:>8.3f} {cs_dir['sigma_tr']:>8.3f} "
                  f"{cs_dir['sigma_vis']:>8.3f} {R_tr_dir:>7.3f} │ "
                  f"{delta_R_tr:>+7.1f}%")

            all_rows.append({
                "bp":               bp_label,
                "m_chi_GeV":        m_chi,
                "m_phi_MeV":        bp["m_phi_MeV"],
                "alpha":            alpha,
                "v_km_s":           v,
                "sigma_T_maj":      cs_maj["sigma_T"],
                "sigma_tr_maj":     cs_maj["sigma_tr"],
                "sigma_vis_maj":    cs_maj["sigma_vis"],
                "sigma_T_dir":      cs_dir["sigma_T"],
                "sigma_tr_dir":     cs_dir["sigma_tr"],
                "sigma_vis_dir":    cs_dir["sigma_vis"],
                "R_tr_maj":         R_tr_maj,
                "R_tr_dir":         R_tr_dir,
                "R_vis_maj":        R_vis_maj,
                "R_vis_dir":        R_vis_dir,
                "delta_R_tr_pct":   delta_R_tr,
            })

        print()

    # ── Summary ───────────────────────────────────────────────────────
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()

    # Average R_tr by velocity
    for v in velocities:
        v_rows = [r for r in all_rows if r["v_km_s"] == v]
        avg_R_maj = np.mean([r["R_tr_maj"] for r in v_rows if math.isfinite(r["R_tr_maj"])])
        avg_R_dir = np.mean([r["R_tr_dir"] for r in v_rows if math.isfinite(r["R_tr_dir"])])
        diff = (avg_R_maj - avg_R_dir) / avg_R_dir * 100 if avg_R_dir != 0 else 0
        print(f"  v={v:>6.0f} km/s: ⟨R_tr^Maj⟩={avg_R_maj:.4f}  "
              f"⟨R_tr^Dir⟩={avg_R_dir:.4f}  Δ={diff:+.1f}%")

    print()
    print("CONCLUSION:")
    print()

    # Check if the ratio R_tr differs significantly
    delta_values = [r["delta_R_tr_pct"] for r in all_rows
                    if math.isfinite(r["delta_R_tr_pct"])]

    if delta_values:
        max_delta = max(abs(d) for d in delta_values)
        rms_delta = np.sqrt(np.mean(np.array(delta_values)**2))

        if max_delta > 5.0:
            print(f"  ✓ SIGNIFICANT DIFFERENCE: σ_tr/σ_T differs by up to {max_delta:.1f}%")
            print(f"    between Majorana and Dirac at the same (m_χ, m_φ, α).")
            print()
            print("  PHYSICS IMPLICATION:")
            print("    SIDM constraints from N-body simulations that use σ_tr")
            print("    map to DIFFERENT σ_T depending on whether DM is Majorana or Dirac.")
            print("    This systematic bias has been largely overlooked in the literature.")
            print()
            print("  OBSERVATIONAL TEST:")
            print("    If a cluster (v~1000 km/s) and dwarf (v~30 km/s) both have σ/m ~ 1,")
            print("    the R_tr(v) shape distinguishes Majorana from Dirac statistics.")
        elif max_delta > 1.0:
            print(f"  ~ MARGINAL DIFFERENCE: σ_tr/σ_T differs by up to {max_delta:.1f}%")
            print(f"    Effect is small but may be measurable in precision SIDM analyses.")
        else:
            print(f"  ✗ NO SIGNIFICANT DIFFERENCE: max |Δ| = {max_delta:.2f}%")
            print(f"    σ_tr/σ_T is approximately independent of fermion type.")

        print(f"\n  RMS deviation: {rms_delta:.2f}%")

    print("=" * 70)
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha", "v_km_s",
        "sigma_T_maj", "sigma_tr_maj", "sigma_vis_maj",
        "sigma_T_dir", "sigma_tr_dir", "sigma_vis_dir",
        "R_tr_maj", "R_tr_dir", "R_vis_maj", "R_vis_dir",
        "delta_R_tr_pct",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T8_transfer_vs_elastic.py",
        stage="T8 — Transfer vs Elastic",
        params={"benchmarks": bp_labels, "n_velocities": len(velocities)},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T8_transfer_vs_elastic",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_rows)

        rl.add_output(str(out_csv))
        n_sig = sum(1 for d in delta_values if abs(d) > 5.0) if delta_values else 0
        rl.set_notes(f"{len(all_rows)} measurements, "
                     f"{n_sig} with |Δ|>5%. "
                     f"RMS deviation = {rms_delta:.2f}%" if delta_values else "no data")
        print(f"Output: {out_csv}")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T8_transfer_vs_elastic"):
        main()
