#!/usr/bin/env python3
r"""
test_T6_annihilation_relic.py
=============================
T6: Resolve the α convention by testing annihilation → relic density.

THE KEY INSIGHT:
  SIDM scattering (VPM) uses α as the coupling in V(r) = -α e^{-mr}/r.
  Both conventions (A: α=y²/16π, B: α=y²/8π) give the SAME σ_T
  because w_ℓ compensates.

  But annihilation ⟨σv⟩ = πα²/(4m²) depends on α SQUARED —
  so the relic density Ωh² discriminates between conventions.

TEST LOGIC:
  For each BP, compute Ωh² under THREE conventions:
    1. α_ann = α_config                    (α=y²/4π, as v27 assumes)
    2. α_ann = α_config/4                  (if α_config=y²/16π → y²=16πα → α_std=4α_config)
    3. α_ann = α_config/2                  (if α_config=y²/8π  → y²=8πα  → α_std=2α_config)

  Actually, the question is: what goes into ⟨σv⟩ = πα²/(4m²)?
  The v27 code uses: sv0 = π·α²/(4m²) with α = α_config directly.
  
  If α_config is the "effective" coupling that already absorbed vertex factors,
  then we need to know: does ⟨σv⟩ = π·α_config²/(4m²) give correct Ωh²?
  Or should it be ⟨σv⟩ = π·(Nα_config)²/(4m²) for some factor N?

  The ANSWER: Ωh² ≈ 0.12 pins the correct convention. Only one factor N works.

Output: data/archive/T6_annihilation_relic_*.csv
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
from core.v27_boltzmann_relic import (
    solve_boltzmann, Y_to_omega_h2, kolb_turner_swave, sigma_v_swave,
)

# ── local config ──────────────────────────────────────────────────────
_LOCAL_CFG_PATH = _HERE / "data" / "config.json"
with open(_LOCAL_CFG_PATH, "r", encoding="utf-8") as _fh:
    _LOCAL_CFG = json.load(_fh)

# ── constants ─────────────────────────────────────────────────────────
_CC = GC.cosmological_constants()
OMEGA_CDM_TARGET = _CC["omega_h2_target"]   # 0.120 ± 0.001


def test_one_bp(bp_label: str) -> list[dict]:
    """Test three α conventions for one BP via relic density."""
    bp = GC.benchmark(bp_label)
    m_chi   = bp["m_chi_GeV"]
    m_phi   = bp["m_phi_MeV"] * 1e-3
    alpha_c = bp["alpha"]

    # Three hypotheses for what α_config means:
    #   H1: α_config = y²/(4π)   →  ⟨σv⟩ uses α_config directly  (v27 assumption)
    #   H2: α_config = y²/(16π)  →  y² = 16π·α_c  →  α_std = 4·α_c  →  ⟨σv⟩ uses 4·α_c
    #   H3: α_config = y²/(8π)   →  y² = 8π·α_c   →  α_std = 2·α_c  →  ⟨σv⟩ uses 2·α_c
    #
    # In all cases: ⟨σv⟩ = π·α_ann²/(4m²), where α_ann = y²/(4π).
    # The question is what α_ann equals in terms of α_config.

    hypotheses = [
        {"name": "H1: α=y²/(4π)",   "label": "y²/(4π)",   "alpha_ann_factor": 1},
        {"name": "H2: α=y²/(16π)",  "label": "y²/(16π)",  "alpha_ann_factor": 4},
        {"name": "H3: α=y²/(8π)",   "label": "y²/(8π)",   "alpha_ann_factor": 2},
    ]

    rows = []
    for hyp in hypotheses:
        factor = hyp["alpha_ann_factor"]
        # Under this hypothesis, α_std = y²/(4π) = factor × α_config
        alpha_ann = factor * alpha_c

        # s-wave annihilation
        sv0 = math.pi * alpha_ann**2 / (4.0 * m_chi**2)  # GeV⁻²

        # Boltzmann solve
        x_arr, Y_arr = solve_boltzmann(m_chi, sv0, x_start=1.0,
                                        x_end=500.0, n_steps=20000)
        Y_inf = Y_arr[-1]
        omega_h2 = Y_to_omega_h2(Y_inf, m_chi)

        # Kolb-Turner analytic
        x_fo_kt, Y_inf_kt = kolb_turner_swave(m_chi, sv0)
        omega_kt = Y_to_omega_h2(Y_inf_kt, m_chi)

        # Compare to target
        delta_pct = (omega_h2 - OMEGA_CDM_TARGET) / OMEGA_CDM_TARGET * 100

        rows.append({
            "bp":               bp_label,
            "m_chi_GeV":        m_chi,
            "m_phi_MeV":        bp["m_phi_MeV"],
            "alpha_config":     alpha_c,
            "hypothesis":       hyp["label"],
            "factor":           factor,
            "alpha_annihilation": alpha_ann,
            "sv0_GeV2":         sv0,
            "Y_inf":            Y_inf,
            "omega_h2_boltzmann": omega_h2,
            "omega_h2_KT":      omega_kt,
            "omega_target":     OMEGA_CDM_TARGET,
            "delta_pct":        delta_pct,
        })

    return rows


def main():
    bp_labels = _LOCAL_CFG["test_benchmarks"]

    print("=" * 75)
    print("T6: ANNIHILATION → RELIC DENSITY — CONVENTION DISCRIMINATOR")
    print("=" * 75)
    print()
    print(f"Target: Ωh² = {OMEGA_CDM_TARGET}")
    print(f"If only one convention gives Ωh² ≈ 0.12, it's the correct one.")
    print()

    all_rows = []

    for bp_label in bp_labels:
        bp = GC.benchmark(bp_label)
        print(f"── {bp_label}: m_χ={bp['m_chi_GeV']} GeV, "
              f"α_config={bp['alpha']:.4e} ──")

        rows = test_one_bp(bp_label)
        all_rows.extend(rows)

        header = (f"  {'Hypothesis':<16} {'factor':>6} {'α_ann':>12} "
                  f"{'⟨σv⟩ [GeV⁻²]':>14} {'Ωh²':>10} {'Δ%':>8}")
        print(header)
        print("  " + "-" * (len(header) - 2))

        for r in rows:
            marker = " ◄" if abs(r["delta_pct"]) < 50 else ""
            print(f"  {r['hypothesis']:<16} {r['factor']:>6} "
                  f"{r['alpha_annihilation']:>12.4e} "
                  f"{r['sv0_GeV2']:>14.4e} "
                  f"{r['omega_h2_boltzmann']:>10.4f} "
                  f"{r['delta_pct']:>+8.1f}%{marker}")
        print()

    # ── Conclusion ────────────────────────────────────────────────────
    print("=" * 75)
    print("CONCLUSION")
    print("=" * 75)

    for hyp_label in ["y²/(4π)", "y²/(16π)", "y²/(8π)"]:
        hyp_rows = [r for r in all_rows if r["hypothesis"] == hyp_label]
        avg_delta = np.mean([abs(r["delta_pct"]) for r in hyp_rows])
        close = all(abs(r["delta_pct"]) < 100 for r in hyp_rows)
        status = "✓ CONSISTENT" if close else "✗ RULED OUT"
        print(f"  {hyp_label:<12}: avg |Δ%| = {avg_delta:>8.1f}%  → {status}")

    # Find winner
    best = min(
        ["y²/(4π)", "y²/(16π)", "y²/(8π)"],
        key=lambda h: np.mean([abs(r["delta_pct"])
                               for r in all_rows if r["hypothesis"] == h])
    )
    print()
    print(f"  ► WINNER: α_config = {best}")
    print(f"    This means the v27 Boltzmann solver and global_config.json")
    print(f"    use the convention α = {best}.")
    print("=" * 75)
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha_config", "hypothesis",
        "factor", "alpha_annihilation", "sv0_GeV2", "Y_inf",
        "omega_h2_boltzmann", "omega_h2_KT", "omega_target", "delta_pct",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T6_annihilation_relic.py",
        stage="T6 — Annihilation Convention",
        params={"benchmarks": bp_labels, "omega_target": OMEGA_CDM_TARGET},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T6_annihilation_relic",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_rows)

        rl.add_output(str(out_csv))
        rl.set_notes(f"Winner: α = {best}. Tested {len(bp_labels)} BPs × 3 conventions.")
        print(f"Output: {out_csv}")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T6_annihilation_relic"):
        main()
