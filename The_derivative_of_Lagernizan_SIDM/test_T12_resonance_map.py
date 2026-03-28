#!/usr/bin/env python3
r"""
test_T12_resonance_map.py
==========================
T12: Map resonances in (λ, κ) space and study how they amplify R(v).

PHYSICS:
  The VPM phase shifts pass through δ_ℓ = nπ at specific values of
  the dimensionless coupling λ = αm_χ/m_φ and momentum κ = k/m_φ.
  Near these resonances, sin²δ_ℓ → 1 and the ℓ-th partial wave
  contributes maximally to σ_T.

  Since Majorana weights w_ℓ = (1, 3) differ from Dirac w_ℓ = 1,
  the ratio R(v) is most sensitive when specific ℓ values dominate.

  KEY QUESTION: Where in parameter space is R(v) maximally enhanced?
  Do our benchmarks sit near resonances?
  Can we find a "sweet spot" where R >> 2?

METHOD:
  1. For each BP, scan velocities densely (50 points)
  2. At each (BP, v), compute individual δ_ℓ values
  3. Identify near-resonance ℓ values (|sin²δ_ℓ| > 0.95)
  4. Compute R(v) and correlate with resonance structure
  5. Map the (λ, κ) plane for one BP with dense grid

Output: data/archive/T12_resonance_map_*.csv + output/T12_resonance_map.txt
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

from _vpm_fast import sigma_T_weighted, get_phase_shifts


def count_resonances(deltas):
    """Count near-resonance partial waves (sin²δ_ℓ > 0.95)."""
    n_res = 0
    res_ells = []
    for l, d in enumerate(deltas):
        sin2 = math.sin(d)**2
        if sin2 > 0.95:
            n_res += 1
            res_ells.append(l)
    return n_res, res_ells


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    bp_labels = _LOCAL_CFG["test_benchmarks"]
    bps = GC.benchmarks(*bp_labels)

    # Dense velocity grid
    velocities = np.logspace(np.log10(15), np.log10(2500), 40)

    print()
    print("╔═══════════════════════════════════════════════════════════════════╗")
    print("║  T12: RESONANCE MAP — Where do resonances amplify R(v)?        ║")
    print("╚═══════════════════════════════════════════════════════════════════╝")
    print()
    print(f"  Velocities: {len(velocities)} points from {velocities[0]:.0f} to {velocities[-1]:.0f} km/s")
    print(f"  Benchmarks: {bp_labels}")
    print()

    all_rows = []
    bp_summary = []

    for bp in bps:
        label = bp["label"]
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3
        alpha = bp["alpha"]
        lam_bp = alpha * m_chi / m_phi

        print(f"  ── {label} (λ={lam_bp:.2f}) ──")

        R_vals = []
        max_R = 0.0
        best_v = 0.0
        total_resonances = 0
        res_even_total = 0
        res_odd_total = 0

        for v in velocities:
            # Compute R(v)
            sigma_maj = sigma_T_weighted(m_chi, m_phi, alpha, v, 1.0, 3.0)
            sigma_dir = sigma_T_weighted(m_chi, m_phi, alpha, v, 1.0, 1.0)
            R = sigma_maj / sigma_dir if sigma_dir > 0 else float('nan')

            # Phase shift analysis
            deltas, kappa, lam = get_phase_shifts(m_chi, m_phi, alpha, v)
            n_res, res_ells = count_resonances(deltas)

            res_even = sum(1 for l in res_ells if l % 2 == 0)
            res_odd = sum(1 for l in res_ells if l % 2 == 1)

            # Dominant ℓ (highest sin²δ_ℓ contribution)
            if deltas:
                contribs = [(l, math.sin(d)**2 * (2*l+1)) for l, d in enumerate(deltas)]
                dom_l, dom_contrib = max(contribs, key=lambda x: x[1])
            else:
                dom_l, dom_contrib = -1, 0.0

            R_vals.append(R)
            if math.isfinite(R) and R > max_R:
                max_R = R
                best_v = v
            total_resonances += n_res
            res_even_total += res_even
            res_odd_total += res_odd

            all_rows.append({
                "bp": label,
                "m_chi_GeV": m_chi,
                "m_phi_MeV": bp["m_phi_MeV"],
                "alpha": alpha,
                "lambda": lam,
                "kappa": kappa,
                "v_km_s": v,
                "sigma_maj": sigma_maj,
                "sigma_dir": sigma_dir,
                "R": R,
                "n_resonances": n_res,
                "res_even": res_even,
                "res_odd": res_odd,
                "dominant_l": dom_l,
                "n_partialwaves": len(deltas),
            })

        R_arr = np.array([r for r in R_vals if math.isfinite(r)])
        print(f"    R(v) range: [{R_arr.min():.3f}, {R_arr.max():.3f}], mean={R_arr.mean():.3f}")
        print(f"    Max R={max_R:.3f} at v={best_v:.0f} km/s")
        print(f"    Resonances: {total_resonances} total "
              f"(even ℓ: {res_even_total}, odd ℓ: {res_odd_total})")
        print()

        bp_summary.append({
            "bp": label,
            "lambda": lam_bp,
            "R_min": R_arr.min(),
            "R_max": R_arr.max(),
            "R_mean": R_arr.mean(),
            "best_v": best_v,
            "total_resonances": total_resonances,
            "res_even": res_even_total,
            "res_odd": res_odd_total,
        })

    # ── Global summary ────────────────────────────────────────────────
    print("=" * 70)
    print("  T12 SUMMARY: RESONANCE MAP")
    print("=" * 70)
    print()
    print(f"  {'BP':>12} {'λ':>8} {'R_min':>8} {'R_max':>8} {'R_mean':>8} "
          f"{'best_v':>8} {'n_res':>6} {'even':>5} {'odd':>5}")
    print("  " + "─" * 75)
    for s in bp_summary:
        print(f"  {s['bp']:>12} {s['lambda']:>8.2f} {s['R_min']:>8.3f} "
              f"{s['R_max']:>8.3f} {s['R_mean']:>8.3f} "
              f"{s['best_v']:>8.0f} {s['total_resonances']:>6} "
              f"{s['res_even']:>5} {s['res_odd']:>5}")

    print()
    # Key insight: odd-ℓ resonances favor Majorana (weight 3)
    all_res_even = sum(s["res_even"] for s in bp_summary)
    all_res_odd = sum(s["res_odd"] for s in bp_summary)
    r_max_all = max(s["R_max"] for s in bp_summary)
    r_min_all = min(s["R_min"] for s in bp_summary)

    print(f"  Total resonances across all BPs: even={all_res_even}, odd={all_res_odd}")
    print(f"  Global R range: [{r_min_all:.3f}, {r_max_all:.3f}]")
    print()

    if all_res_odd > all_res_even:
        print("  ► ODD-ℓ resonances dominate — Majorana weight w=3 enhances σ_T.")
        print("    This explains R > 2: resonant odd-ℓ contributions are tripled.")
    elif all_res_even > all_res_odd:
        print("  ► EVEN-ℓ resonances dominate — weight w=1, less Majorana enhancement.")
    else:
        print("  ► Even/odd resonances comparable — R governed by non-resonant ℓ.")

    if r_max_all > 2.5:
        print(f"  ► R_max = {r_max_all:.3f} > 2.5: resonance-enhanced Majorana signature!")
    elif r_max_all > 2.0:
        print(f"  ► R_max = {r_max_all:.3f}: moderate enhancement over Born limit R=2.")
    print("=" * 70)
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha", "lambda", "kappa",
        "v_km_s", "sigma_maj", "sigma_dir", "R",
        "n_resonances", "res_even", "res_odd",
        "dominant_l", "n_partialwaves",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T12_resonance_map.py",
        stage="T12 — Resonance Map",
        params={"benchmarks": bp_labels, "n_velocities": len(velocities)},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T12_resonance_map",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_rows)

        rl.add_output(str(out_csv))
        rl.set_notes(f"R range: [{r_min_all:.3f}, {r_max_all:.3f}]. "
                     f"Resonances: even={all_res_even}, odd={all_res_odd}.")
        print(f"Output: {out_csv}")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T12_resonance_map"):
        main()
