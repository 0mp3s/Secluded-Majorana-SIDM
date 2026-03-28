#!/usr/bin/env python3
r"""
test_T11_thermal_averaging.py
==============================
T11: Does the Majorana vs Dirac difference survive Maxwell-Boltzmann averaging?

PHYSICS:
  Real SIDM halos have a velocity DISTRIBUTION, not a single velocity.
  The thermally averaged cross section is:

    ⟨σ/m⟩ = ∫₀^∞ f_MB(v; v₀) · [σ_T(v)/m] dv

  where f_MB(v; v₀) = (4/√π)(v/v₀)² exp(-(v/v₀)²) × (1/v₀)
  is the Maxwell-Boltzmann speed distribution with most-probable speed v₀.

  KEY QUESTION: T9 found R(v) = σ^Maj/σ^Dir varies from ~1 to ~2.18
  at specific velocities. But once averaged over a thermal distribution,
  does the signal wash out?

  If ⟨R⟩ = ⟨σ^Maj⟩/⟨σ^Dir⟩ remains significantly > 1 at all v₀,
  the Majorana signature SURVIVES thermal averaging.

TEST:
  For each benchmark, compute ⟨σ/m⟩ for both Majorana and Dirac
  at characteristic v₀ for:
    - Dwarfs:   v₀ ~ 20-50 km/s
    - LSBs:     v₀ ~ 50-100 km/s
    - Galaxies: v₀ ~ 100-300 km/s
    - Groups:   v₀ ~ 300-700 km/s
    - Clusters: v₀ ~ 700-2000 km/s

Output: data/archive/T11_thermal_averaging_*.csv + output/T11_thermal_averaging.txt
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

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

from _vpm_fast import sigma_T_weighted


# ═══════════════════════════════════════════════════════════════════════
# MAXWELL-BOLTZMANN THERMAL AVERAGING
# ═══════════════════════════════════════════════════════════════════════

def f_mb(v, v0):
    """Maxwell-Boltzmann speed distribution f(v; v₀).
    Normalized so ∫₀^∞ f(v; v₀) dv = 1.
    v₀ = most probable speed = √(2kT/m).
    """
    u = v / v0
    return (4.0 / math.sqrt(math.pi)) * u * u * math.exp(-u * u) / v0


def thermal_average(m_chi, m_phi, alpha, v0_km_s, w_even, w_odd, n_gauss=30):
    """Compute ⟨σ_T/m⟩_MB at characteristic speed v₀.
    Uses Gauss-Legendre quadrature on [0, 4*v₀] (captures >99.9% of integral).
    """
    v_lo = max(5.0, v0_km_s * 0.01)   # avoid v → 0 divergence
    v_hi = v0_km_s * 4.0               # tail cutoff

    def integrand(v):
        s = sigma_T_weighted(m_chi, m_phi, alpha, v, w_even, w_odd)
        return f_mb(v, v0_km_s) * s

    result, _ = quad(integrand, v_lo, v_hi, limit=80, epsrel=1e-3)
    return result


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    bp_labels = _LOCAL_CFG["test_benchmarks"]
    bps = GC.benchmarks(*bp_labels)

    # Characteristic speeds for different halo types
    halo_types = [
        ("dwarf",   30.0),
        ("LSB",     70.0),
        ("MW-size", 200.0),
        ("group",   500.0),
        ("cluster", 1200.0),
    ]

    print()
    print("╔═══════════════════════════════════════════════════════════════════╗")
    print("║  T11: THERMAL AVERAGING — Does Majorana signal survive ⟨·⟩_MB? ║")
    print("╚═══════════════════════════════════════════════════════════════════╝")
    print()
    print(f"  Halo types: {[h[0] for h in halo_types]}")
    print(f"  v₀ values:  {[h[1] for h in halo_types]} km/s")
    print(f"  Benchmarks: {bp_labels}")
    print()

    all_rows = []
    all_R_thermal = []

    for bp in bps:
        label = bp["label"]
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3
        alpha = bp["alpha"]

        print(f"  ── {label} (m_χ={m_chi:.3f} GeV, m_φ={bp['m_phi_MeV']:.3f} MeV, α={alpha:.4e}) ──")

        for halo_name, v0 in halo_types:
            # Single-velocity R(v₀) for comparison
            sigma_maj_v0 = sigma_T_weighted(m_chi, m_phi, alpha, v0, 1.0, 3.0)
            sigma_dir_v0 = sigma_T_weighted(m_chi, m_phi, alpha, v0, 1.0, 1.0)
            R_single = sigma_maj_v0 / sigma_dir_v0 if sigma_dir_v0 > 0 else float('nan')

            # Thermally averaged
            avg_maj = thermal_average(m_chi, m_phi, alpha, v0, 1.0, 3.0)
            avg_dir = thermal_average(m_chi, m_phi, alpha, v0, 1.0, 1.0)
            R_thermal = avg_maj / avg_dir if avg_dir > 0 else float('nan')

            # How much does averaging change R?
            wash_pct = (R_thermal - R_single) / R_single * 100 if R_single > 0 else float('nan')

            all_R_thermal.append(R_thermal)

            survived = R_thermal > 1.05  # > 5% effect

            print(f"    {halo_name:>8} (v₀={v0:>6.0f} km/s): "
                  f"R_single={R_single:.3f}, ⟨R⟩_MB={R_thermal:.3f}, "
                  f"wash={wash_pct:+.1f}%  {'✓' if survived else '~'}")

            all_rows.append({
                "bp": label,
                "m_chi_GeV": m_chi,
                "m_phi_MeV": bp["m_phi_MeV"],
                "alpha": alpha,
                "halo_type": halo_name,
                "v0_km_s": v0,
                "sigma_maj_single": sigma_maj_v0,
                "sigma_dir_single": sigma_dir_v0,
                "R_single": R_single,
                "sigma_maj_thermal": avg_maj,
                "sigma_dir_thermal": avg_dir,
                "R_thermal": R_thermal,
                "washout_pct": wash_pct,
            })

        print()

    # ── Global summary ────────────────────────────────────────────────
    R_arr = np.array([r for r in all_R_thermal if math.isfinite(r)])

    print("=" * 70)
    print("  T11 SUMMARY: THERMAL AVERAGING")
    print("=" * 70)
    print()
    print(f"  ⟨R⟩_MB range: [{R_arr.min():.3f}, {R_arr.max():.3f}]")
    print(f"  ⟨R⟩_MB mean:  {R_arr.mean():.3f}")
    print(f"  ⟨R⟩_MB > 1.05 (5% effect): {np.sum(R_arr > 1.05)}/{len(R_arr)} measurements")
    print(f"  ⟨R⟩_MB > 1.10 (10% effect): {np.sum(R_arr > 1.10)}/{len(R_arr)} measurements")
    print(f"  ⟨R⟩_MB > 1.50 (50% effect): {np.sum(R_arr > 1.50)}/{len(R_arr)} measurements")
    print()

    if R_arr.min() > 1.05:
        print("  ► CONCLUSION: Majorana signal SURVIVES thermal averaging at ALL v₀!")
        print("    ⟨R⟩_MB never drops below 1.05 — the effect is robust.")
    elif R_arr.mean() > 1.10:
        print("  ► CONCLUSION: Majorana signal PARTIALLY survives.")
        print("    Mean ⟨R⟩_MB > 1.10 but washed out at some velocities.")
    else:
        print("  ► CONCLUSION: Thermal averaging WASHES OUT the signal.")
        print("    ⟨R⟩_MB ≈ 1 — cannot distinguish Majorana from Dirac via ⟨σ/m⟩.")

    # ── Washout by halo type ──────────────────────────────────────────
    print()
    print("  Breakdown by halo type:")
    for halo_name, v0 in halo_types:
        rows_halo = [r for r in all_rows if r["halo_type"] == halo_name]
        mean_R = np.mean([r["R_thermal"] for r in rows_halo])
        mean_wash = np.mean([r["washout_pct"] for r in rows_halo
                            if math.isfinite(r["washout_pct"])])
        print(f"    {halo_name:>8}: ⟨R⟩_MB={mean_R:.3f}, mean washout={mean_wash:+.1f}%")

    print("=" * 70)
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha", "halo_type", "v0_km_s",
        "sigma_maj_single", "sigma_dir_single", "R_single",
        "sigma_maj_thermal", "sigma_dir_thermal", "R_thermal", "washout_pct",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T11_thermal_averaging.py",
        stage="T11 — Thermal Averaging",
        params={"benchmarks": bp_labels, "halo_types": [h[0] for h in halo_types]},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T11_thermal_averaging",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_rows)

        rl.add_output(str(out_csv))
        rl.set_notes(f"⟨R⟩_MB range=[{R_arr.min():.3f}, {R_arr.max():.3f}], "
                     f"mean={R_arr.mean():.3f}. "
                     f"{len(bp_labels)} BPs × {len(halo_types)} halo types.")
        print(f"Output: {out_csv}")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T11_thermal_averaging"):
        main()
