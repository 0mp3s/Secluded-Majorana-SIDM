#!/usr/bin/env python3
r"""
test_T7_boltzmann_all_bps.py
==============================
T7: Extend T6 relic density calculation to ALL benchmarks + p-wave correction.

PHYSICS:
  T6 used 4 BPs to identify the α convention via relic density Ωh².
  Now we extend to ALL available BPs and add:

  1. s-wave Boltzmann for all BPs (using verified α = y²/4π convention)
  2. p-wave correction: ⟨σv⟩ = a + b/x_fo (b = p-wave coefficient)
     For Majorana fermions: a = π α²/(4m²), b ≈ 3a/x_fo (velocity suppressed)
  3. Compute freeze-out temperature x_fo = m/T_fo for each BP
  4. Compare numerical Boltzmann vs Kolb-Turner analytic formula

  KEY PREDICTIONS:
    - Only BP1 (MAP_relic) is tuned to give Ωh² ≈ 0.12
    - Other BPs were chosen for SIDM σ_T fits, not relic → expect Ωh² ≠ 0.12
    - p-wave correction is small (< few %) for s-wave-dominated annihilation

Output: data/archive/T7_boltzmann_all_bps_*.csv + output/T7_boltzmann_all_bps.txt
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

# ── import Boltzmann solver from v27 ──────────────────────────────────
from v27_boltzmann_relic import (
    solve_boltzmann,
    sigma_v_swave,
    Y_to_omega_h2,
    kolb_turner_swave,
)

# ── local config ──────────────────────────────────────────────────────
_LOCAL_CFG_PATH = _HERE / "data" / "config.json"
with open(_LOCAL_CFG_PATH, "r", encoding="utf-8") as _fh:
    _LOCAL_CFG = json.load(_fh)

_CC = GC.cosmological_constants()
OMEGA_CDM_TARGET = _CC.get("omega_h2_target", 0.120)


# ═══════════════════════════════════════════════════════════════════════
# P-WAVE BOLTZMANN SOLVER (extends v27's s-wave)
# ═══════════════════════════════════════════════════════════════════════

def sigma_v_with_pwave(alpha_d, m_chi, x):
    """⟨σv⟩ including leading p-wave correction.
    
    s-wave:  a = π α² / (4 m²)
    p-wave:  b ≈ 3a (for Majorana, scalar mediator)
    Total:   ⟨σv⟩ ≈ a + b·(6/x) where x = m/T
    
    The factor 6/x comes from ⟨v²⟩ = 6T/m for MB distribution.
    """
    a = math.pi * alpha_d**2 / (4.0 * m_chi**2)
    b = 3.0 * a  # leading p-wave for Majorana scalar mediator
    return a + b * 6.0 / x


def solve_boltzmann_pwave(m_chi, alpha_d, x_start=1.0, x_end=1000.0,
                          n_steps=10000, g_chi=2):
    """Boltzmann ODE with velocity-dependent ⟨σv⟩(x).
    Same structure as v27 solve_boltzmann but with p-wave correction.
    """
    from v27_boltzmann_relic import g_star_S, Y_eq_full, g_star_rho

    _PC = GC.physical_constants()
    M_PL = _PC["m_pl_GeV"]

    x_arr = np.linspace(x_start, x_end, n_steps)
    h = x_arr[1] - x_arr[0]

    Y = Y_eq_full(x_start, m_chi, g_chi)
    Y_arr = np.zeros(n_steps)
    Y_arr[0] = Y

    for i in range(n_steps - 1):
        x = x_arr[i]
        T = m_chi / x
        sv = sigma_v_with_pwave(alpha_d, m_chi, x)

        g_s = g_star_S(T)
        g_rho = g_star_rho(T)

        # Hubble parameter
        H = math.sqrt(8.0 * math.pi**3 * g_rho / 90.0) * T**2 / M_PL
        s = (2.0 * math.pi**2 / 45.0) * g_s * T**3

        Y_eq = Y_eq_full(x, m_chi, g_chi)

        # dY/dx = -(s ⟨σv⟩ / (x H)) (Y² - Y_eq²)
        dYdx = -(s * sv / (x * H)) * (Y**2 - Y_eq**2)

        # RK4
        k1 = h * dYdx

        x2 = x + 0.5 * h
        T2 = m_chi / x2
        sv2 = sigma_v_with_pwave(alpha_d, m_chi, x2)
        H2 = math.sqrt(8.0 * math.pi**3 * g_star_rho(T2) / 90.0) * T2**2 / M_PL
        s2 = (2.0 * math.pi**2 / 45.0) * g_star_S(T2) * T2**3
        Y_eq2 = Y_eq_full(x2, m_chi, g_chi)
        Y2 = Y + 0.5 * k1
        k2 = h * (-(s2 * sv2 / (x2 * H2)) * (Y2**2 - Y_eq2**2))

        Y3 = Y + 0.5 * k2
        k3 = h * (-(s2 * sv2 / (x2 * H2)) * (Y3**2 - Y_eq2**2))

        x4 = x + h
        T4 = m_chi / x4
        sv4 = sigma_v_with_pwave(alpha_d, m_chi, x4)
        H4 = math.sqrt(8.0 * math.pi**3 * g_star_rho(T4) / 90.0) * T4**2 / M_PL
        s4 = (2.0 * math.pi**2 / 45.0) * g_star_S(T4) * T4**3
        Y_eq4 = Y_eq_full(x4, m_chi, g_chi)
        Y4 = Y + k3
        k4 = h * (-(s4 * sv4 / (x4 * H4)) * (Y4**2 - Y_eq4**2))

        Y = Y + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        if Y < 0:
            Y = 0.0
        Y_arr[i + 1] = Y

    return x_arr, Y_arr


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    # Use ALL available BPs, not just the 4 from local config
    all_bps = GC.all_benchmarks()
    bp_labels = [bp["label"] for bp in all_bps]

    print()
    print("╔═══════════════════════════════════════════════════════════════════╗")
    print("║  T7: BOLTZMANN ALL BPs — Relic density + p-wave for all BPs    ║")
    print("╚═══════════════════════════════════════════════════════════════════╝")
    print()
    print(f"  Target: Ωh² = {OMEGA_CDM_TARGET}")
    print(f"  Benchmarks: {bp_labels}")
    print(f"  Convention: α = y²/(4π) — verified by T6")
    print()

    all_rows = []

    for bp in all_bps:
        label = bp["label"]
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3
        alpha = bp["alpha"]

        print(f"  ── {label} (m_χ={m_chi:.3f} GeV, m_φ={bp['m_phi_MeV']:.3f} MeV, α={alpha:.4e}) ──")

        # ── s-wave ⟨σv⟩ ──
        sv0 = sigma_v_swave(alpha, m_chi)
        print(f"    ⟨σv⟩_s = {sv0:.4e} GeV⁻²")

        # ── Numerical Boltzmann (s-wave only) ──
        x_arr, Y_arr = solve_boltzmann(m_chi, sv0)
        Y_inf_swave = Y_arr[-1]
        omega_swave = Y_to_omega_h2(Y_inf_swave, m_chi)

        # ── Kolb-Turner analytic (s-wave) ──
        x_fo_kt, Y_inf_kt = kolb_turner_swave(m_chi, sv0)
        omega_kt = Y_to_omega_h2(Y_inf_kt, m_chi)

        # ── Numerical Boltzmann (s+p-wave) ──
        x_arr_pw, Y_arr_pw = solve_boltzmann_pwave(m_chi, alpha)
        Y_inf_pwave = Y_arr_pw[-1]
        omega_pwave = Y_to_omega_h2(Y_inf_pwave, m_chi)

        # ── p-wave correction ──
        pwave_correction_pct = (omega_pwave - omega_swave) / omega_swave * 100 \
            if omega_swave > 0 else 0.0

        # ── Deviation from target ──
        delta_swave = (omega_swave - OMEGA_CDM_TARGET) / OMEGA_CDM_TARGET * 100
        delta_pwave = (omega_pwave - OMEGA_CDM_TARGET) / OMEGA_CDM_TARGET * 100
        delta_kt = (omega_kt - OMEGA_CDM_TARGET) / OMEGA_CDM_TARGET * 100

        # ── Freeze-out temperature ──
        # Find x where Y departs from Y_eq (Y > 3·Y_eq)
        from v27_boltzmann_relic import Y_eq_full
        x_fo_num = 20.0  # default
        for i, x in enumerate(x_arr):
            y_eq = Y_eq_full(x, m_chi)
            if y_eq > 0 and Y_arr[i] > 3.0 * y_eq and x > 5.0:
                x_fo_num = x
                break

        T_fo = m_chi / x_fo_num

        match_target = abs(delta_swave) < 20  # within 20% of Ωh² = 0.12

        print(f"    x_fo = {x_fo_num:.1f} (T_fo = {T_fo*1e3:.1f} MeV)")
        print(f"    Ωh²_swave = {omega_swave:.4f} (Δ = {delta_swave:+.1f}%)")
        print(f"    Ωh²_pwave = {omega_pwave:.4f} (Δ = {delta_pwave:+.1f}%)")
        print(f"    Ωh²_KT    = {omega_kt:.4f} (Δ = {delta_kt:+.1f}%)")
        print(f"    p-wave correction: {pwave_correction_pct:+.2f}%")
        print(f"    {'✓ matches target' if match_target else '✗ does NOT match target'}")
        print()

        all_rows.append({
            "bp": label,
            "m_chi_GeV": m_chi,
            "m_phi_MeV": bp["m_phi_MeV"],
            "alpha": alpha,
            "sv0_GeV2": sv0,
            "x_fo_numerical": x_fo_num,
            "x_fo_KT": x_fo_kt,
            "T_fo_MeV": T_fo * 1e3,
            "Y_inf_swave": Y_inf_swave,
            "Y_inf_pwave": Y_inf_pwave,
            "Y_inf_KT": Y_inf_kt,
            "omega_h2_swave": omega_swave,
            "omega_h2_pwave": omega_pwave,
            "omega_h2_KT": omega_kt,
            "omega_target": OMEGA_CDM_TARGET,
            "delta_swave_pct": delta_swave,
            "delta_pwave_pct": delta_pwave,
            "delta_KT_pct": delta_kt,
            "pwave_correction_pct": pwave_correction_pct,
        })

    # ── Global summary ────────────────────────────────────────────────
    print("=" * 70)
    print("  T7 SUMMARY: BOLTZMANN + p-WAVE FOR ALL BPs")
    print("=" * 70)
    print()
    print(f"  {'BP':>12} {'Ωh²_s':>10} {'Ωh²_pw':>10} {'Ωh²_KT':>10} "
          f"{'Δ_s%':>8} {'pw_corr%':>9} {'x_fo':>6}")
    print("  " + "─" * 68)

    n_match = 0
    for r in all_rows:
        marker = "✓" if abs(r["delta_swave_pct"]) < 20 else " "
        if abs(r["delta_swave_pct"]) < 20:
            n_match += 1
        print(f"  {r['bp']:>12} {r['omega_h2_swave']:>10.4f} "
              f"{r['omega_h2_pwave']:>10.4f} {r['omega_h2_KT']:>10.4f} "
              f"{r['delta_swave_pct']:>+8.1f} {r['pwave_correction_pct']:>+9.2f} "
              f"{r['x_fo_numerical']:>6.1f} {marker}")

    print()
    print(f"  BPs matching Ωh²_target (±20%): {n_match}/{len(all_rows)}")
    print()

    # ── p-wave assessment ──
    pw_corrections = [abs(r["pwave_correction_pct"]) for r in all_rows]
    max_pw = max(pw_corrections)
    mean_pw = sum(pw_corrections) / len(pw_corrections)

    print(f"  p-wave correction: |max| = {max_pw:.2f}%, mean = {mean_pw:.2f}%")
    if max_pw < 5.0:
        print("  ► p-wave correction is NEGLIGIBLE (< 5%) — s-wave dominates.")
        print("    This confirms Majorana annihilation is s-wave for scalar mediator.")
    else:
        print("  ► p-wave correction is SIGNIFICANT — must include for precision.")

    print()

    # ── KT vs numerical ──
    kt_deltas = [abs(r["omega_h2_swave"] - r["omega_h2_KT"]) / r["omega_h2_swave"] * 100
                 for r in all_rows if r["omega_h2_swave"] > 0]
    print(f"  KT analytic vs numerical: max deviation = {max(kt_deltas):.1f}%, "
          f"mean = {sum(kt_deltas)/len(kt_deltas):.1f}%")
    print("=" * 70)
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha", "sv0_GeV2",
        "x_fo_numerical", "x_fo_KT", "T_fo_MeV",
        "Y_inf_swave", "Y_inf_pwave", "Y_inf_KT",
        "omega_h2_swave", "omega_h2_pwave", "omega_h2_KT",
        "omega_target", "delta_swave_pct", "delta_pwave_pct",
        "delta_KT_pct", "pwave_correction_pct",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T7_boltzmann_all_bps.py",
        stage="T7 — Boltzmann All BPs",
        params={"benchmarks": bp_labels, "omega_target": OMEGA_CDM_TARGET},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T7_boltzmann_all_bps",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_rows)

        rl.add_output(str(out_csv))
        rl.set_notes(f"{n_match}/{len(all_rows)} BPs within 20% of Ωh²={OMEGA_CDM_TARGET}. "
                     f"Max p-wave correction: {max_pw:.2f}%.")
        print(f"Output: {out_csv}")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T7_boltzmann_all_bps"):
        main()
