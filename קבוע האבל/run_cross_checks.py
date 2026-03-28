#!/usr/bin/env python3
"""
Layer 8 cross-checks: θ_i, f/M_Pl, m_σ/H₀ sensitivity scans.
Writes each result to CSV immediately — safe to interrupt and resume.
"""
import csv
import sys
import os
from pathlib import Path

# ── path setup ────────────────────────────────────────────────────────
_HERE = Path(__file__).resolve().parent
_SIDM = _HERE.parent
sys.path.insert(0, str(_SIDM))
sys.path.insert(0, str(_SIDM / "core"))
sys.path.insert(0, str(_SIDM / "The_Lagernizant_integral_SIDM"))

from lagrangian_path_integral import (
    solve_friedmann_sigma, H0_KM_S_MPC, RHO_LAMBDA,
    LAMBDA_D_GEV, F_OVER_M_PL, M_PL_GEV, THETA_I_DEFAULT,
)
import numpy as np

CSV_PATH = _HERE / "cross_checks.csv"
FIELDS = [
    "scan", "param_name", "param_value",
    "H0_km_s_Mpc", "delta_H0_pct", "theta_final_rad", "theta_final_deg",
    "V_sigma_GeV4", "V_over_rhoL", "m_sigma_GeV", "m_sigma_over_H0",
    "converged", "iterations",
]
LD_BEST = 2.28e-12

# ── resume logic ──────────────────────────────────────────────────────
def load_done():
    """Return set of (scan, param_value) already in CSV."""
    done = set()
    if CSV_PATH.exists():
        with open(CSV_PATH, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                done.add((row["scan"], row["param_value"]))
    return done


def append_row(row):
    """Append one row to CSV, create header if needed."""
    write_header = not CSV_PATH.exists() or CSV_PATH.stat().st_size == 0
    with open(CSV_PATH, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def run_one(scan_name, param_name, param_value, **kwargs):
    """Run one solve_friedmann_sigma call and return a CSV row dict."""
    f_gev = kwargs.get("f_gev", F_OVER_M_PL * M_PL_GEV)
    ld = kwargs.get("Lambda_d", LD_BEST)
    m_sigma = ld**2 / f_gev

    res = solve_friedmann_sigma(
        m_chi_gev=kwargs.get("m_chi_gev", 30.0),
        alpha=kwargs.get("alpha", 0.01),
        theta_i=kwargs.get("theta_i", THETA_I_DEFAULT),
        f_gev=f_gev,
        Lambda_d=ld,
    )

    conv = res.get("converged", False)
    h0 = res.get("H0_km_s_Mpc", 0)
    delta = (h0 - H0_KM_S_MPC) / H0_KM_S_MPC * 100 if h0 > 0 else -999
    tf = res.get("theta_final", 0)
    Vf = res.get("V_sigma_final_GeV4", 0)
    H0g = res.get("H0_GeV", 1e-99)

    return {
        "scan": scan_name,
        "param_name": param_name,
        "param_value": str(param_value),
        "H0_km_s_Mpc": f"{h0:.4f}",
        "delta_H0_pct": f"{delta:.4f}",
        "theta_final_rad": f"{tf:.8f}",
        "theta_final_deg": f"{np.degrees(tf):.4f}",
        "V_sigma_GeV4": f"{Vf:.6e}",
        "V_over_rhoL": f"{Vf / RHO_LAMBDA:.6f}",
        "m_sigma_GeV": f"{m_sigma:.6e}",
        "m_sigma_over_H0": f"{m_sigma / H0g:.4f}" if conv else "N/A",
        "converged": str(conv),
        "iterations": str(res.get("iterations", "N/A")),
    }


def main():
    done = load_done()
    total = 0
    skipped = 0

    # ── 6.1: θ_i sensitivity ─────────────────────────────────────────
    print("=" * 70)
    print("  SCAN 6.1: theta_i sensitivity  (Lambda_d = %.2e)" % LD_BEST)
    print("=" * 70)

    for ti in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.1]:
        key = ("theta_i", str(ti))
        if key in done:
            print("  [SKIP] theta_i = %.1f (already done)" % ti)
            skipped += 1
            continue
        print("  [RUN]  theta_i = %.1f ..." % ti, end="", flush=True)
        row = run_one("theta_i", "theta_i", ti, theta_i=ti)
        append_row(row)
        total += 1
        print("  H0 = %s km/s/Mpc  (%s)" % (row["H0_km_s_Mpc"], row["delta_H0_pct"] + "%"))

    # ── 6.2: f/M_Pl sensitivity ──────────────────────────────────────
    print()
    print("=" * 70)
    print("  SCAN 6.2: f/M_Pl sensitivity")
    print("=" * 70)

    for f_ratio in [0.05, 0.10, 0.15, 0.20, 0.24, 0.30, 0.40, 0.50, 0.75, 1.00]:
        key = ("f_over_mpl", str(f_ratio))
        if key in done:
            print("  [SKIP] f/M_Pl = %.2f (already done)" % f_ratio)
            skipped += 1
            continue
        f_gev = f_ratio * M_PL_GEV
        print("  [RUN]  f/M_Pl = %.2f ..." % f_ratio, end="", flush=True)
        row = run_one("f_over_mpl", "f/M_Pl", f_ratio, f_gev=f_gev)
        append_row(row)
        total += 1
        print("  H0 = %s km/s/Mpc" % row["H0_km_s_Mpc"])

    # ── 6.3: Lambda_d fine scan ───────────────────────────────────────
    print()
    print("=" * 70)
    print("  SCAN 6.3: Lambda_d fine scan (f/M_Pl = 0.24, theta_i = 2.0)")
    print("=" * 70)

    for ld in np.linspace(1.5e-12, 3.0e-12, 31):
        ld_str = "%.6e" % ld
        key = ("lambda_d", ld_str)
        if key in done:
            skipped += 1
            continue
        print("  [RUN]  Ld = %s ..." % ld_str, end="", flush=True)
        row = run_one("lambda_d", "Lambda_d_GeV", float(ld_str), Lambda_d=ld)
        append_row(row)
        total += 1
        print("  H0 = %s km/s/Mpc" % row["H0_km_s_Mpc"])

    # ── 6.4: m_σ / H₀ structural check ──────────────────────────────
    print()
    print("=" * 70)
    print("  SCAN 6.4: m_sigma / H0 ratio vs f/M_Pl")
    print("=" * 70)

    for f_ratio in [0.05, 0.10, 0.20, 0.24, 0.30, 0.50, 0.75, 1.00, 2.00]:
        key = ("msig_H0_ratio", str(f_ratio))
        if key in done:
            skipped += 1
            continue
        f_gev = f_ratio * M_PL_GEV
        print("  [RUN]  f/M_Pl = %.2f ..." % f_ratio, end="", flush=True)
        row = run_one("msig_H0_ratio", "f/M_Pl", f_ratio, f_gev=f_gev)
        append_row(row)
        total += 1
        print("  m_sig/H0 = %s" % row["m_sigma_over_H0"])

    # ── 6.5: Hubble tension — SH0ES vs Planck ────────────────────────
    print()
    print("=" * 70)
    print("  SCAN 6.5: Lambda_d for SH0ES H0=73.0 vs Planck H0=67.4")
    print("=" * 70)

    for ld in np.linspace(2.20e-12, 2.50e-12, 31):
        ld_str = "%.6e" % ld
        key = ("hubble_tension", ld_str)
        if key in done:
            skipped += 1
            continue
        print("  [RUN]  Ld = %s ..." % ld_str, end="", flush=True)
        row = run_one("hubble_tension", "Lambda_d_GeV", float(ld_str), Lambda_d=ld)
        append_row(row)
        total += 1
        h0v = float(row["H0_km_s_Mpc"])
        tag = ""
        if abs(h0v - 67.4) < 0.5:
            tag = " <-- PLANCK"
        elif abs(h0v - 73.0) < 0.5:
            tag = " <-- SH0ES"
        print("  H0 = %s%s" % (row["H0_km_s_Mpc"], tag))

    print()
    print("=" * 70)
    print("  DONE: %d new runs, %d skipped (resumed)" % (total, skipped))
    print("  Results: %s" % CSV_PATH)
    print("=" * 70)


if __name__ == "__main__":
    main()
