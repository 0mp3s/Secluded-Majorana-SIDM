#!/usr/bin/env python3
r"""
test_T10_dirac_mimicry.py
==========================
T10: Can a Dirac fermion with DIFFERENT parameters mimic a Majorana σ_T(v)?

PHYSICS:
  T9 showed R(v) = σ_T^Maj/σ_T^Dir ∈ [1.0, 2.18] at SAME parameters.
  But an experimentalist doesn't know the underlying parameters — they
  observe σ_T/m at several velocities and fit to a model.

  QUESTION: Given a Majorana benchmark (m_χ, m_φ, α), can we find
  Dirac parameters (m'_χ, m'_φ, α') that reproduce the SAME σ_T(v)
  curve across the observational velocity range?

  If NO good fit exists → Majorana is observationally distinguishable
  from Dirac even with full parameter freedom.

  If YES → the degeneracy is exact and we need something else to break it.

METHOD:
  1. Compute Majorana σ_T(v) at 10 velocities spanning dwarfs → clusters
  2. Use scipy.optimize.minimize to find best-fit Dirac (m'_χ, m'_φ, α')
     minimizing χ² = Σ_i [(σ_Dir(v_i) - σ_Maj(v_i))² / σ_Maj(v_i)²]
  3. Constraint: m'_χ > 0, m'_φ > 0, α' > 0 + perturbativity α' < 1
  4. Report best-fit χ², relative deviations, parameter shifts

Output: data/archive/T10_dirac_mimicry_*.csv + output/T10_dirac_mimicry.txt
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

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
# DIRAC MIMICRY FIT
# ═══════════════════════════════════════════════════════════════════════

def compute_majorana_curve(m_chi, m_phi, alpha, velocities):
    """Compute Majorana σ_T/m at each velocity."""
    return np.array([sigma_T_weighted(m_chi, m_phi, alpha, v, 1.0, 3.0)
                     for v in velocities])


def chi2_dirac(params, velocities, sigma_maj_target):
    """χ² between Dirac(m'_χ, m'_φ, α') and Majorana target.
    params = [log10(m_chi'), log10(m_phi'), log10(alpha')]
    """
    m_chi_d = 10.0 ** params[0]
    m_phi_d = 10.0 ** params[1]
    alpha_d = 10.0 ** params[2]

    chi2 = 0.0
    for i, v in enumerate(velocities):
        s_dir = sigma_T_weighted(m_chi_d, m_phi_d, alpha_d, v, 1.0, 1.0)
        s_maj = sigma_maj_target[i]
        if s_maj > 0:
            chi2 += ((s_dir - s_maj) / s_maj) ** 2
        elif s_dir > 0:
            chi2 += 1e6

    return chi2


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    bp_labels = _LOCAL_CFG["test_benchmarks"]
    bps = GC.benchmarks(*bp_labels)

    # Velocity grid: dwarfs → clusters (5 points, log-spaced for speed)
    velocities = np.logspace(np.log10(30), np.log10(1500), 5)  # km/s

    print()
    print("╔═══════════════════════════════════════════════════════════════════╗")
    print("║  T10: DIRAC MIMICRY — Can Dirac mimic Majorana σ_T(v)?         ║")
    print("╚═══════════════════════════════════════════════════════════════════╝")
    print()
    print(f"  Velocities: {len(velocities)} points from {velocities[0]:.0f} to {velocities[-1]:.0f} km/s")
    print(f"  Benchmarks: {bp_labels}")
    print()

    all_rows = []
    results_summary = []

    for bp in bps:
        label = bp["label"]
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] * 1e-3  # → GeV
        alpha = bp["alpha"]

        print(f"  ── {label} (m_χ={m_chi:.3f} GeV, m_φ={bp['m_phi_MeV']:.3f} MeV, α={alpha:.4e}) ──")

        # Step 1: compute Majorana target
        sigma_maj = compute_majorana_curve(m_chi, m_phi, alpha, velocities)
        print(f"    Majorana σ/m range: [{sigma_maj.min():.4e}, {sigma_maj.max():.4e}] cm²/g")

        # Step 2: fit Dirac parameters via Nelder-Mead minimization
        # Initial guess: same params but double alpha (since R ≈ 2)
        x0 = np.array([np.log10(m_chi), np.log10(m_phi),
                        np.log10(min(alpha * 2.0, 0.5))])

        result = minimize(
            chi2_dirac,
            x0,
            args=(velocities, sigma_maj),
            method='Nelder-Mead',
            options={'maxiter': 300, 'xatol': 0.01, 'fatol': 1e-8},
        )

        m_chi_fit = 10.0 ** result.x[0]
        m_phi_fit = 10.0 ** result.x[1]
        alpha_fit = 10.0 ** result.x[2]
        chi2_best = result.fun
        n_velocities = len(velocities)
        chi2_per_dof = chi2_best / max(1, n_velocities - 3)

        print(f"    Best-fit Dirac: m'_χ={m_chi_fit:.3f} GeV, "
              f"m'_φ={m_phi_fit*1e3:.3f} MeV, α'={alpha_fit:.4e}")
        print(f"    χ²={chi2_best:.6f}, χ²/dof={chi2_per_dof:.6f}")
        print(f"    Parameter shifts: Δm_χ={m_chi_fit/m_chi:.3f}×, "
              f"Δm_φ={m_phi_fit/m_phi:.3f}×, Δα={alpha_fit/alpha:.3f}×")

        # Step 3: compute residuals at each velocity
        sigma_dir_fit = np.array([
            sigma_T_weighted(m_chi_fit, m_phi_fit, alpha_fit, v, 1.0, 1.0)
            for v in velocities
        ])

        max_dev = 0.0
        for i, v in enumerate(velocities):
            if sigma_maj[i] > 0:
                dev = (sigma_dir_fit[i] - sigma_maj[i]) / sigma_maj[i] * 100
            else:
                dev = 0.0
            max_dev = max(max_dev, abs(dev))
            all_rows.append({
                "bp": label,
                "m_chi_GeV": m_chi,
                "m_phi_MeV": bp["m_phi_MeV"],
                "alpha": alpha,
                "v_km_s": v,
                "sigma_maj": sigma_maj[i],
                "sigma_dir_fit": sigma_dir_fit[i],
                "deviation_pct": dev,
                "m_chi_fit": m_chi_fit,
                "m_phi_fit_MeV": m_phi_fit * 1e3,
                "alpha_fit": alpha_fit,
                "chi2": chi2_best,
                "chi2_per_dof": chi2_per_dof,
            })

        distinguishable = chi2_per_dof > 0.01  # > 1% residuals per dof
        status = "DISTINGUISHABLE" if distinguishable else "DEGENERATE"
        print(f"    Max deviation: {max_dev:.2f}%")
        print(f"    → {status}")
        print()

        results_summary.append({
            "bp": label,
            "chi2_per_dof": chi2_per_dof,
            "max_dev_pct": max_dev,
            "status": status,
            "m_chi_fit": m_chi_fit,
            "m_phi_fit_MeV": m_phi_fit * 1e3,
            "alpha_fit": alpha_fit,
        })

    # ── Global summary ────────────────────────────────────────────────
    print("=" * 70)
    print("  T10 SUMMARY: DIRAC MIMICRY TEST")
    print("=" * 70)
    print()
    n_dist = sum(1 for r in results_summary if r["status"] == "DISTINGUISHABLE")
    n_degen = sum(1 for r in results_summary if r["status"] == "DEGENERATE")
    print(f"  DISTINGUISHABLE: {n_dist}/{len(results_summary)} BPs")
    print(f"  DEGENERATE:      {n_degen}/{len(results_summary)} BPs")
    print()

    for r in results_summary:
        print(f"  {r['bp']:>12}: χ²/dof={r['chi2_per_dof']:.6f}, "
              f"max_dev={r['max_dev_pct']:.2f}%, {r['status']}")
        print(f"               fit: m'_χ={r['m_chi_fit']:.3f} GeV, "
              f"m'_φ={r['m_phi_fit_MeV']:.3f} MeV, α'={r['alpha_fit']:.4e}")

    print()
    if n_dist > n_degen:
        print("  ► CONCLUSION: Majorana σ_T(v) is NOT fully mimicked by Dirac")
        print("    even with free parameters. The velocity dependence differs.")
    elif n_degen > n_dist:
        print("  ► CONCLUSION: Dirac CAN approximate Majorana σ_T(v) via parameter shifts.")
        print("    Need additional observables (relic, transfer CS) to break degeneracy.")
    else:
        print("  ► CONCLUSION: Mixed results — some BPs distinguishable, some not.")
        print("    Regime-dependent: resonances help, Born limit harder.")
    print("=" * 70)
    print()

    # ── CSV output ────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha", "v_km_s",
        "sigma_maj", "sigma_dir_fit", "deviation_pct",
        "m_chi_fit", "m_phi_fit_MeV", "alpha_fit",
        "chi2", "chi2_per_dof",
    ]

    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T10_dirac_mimicry.py",
        stage="T10 — Dirac Mimicry",
        params={"benchmarks": bp_labels, "n_velocities": len(velocities)},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T10_dirac_mimicry",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_rows)

        rl.add_output(str(out_csv))
        rl.set_notes(f"Distinguishable: {n_dist}/{len(results_summary)}. "
                     f"Degenerate: {n_degen}/{len(results_summary)}.")
        print(f"Output: {out_csv}")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T10_dirac_mimicry"):
        main()
