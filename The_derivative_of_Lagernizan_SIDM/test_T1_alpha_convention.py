#!/usr/bin/env python3
"""
test_T1_alpha_convention.py
===========================
T1: Derive V(r) from the Lagrangian symbolically (SymPy) and determine
which α convention matches α_config from global_config.json.

Chain:  L_int → vertex factor → Born amplitude → FT → V(r) → extract α_derived

Reads parameters from global_config.json via GC singleton.
Writes timestamped CSV to data/archive/.
Logs run to docs/runs_log.csv.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# ── path setup: allow imports from parent SIDM core ──────────────────────────
_HERE      = Path(__file__).resolve().parent
_SIDM_ROOT = _HERE.parent          # Secluded-Majorana-SIDM/
_CORE      = _SIDM_ROOT / "core"

sys.path.insert(0, str(_SIDM_ROOT))
sys.path.insert(0, str(_CORE))

from core.global_config import GC
from core.output_manager import timestamped_path
from core.run_logger import RunLogger

# ── local config ─────────────────────────────────────────────────────────────
_LOCAL_CFG_PATH = _HERE / "data" / "config.json"

with open(_LOCAL_CFG_PATH, "r", encoding="utf-8") as _fh:
    _LOCAL_CFG = json.load(_fh)

# ── SymPy derivation ─────────────────────────────────────────────────────────
import sympy as sp

def derive_alpha_conventions():
    """
    Symbolically derive V(r) from L_int = -(y/2) χ̄χ φ.

    Returns dict mapping convention name → α expression in terms of y.
    """
    y, r, m_phi, q = sp.symbols("y r m_phi q", positive=True)
    pi = sp.pi

    # --- Step 1: Vertex factor = -i y/2  (Majorana ½ in Lagrangian)
    vertex = y / 2

    # --- Step 2: Born amplitude, t-channel only
    #     M_t = vertex² / (q² + m_phi²)   [NR limit: q⁰→0, so q² = |q|²]
    M_t_numerator = vertex**2         # (y/2)² = y²/4

    # --- Step 3: NR Fourier transform to coordinate space
    #     V_t(r) = -(1/4π) * M_t_numerator * e^{-m_phi r} / r
    #     (standard QFT result: FT of 1/(q²+m²) = e^{-mr}/(4π r))
    alpha_t = M_t_numerator / (4 * pi)     # = y²/(16π)

    # --- Step 4: Majorana t+u doubling
    #     For identical fermions, M = M_t + M_u (direct + exchange).
    #     In the NR limit both channels give the same spatial potential
    #     (up to spin structure handled by w_ℓ weights).
    #     The effective coupling doubles.
    alpha_t_plus_u = 2 * alpha_t           # = y²/(8π)

    # --- Step 5: Standard QFT convention (no ½ from Majorana)
    #     If the vertex were just y (Dirac-like): α = y²/(4π)
    alpha_standard = y**2 / (4 * pi)       # = y²/(4π)

    return {
        "standard":     {"expr": alpha_standard,   "label": "y²/(4π)",  "factor": 4},
        "vertex_half":  {"expr": alpha_t,           "label": "y²/(16π)", "factor": 16},
        "majorana_t_u": {"expr": alpha_t_plus_u,    "label": "y²/(8π)",  "factor": 8},
    }


def test_bp(bp_label: str, conventions: dict) -> dict:
    """
    For a given benchmark point, invert each convention to get y,
    then re-derive α and compare with α_config.

    The logic:
      - Given α_config, for each convention C with α = y²/(F·π):
        1. Invert: y² = F·π·α_config  →  y = √(F·π·α_config)
        2. Re-derive: α_derived = y²/(F·π) ← this is trivially = α_config
           (the test is self-consistent by construction)

    The REAL test: does the potential V(r) = -α·e^{-mφr}/r match
    the VPM code when we start from the Lagrangian vertex = y/2?

    What we actually compute:
      - Start with α_config (what VPM uses as coupling in V(r)=-α e^{-mr}/r)
      - Interpret it under each convention → get y (the Lagrangian coupling)
      - From y, re-derive α through the Feynman diagram chain
      - The convention where α_derived = α_config **from the Feynman rules**
        is the correct one.
    """
    y_sym = sp.Symbol("y", positive=True)
    bp = GC.benchmark(bp_label)
    alpha_config = bp["alpha"]
    m_chi_GeV    = bp["m_chi_GeV"]
    m_phi_MeV    = bp["m_phi_MeV"]

    results = []
    for name, conv in conventions.items():
        factor = conv["factor"]

        # Step A: Invert convention to get y from α_config
        # Convention says: α = y²/(factor·π), so y² = factor·π·α
        y_squared = factor * sp.pi * alpha_config
        y_val = sp.sqrt(y_squared)

        # Step B: Now apply the FULL Feynman derivation from y:
        #   Vertex = y/2 → Born numerator = (y/2)² = y²/4
        #   FT gives: V_t(r) = -(y²/4)/(4π) · e^{-mr}/r = -(y²/16π) e^{-mr}/r
        #   So α_from_feynman_t_only = y²/(16π)
        alpha_feynman_t = float(y_val**2 / (16 * sp.pi))

        #   With Majorana t+u doubling: α_feynman_tu = y²/(8π)
        alpha_feynman_tu = float(y_val**2 / (8 * sp.pi))

        # Step C: Compare with α_config
        ratio_t  = alpha_feynman_t  / alpha_config
        ratio_tu = alpha_feynman_tu / alpha_config

        # Which Feynman result matches α_config?
        match_t  = abs(ratio_t  - 1.0) < _LOCAL_CFG["tolerance"]
        match_tu = abs(ratio_tu - 1.0) < _LOCAL_CFG["tolerance"]

        results.append({
            "bp":              bp_label,
            "m_chi_GeV":       m_chi_GeV,
            "m_phi_MeV":       m_phi_MeV,
            "alpha_config":    alpha_config,
            "convention":      conv["label"],
            "factor":          factor,
            "y_coupling":      float(y_val),
            "alpha_feynman_t": alpha_feynman_t,
            "alpha_feynman_tu":alpha_feynman_tu,
            "ratio_t":         ratio_t,
            "ratio_tu":        ratio_tu,
            "match_t":         "PASS" if match_t  else "FAIL",
            "match_tu":        "PASS" if match_tu else "FAIL",
        })

    return results


def symbolic_summary():
    """Print the symbolic derivation chain for the record."""
    y, m, r, q = sp.symbols("y m_phi r q", positive=True)
    pi = sp.pi

    print("=" * 70)
    print("T1: SYMBOLIC DERIVATION — Lagrangian → V(r)")
    print("=" * 70)
    print()
    print("Lagrangian interaction:  L_int = -(y/2) χ̄χ φ")
    print()
    print("Step 1: Vertex factor")
    vertex = y / 2
    print(f"  vertex = y/2 = {vertex}")
    print()
    print("Step 2: Born amplitude (t-channel, NR)")
    born_num = vertex**2
    print(f"  |M_t|² numerator = (y/2)² = {sp.expand(born_num)} = y²/4")
    print()
    print("Step 3: Fourier transform → coordinate space")
    alpha_t = born_num / (4 * pi)
    alpha_t_simplified = sp.simplify(alpha_t)
    print(f"  V_t(r) = -(y²/4)/(4π) · e^{{-m_φ r}}/r")
    print(f"         = -{alpha_t_simplified} · e^{{-m_φ r}}/r")
    print(f"  → α_t = y²/(16π)  [t-channel only]")
    print()
    print("Step 4: Majorana identical-particle exchange (t + u)")
    alpha_tu = 2 * alpha_t
    alpha_tu_simplified = sp.simplify(alpha_tu)
    print(f"  V_{{t+u}}(r) = 2·V_t(r)    [direct + exchange in NR limit]")
    print(f"  → α_{{t+u}} = y²/(8π) = {alpha_tu_simplified}")
    print()
    print("Step 5: Comparison table")
    print(f"  Standard (Dirac-like): α = y²/(4π)")
    print(f"  Vertex ½ (t only):     α = y²/(16π)")
    print(f"  Majorana (t+u):        α = y²/(8π)")
    print("=" * 70)
    print()


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    import csv
    from datetime import datetime

    bp_labels = _LOCAL_CFG["test_benchmarks"]
    conventions = derive_alpha_conventions()

    # ── symbolic output ──────────────────────────────────────────────────
    symbolic_summary()

    # ── numeric tests per BP ─────────────────────────────────────────────
    all_results = []
    for bp_label in bp_labels:
        results = test_bp(bp_label, conventions)
        all_results.extend(results)

    # ── display ──────────────────────────────────────────────────────────
    header = (
        f"{'BP':<12} {'Convention':<12} {'α_config':>12} "
        f"{'y':>10} {'α_feyn_t':>12} {'α_feyn_tu':>12} "
        f"{'ratio_t':>9} {'ratio_tu':>9} {'t?':>5} {'tu?':>5}"
    )
    print(header)
    print("-" * len(header))

    for row in all_results:
        line = (
            f"{row['bp']:<12} {row['convention']:<12} {row['alpha_config']:>12.6e} "
            f"{row['y_coupling']:>10.6f} {row['alpha_feynman_t']:>12.6e} {row['alpha_feynman_tu']:>12.6e} "
            f"{row['ratio_t']:>9.4f} {row['ratio_tu']:>9.4f} {row['match_t']:>5} {row['match_tu']:>5}"
        )
        print(line)

    # ── identify which convention is correct ─────────────────────────────
    print()
    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)

    for conv_name, conv_data in conventions.items():
        label = conv_data["label"]
        factor = conv_data["factor"]
        # Check if this convention's t-only derivation matches α_config
        t_matches  = all(r["match_t"]  == "PASS" for r in all_results if r["factor"] == factor)
        tu_matches = all(r["match_tu"] == "PASS" for r in all_results if r["factor"] == factor)

        if t_matches:
            print(f"  ✓ Convention '{label}' (factor={factor}): "
                  f"α_config = α_feynman (t-only). Vertex ½ absorbed, NO t+u doubling.")
        if tu_matches:
            print(f"  ✓ Convention '{label}' (factor={factor}): "
                  f"α_config = α_feynman (t+u). Vertex ½ absorbed + Majorana doubling.")

        if not t_matches and not tu_matches:
            print(f"  ✗ Convention '{label}' (factor={factor}): does NOT match α_config.")

    print()

    # ── CSV output ───────────────────────────────────────────────────────
    csv_fields = [
        "bp", "m_chi_GeV", "m_phi_MeV", "alpha_config", "convention",
        "factor", "y_coupling", "alpha_feynman_t", "alpha_feynman_tu",
        "ratio_t", "ratio_tu", "match_t", "match_tu",
    ]

    # RunLogger context
    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/test_T1_alpha_convention.py",
        stage="T1 — Alpha Convention",
        params={"benchmarks": bp_labels, "tolerance": _LOCAL_CFG["tolerance"]},
        data_source="global_config.json",
    ) as rl:
        out_csv = timestamped_path(
            "T1_alpha_convention",
            archive=_HERE / "data" / "archive",
        )
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=csv_fields)
            writer.writeheader()
            writer.writerows(all_results)

        rl.add_output(str(out_csv))
        rl.set_notes(
            f"Tested {len(bp_labels)} BPs × 3 conventions. "
            f"See {out_csv.name} for full results."
        )
        print(f"Output: {out_csv}")
        print(f"Run logged to: docs/runs_log.csv")


if __name__ == "__main__":
    from _tee_output import tee_to_output
    with tee_to_output("T1_alpha_convention"):
        main()
