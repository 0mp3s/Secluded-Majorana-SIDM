#!/usr/bin/env python3
"""
model_validations/cp_separation/cp_separation_table.py
======================================================
§7.1 — CP-separation signature (unique prediction of Mixed Majorana SIDM).

The key insight: relic density depends on α_s × α_p (s-wave annihilation),
but SIDM σ_T/m depends on α_s only (elastic scattering). This produces a
1D band of viable (α_s, α_p) pairs where α_s/α_p can span ~1.7 decades
while maintaining identical SIDM phenomenology at each α_s.

This script:
  1. Reads the 13 viable points from condition2 output
  2. Computes additional velocities (v = 10, 50, 200, 500 km/s)
  3. Builds the CP-separation table for the paper
  4. Saves CSV + formatted table

Output: model_validations/cp_separation/output/
"""
import sys, os, csv, math, json
import numpy as np

# ── path bootstrap ──
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_SCRIPT_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ── load config ──
with open(os.path.join(_SCRIPT_DIR, 'config.json')) as f:
    _cfg = json.load(f)

from global_config import GC
_gc_bp1 = GC.benchmark("BP1")
M_CHI = _gc_bp1['m_chi_GeV']          # GeV
M_PHI = _gc_bp1['m_phi_MeV'] * 1e-3   # GeV (VPM needs GeV)
RELIC_PRODUCT = _cfg['BP1']['relic_product']  # α_s × α_p

VIABLE_ALPHA_S = _cfg['viable_alpha_s_BP1']
VELOCITIES = _cfg['velocities_km_s']

# ── output directory ──
OUT_DIR = os.path.join(_SCRIPT_DIR, _cfg.get('output_dir', 'output'))
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    print("=" * 80)
    print("§7.1 — CP-SEPARATION TABLE")
    print("=" * 80)
    print(f"  m_χ = {M_CHI} GeV,  m_φ = {M_PHI*1e3:.2f} MeV")
    print(f"  Relic constraint: α_s × α_p = {RELIC_PRODUCT:.6e}")
    print(f"  Viable points: {len(VIABLE_ALPHA_S)}")
    print(f"  Velocities: {VELOCITIES} km/s")
    print()

    rows = []
    for alpha_s in VIABLE_ALPHA_S:
        alpha_p = RELIC_PRODUCT / alpha_s
        ratio = alpha_s / alpha_p
        lam = alpha_s * M_CHI / M_PHI  # λ = α_s m_χ / m_φ (NO factor 2)

        sigmas = {}
        for v in VELOCITIES:
            sigmas[v] = sigma_T_vpm(M_CHI, M_PHI, alpha_s, float(v))

        rows.append({
            'alpha_s': alpha_s,
            'alpha_p': alpha_p,
            'ratio': ratio,
            'lambda': lam,
            **{f'sigma_{v}': sigmas[v] for v in VELOCITIES},
        })

    # ── Print formatted table ──
    header = (
        f"{'α_s':>12s}  {'α_p':>12s}  {'α_s/α_p':>10s}  {'λ':>8s}  " +
        "  ".join(f"{'σ/m(' + str(v) + ')':>10s}" for v in VELOCITIES)
    )
    print(header)
    print("-" * len(header))

    for row in rows:
        line = (
            f"{row['alpha_s']:12.4e}  {row['alpha_p']:12.4e}  "
            f"{row['ratio']:10.2f}  {row['lambda']:8.3f}  " +
            "  ".join(f"{row[f'sigma_{v}']:10.4f}" for v in VELOCITIES)
        )
        print(line)

    # ── Key observables ──
    print()
    print("=" * 80)
    print("KEY OBSERVABLES:")
    print("=" * 80)

    ratios = [r['ratio'] for r in rows]
    print(f"  α_s/α_p range:   {min(ratios):.1f} — {max(ratios):.1f}")
    print(f"  Dynamic range:   {max(ratios)/min(ratios):.1f}× ({math.log10(max(ratios)/min(ratios)):.2f} decades)")

    for v in [30, 100, 1000]:
        sigs = [r[f'sigma_{v}'] for r in rows]
        print(f"  σ/m({v:>4d}):      {min(sigs):.4f} — {max(sigs):.4f} cm²/g")

    lambdas = [r['lambda'] for r in rows]
    print(f"  λ range:         {min(lambdas):.2f} — {max(lambdas):.2f}")
    print(f"  λ < π = {math.pi:.2f}? Below first resonance: "
          f"{'ALL' if max(lambdas) < math.pi else 'NOT ALL'}")

    # ── Unique prediction: collider fingerprint ──
    print()
    print("UNIQUE PREDICTION — CP-separation collider signature:")
    print(f"  If α_s ≠ α_p (ratio ≠ 1), the dark sector has CP-violating structure.")
    print(f"  Mono-jet + MET cross section ∝ (y_s² + y_p²) = 4π(α_s + α_p)")
    print(f"  But relic Ωh² depends on α_s × α_p → independent constraints.")
    print(f"  Example: α_s/α_p = 200 → y_s ≈ 0.26, y_p ≈ 0.018 → pseudo-scalar coupling")
    print(f"           suppressed by factor ~14 relative to scalar.")
    print(f"  A collider measuring both couplings independently can test this.")

    # ── Save CSV ──
    csv_path = os.path.join(OUT_DIR, 'cp_separation_table.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['alpha_s', 'alpha_p', 'alpha_s_over_alpha_p', 'lambda'] +
                   [f'sigma_m_{v}' for v in VELOCITIES])
        for row in rows:
            w.writerow([
                f"{row['alpha_s']:.6e}",
                f"{row['alpha_p']:.6e}",
                f"{row['ratio']:.3f}",
                f"{row['lambda']:.4f}",
            ] + [f"{row[f'sigma_{v}']:.6f}" for v in VELOCITIES])
    print(f"\n  CSV saved: {csv_path}")

    # ── Save markdown table ──
    md_path = os.path.join(OUT_DIR, 'cp_separation_table.md')
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write("# §7.1 — CP-Separation Table\n\n")
        f.write(f"m_χ = {M_CHI} GeV, m_φ = {M_PHI*1e3:.2f} MeV\n\n")
        f.write(f"Relic constraint: α_s × α_p = {RELIC_PRODUCT:.6e}\n\n")

        # Table header
        f.write("| α_s | α_p | α_s/α_p | λ | σ/m(30) | σ/m(100) | σ/m(1000) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for row in rows:
            f.write(f"| {row['alpha_s']:.4e} | {row['alpha_p']:.4e} | "
                    f"{row['ratio']:.1f} | {row['lambda']:.2f} | "
                    f"{row['sigma_30']:.3f} | {row['sigma_100']:.3f} | "
                    f"{row['sigma_1000']:.3f} |\n")

        f.write(f"\nα_s/α_p spans {min(ratios):.1f} — {max(ratios):.1f} "
                f"({math.log10(max(ratios)/min(ratios)):.2f} decades)\n")
    print(f"  Markdown saved: {md_path}")

    print("\n✓ §7.1 CP-separation DONE.")


if __name__ == '__main__':
    main()
