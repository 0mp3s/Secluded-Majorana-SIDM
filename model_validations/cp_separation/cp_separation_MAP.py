#!/usr/bin/env python3
"""
model_validations/cp_separation/cp_separation_MAP.py
====================================================
§7.1 (continued) — CP-separation table for MAP's mass point.

Same logic as cp_separation_table.py (BP1 masses), but using MAP's
(m_χ, m_φ) = (94.07 GeV, 11.10 MeV).  Same relic product α_s × α_p.

The viable band is found by scanning α_s ∈ [5e-4, 0.04] and keeping
points where σ/m(30 km/s) ∈ [0.1, 10] cm²/g AND σ/m(1000) < 1.

Output: model_validations/cp_separation/output/cp_separation_MAP.csv
"""
import sys, os, csv, math
import numpy as np

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
from v22_raw_scan import sigma_T_vpm

sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)  # JIT warmup

# ── MAP mass point ──
M_CHI = 94.07       # GeV
M_PHI = 11.10e-3    # GeV
RELIC_PRODUCT = 1.387474e-7

# SIDM viability window (same as condition2)
SIGMA_30_MIN  = 0.1    # cm²/g
SIGMA_30_MAX  = 10.0   # cm²/g
SIGMA_1000_MAX = 1.0   # cm²/g (Bullet Cluster / Harvey+15)

VELOCITIES = [10, 30, 50, 100, 200, 500, 1000]

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    print("=" * 80)
    print("§7.1 — CP-SEPARATION TABLE (MAP MASSES)")
    print("=" * 80)
    print(f"  m_χ = {M_CHI} GeV,  m_φ = {M_PHI*1e3:.2f} MeV")
    print(f"  Relic constraint: α_s × α_p = {RELIC_PRODUCT:.6e}")
    print(f"  λ = α_s × m_χ / m_φ = {M_CHI/M_PHI:.1f} × α_s")
    print()

    # Scan α_s in log-space
    alpha_s_scan = np.geomspace(5e-4, 4e-2, 500)

    viable = []
    for alpha_s in alpha_s_scan:
        alpha_p = RELIC_PRODUCT / alpha_s
        # perturbativity: α < 4π
        if alpha_s > 4 * math.pi or alpha_p > 4 * math.pi:
            continue

        lam = alpha_s * M_CHI / M_PHI

        sigmas = {}
        for v in VELOCITIES:
            sigmas[v] = sigma_T_vpm(M_CHI, M_PHI, alpha_s, float(v))

        # SIDM viability check
        if sigmas[30] < SIGMA_30_MIN or sigmas[30] > SIGMA_30_MAX:
            continue
        if sigmas[1000] > SIGMA_1000_MAX:
            continue

        viable.append({
            'alpha_s': float(alpha_s),
            'alpha_p': alpha_p,
            'ratio': alpha_s / alpha_p,
            'lambda': lam,
            **{f'sigma_{v}': sigmas[v] for v in VELOCITIES},
        })

    print(f"  Scanned: {len(alpha_s_scan)} points")
    print(f"  Viable:  {len(viable)} points")
    print()

    if not viable:
        print("  *** NO VIABLE POINTS FOUND ***")
        print("  This may mean MAP's α_s is already at the edge of the viable band,")
        print("  or that the VPM resonance structure creates gaps.")
        return

    # Print table
    header = (
        f"{'α_s':>12s}  {'α_p':>12s}  {'α_s/α_p':>10s}  {'λ':>8s}  " +
        "  ".join(f"{'σ/m(' + str(v) + ')':>10s}" for v in VELOCITIES)
    )
    print(header)
    print("-" * len(header))
    for row in viable:
        line = (
            f"{row['alpha_s']:12.4e}  {row['alpha_p']:12.4e}  "
            f"{row['ratio']:10.2f}  {row['lambda']:8.2f}  " +
            "  ".join(f"{row[f'sigma_{v}']:10.4f}" for v in VELOCITIES)
        )
        print(line)

    # Summary
    print()
    print("=" * 80)
    print("KEY OBSERVABLES (MAP masses):")
    print("=" * 80)
    ratios = [r['ratio'] for r in viable]
    lambdas = [r['lambda'] for r in viable]
    print(f"  α_s/α_p range:  {min(ratios):.1f} — {max(ratios):.1f}")
    print(f"  Dynamic range:  {max(ratios)/min(ratios):.1f}× ({math.log10(max(ratios)/min(ratios)):.2f} decades)")
    print(f"  λ range:        {min(lambdas):.1f} — {max(lambdas):.1f}")
    print(f"  All λ > π?      {min(lambdas) > math.pi}")

    for v in [30, 100, 1000]:
        sigs = [r[f'sigma_{v}'] for r in viable]
        print(f"  σ/m({v:>4d}):     {min(sigs):.4f} — {max(sigs):.4f} cm²/g")

    # Compare with MAP's own α_s
    map_alpha_s = 5.734e-3
    map_lam = map_alpha_s * M_CHI / M_PHI
    map_in_band = any(abs(r['alpha_s'] - map_alpha_s)/map_alpha_s < 0.05 for r in viable)
    print(f"\n  MAP benchmark (α_s = {map_alpha_s}): λ = {map_lam:.1f}")
    print(f"  MAP in viable band? {map_in_band}")

    # CSV
    csv_path = os.path.join(OUT_DIR, 'cp_separation_MAP.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['alpha_s', 'alpha_p', 'alpha_s_over_alpha_p', 'lambda'] +
                   [f'sigma_m_{v}' for v in VELOCITIES])
        for row in viable:
            w.writerow([
                f"{row['alpha_s']:.6e}",
                f"{row['alpha_p']:.6e}",
                f"{row['ratio']:.3f}",
                f"{row['lambda']:.4f}",
            ] + [f"{row[f'sigma_{v}']:.6f}" for v in VELOCITIES])
    print(f"\n  CSV saved: {csv_path}")

    # Markdown summary
    md_path = os.path.join(OUT_DIR, 'cp_separation_MAP.md')
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write("# §7.1 — CP-Separation Table (MAP Masses)\n\n")
        f.write(f"m_χ = {M_CHI} GeV, m_φ = {M_PHI*1e3:.2f} MeV\n\n")
        f.write(f"Relic constraint: α_s × α_p = {RELIC_PRODUCT:.6e}\n\n")
        f.write(f"Viable points: {len(viable)} (from {len(alpha_s_scan)} scanned)\n\n")

        f.write("| α_s | α_p | α_s/α_p | λ | σ/m(30) | σ/m(100) | σ/m(1000) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for row in viable:
            f.write(f"| {row['alpha_s']:.4e} | {row['alpha_p']:.4e} | "
                    f"{row['ratio']:.1f} | {row['lambda']:.1f} | "
                    f"{row['sigma_30']:.3f} | {row['sigma_100']:.3f} | "
                    f"{row['sigma_1000']:.3f} |\n")

        f.write(f"\n## Summary\n")
        f.write(f"- α_s/α_p: {min(ratios):.1f} — {max(ratios):.1f}\n")
        f.write(f"- Dynamic range: {math.log10(max(ratios)/min(ratios)):.2f} decades\n")
        f.write(f"- λ range: {min(lambdas):.1f} — {max(lambdas):.1f}\n")
        f.write(f"- MAP benchmark in band: {map_in_band}\n")
    print(f"  Markdown saved: {md_path}")

    print(f"\n✓ MAP CP-separation scan DONE.")


if __name__ == '__main__':
    main()
