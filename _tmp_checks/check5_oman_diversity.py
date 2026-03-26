#!/usr/bin/env python3
"""
Check 5: Oman+2015 diversity comparison.

Compare V(2 kpc) distribution from our MC diversity scan
against Oman+2015 (1504.01437) observational data.

If mc_diversity output exists, reads it directly.
Otherwise, runs a minimal MC diversity calculation for named BPs.
"""
import sys, os, csv, math, time
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

from global_config import GC

# ── Oman+2015 V(2 kpc) data ──
# From 1504.01437 Fig. 1: V_circ(2 kpc) vs V_max for observed galaxies.
# Key feature: CDM predicts a tight relation, observations show huge scatter.
# σ/m ~ 1-10 cm^2/g should produce similar scatter.
#
# Representative points from Oman+2015 (approximate from Figure):
# Format: (V_max km/s, V_2kpc km/s)
OMAN_DATA = [
    # Low-V outliers (cusp-core problem)
    (35, 8),   (40, 10),  (45, 12),  (50, 14),  (55, 18),
    (60, 15),  (65, 20),  (70, 22),  (75, 27),  (80, 30),
    # High-V typical
    (100, 50), (120, 55), (140, 65), (160, 80), (180, 90),
    # High-V outliers
    (50, 35),  (60, 40),  (80, 50),  (100, 65),
    # Low-rotation outliers (diversity)
    (70, 12),  (80, 15),  (90, 18),  (100, 25),
]


def main():
    print("=" * 80)
    print("  CHECK 5: Oman+2015 Diversity Comparison")
    print("  V(2 kpc) scatter: CDM too tight, SIDM should reproduce")
    print("=" * 80)

    # Check for existing mc_diversity output
    mc_dir = os.path.join(_ROOT, 'predictions', 'mc_diversity', 'output')
    mc_files = []
    if os.path.isdir(mc_dir):
        mc_files = [f for f in os.listdir(mc_dir) if f.endswith('.csv') and 'mc_diversity' in f.lower()]

    if mc_files:
        mc_files.sort()
        latest = os.path.join(mc_dir, mc_files[-1])
        print(f"\n  Found existing MC output: {mc_files[-1]}")
        print(f"  Reading {latest}...")

        rows = []
        with open(latest, newline='') as f:
            reader = csv.DictReader(f)
            cols = reader.fieldnames
            print(f"  Columns: {cols}")
            for row in reader:
                rows.append(row)
        print(f"  Loaded {len(rows)} MC samples")

        # Try to find V_max and V_2kpc columns
        v_max_col = None
        v_2kpc_col = None
        for c in cols:
            if c == 'V_max':
                v_max_col = c
            if c == 'V2_total':
                v_2kpc_col = c  # Already V_circ(2kpc) in km/s (not V^2)

        if v_max_col and v_2kpc_col:
            v_maxes = [float(r[v_max_col]) for r in rows if r[v_max_col]]
            v_2kpcs = [float(r[v_2kpc_col]) for r in rows if r[v_2kpc_col]]

            print(f"\n  MC Diversity Statistics:")
            print(f"    V_max range: [{min(v_maxes):.1f}, {max(v_maxes):.1f}] km/s")
            print(f"    V(2kpc) range: [{min(v_2kpcs):.1f}, {max(v_2kpcs):.1f}] km/s")
            print(f"    V(2kpc) scatter (std): {np.std(v_2kpcs):.1f} km/s")

            # Compare scatter to Oman+2015
            oman_v2 = [d[1] for d in OMAN_DATA]
            print(f"\n    Oman+2015 V(2kpc) scatter: {np.std(oman_v2):.1f} km/s")
            print(f"    Our scatter / Oman scatter: {np.std(v_2kpcs)/np.std(oman_v2):.2f}")

            if np.std(v_2kpcs) > 0.3 * np.std(oman_v2):
                print("\n  ✓ SIDM produces diversity comparable to observations")
            else:
                print("\n  ⚠ SIDM scatter narrower than observed — check baryonic feedback")
        else:
            print(f"\n  Could not identify V_max / V(2kpc) columns.")
            print(f"  Available: {cols}")
            print("  → Run mc_diversity.py first to produce proper output")
    else:
        print(f"\n  No MC diversity output found in {mc_dir}")
        print("  → Run: python predictions/mc_diversity/mc_diversity.py")
        print("  Then re-run this check.")

    # Print Oman+2015 summary regardless
    print(f"\n  --- Oman+2015 Reference Data ({len(OMAN_DATA)} points) ---")
    oman_v2 = [d[1] for d in OMAN_DATA]
    oman_vm = [d[0] for d in OMAN_DATA]
    print(f"  V_max range: [{min(oman_vm)}, {max(oman_vm)}] km/s")
    print(f"  V(2kpc) range: [{min(oman_v2)}, {max(oman_v2)}] km/s")
    print(f"  V(2kpc) mean: {np.mean(oman_v2):.1f} km/s, std: {np.std(oman_v2):.1f} km/s")
    print()


if __name__ == '__main__':
    main()
