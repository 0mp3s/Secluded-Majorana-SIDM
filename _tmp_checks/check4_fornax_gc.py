#!/usr/bin/env python3
"""
Check 4: Fornax GC timing constraint (Read+2019).

Runs the existing predict_fornax_gc.py and compares σ/m at v ~ 12 km/s
against Read+2019 bound: σ/m < 1.5 cm^2/g at the Fornax GC velocity.
"""
import sys, os, math, time

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

from v22_raw_scan import sigma_T_vpm
from global_config import GC

# Read+2019 (1808.06634): GC orbital decay constrains σ/m at v ~ 10-20 km/s
# Conservative bound: σ/m < 1.5 cm^2/g at the internal velocity of Fornax
READ2019_BOUND = 1.5  # cm^2/g
GC_VELOCITIES = [10, 12, 15, 20]  # km/s — internal velocity of Fornax


def main():
    print("=" * 80)
    print("  CHECK 4: Fornax GC Timing — Read+2019 Constraint")
    print(f"  Bound: σ/m < {READ2019_BOUND} cm²/g at v ~ 10-20 km/s")
    print("=" * 80)

    print("\n  Warming up JIT...")
    t0 = time.time()
    _ = sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)
    print(f"  JIT: {time.time()-t0:.1f}s\n")

    named_bps = GC.all_benchmarks()

    print(f"  {'BP':>12}", end="")
    for v in GC_VELOCITIES:
        print(f"  {'σ/m@'+str(v):>10}", end="")
    print(f"  {'Worst':>10} {'Margin':>10} {'Status':>8}")
    print("  " + "-" * (12 + 10*len(GC_VELOCITIES) + 30))

    all_pass = True
    for bp in named_bps:
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] * 1e-3
        alpha = bp['alpha']

        sms = []
        print(f"  {bp['label']:>12}", end="")
        for v in GC_VELOCITIES:
            sm = sigma_T_vpm(m_chi, m_phi, alpha, float(v))
            sms.append(sm)
            print(f"  {sm:10.4f}", end="")

        worst = max(sms)
        margin = READ2019_BOUND - worst
        ok = worst < READ2019_BOUND

        if not ok:
            all_pass = False

        print(f"  {worst:10.4f} {margin:+10.4f} {'PASS' if ok else 'FAIL':>8}")

    print()
    if all_pass:
        print("  ✓ ALL BPs SURVIVE Fornax GC constraint")
    else:
        print("  ✗ SOME BPs EXCEED Read+2019 bound!")

    # Also check from the relic-viable CSV
    print("\n  --- Relic-viable points: worst-case at v=12 km/s ---")
    from output_manager import get_latest
    import csv as _csv
    BP_CSV = str(get_latest("v31_true_viable_points"))
    with open(BP_CSV, newline='') as f:
        rows = list(_csv.DictReader(f))

    worst_sm = 0
    worst_bp = ""
    n_fail = 0
    for i, row in enumerate(rows):
        m_chi = float(row['m_chi_GeV'])
        m_phi = float(row['m_phi_MeV']) * 1e-3
        alpha = float(row['alpha'])
        sm12 = sigma_T_vpm(m_chi, m_phi, alpha, 12.0)
        if sm12 > worst_sm:
            worst_sm = sm12
            worst_bp = f"BP{i+1}"
        if sm12 > READ2019_BOUND:
            n_fail += 1

    print(f"  {len(rows)} points checked at v=12 km/s")
    print(f"  Worst: {worst_bp} with σ/m = {worst_sm:.4f} cm²/g (bound: {READ2019_BOUND})")
    print(f"  Failures: {n_fail}/{len(rows)}")

    return all_pass


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)
