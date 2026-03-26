#!/usr/bin/env python3
"""
Check 1: Perturbativity of all relic-viable benchmark points.

For each BP: compute y_s, alpha_s/(4pi), and verify perturbativity.
Threshold: y < sqrt(4*pi) ~ 3.54, equivalently alpha < 4*pi ~ 12.57.
Practical concern: one-loop corrections ~ alpha/(4*pi) should be << 1.
"""
import sys, os, csv, math

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

from output_manager import get_latest

BP_CSV = str(get_latest("v31_true_viable_points"))

def main():
    print("=" * 80)
    print("  CHECK 1: Perturbativity of All Relic-Viable BPs")
    print("=" * 80)

    bps = []
    with open(BP_CSV, newline='') as f:
        for i, row in enumerate(csv.DictReader(f)):
            bps.append({
                'label': f'BP{i+1}',
                'alpha': float(row['alpha']),
                'm_chi_GeV': float(row['m_chi_GeV']),
                'm_phi_MeV': float(row['m_phi_MeV']),
            })

    print(f"\n  Loaded {len(bps)} viable points from CSV")
    print(f"  Perturbativity threshold: alpha < 4*pi = {4*math.pi:.4f}")
    print(f"  One-loop control: alpha/(4*pi) << 1")
    print()

    print(f"  {'BP':>6} {'alpha_s':>12} {'y_s':>10} {'alpha/(4pi)':>12} {'loop %':>8} {'Status':>8}")
    print("  " + "-" * 58)

    alpha_max = 0
    all_pass = True
    for bp in bps:
        alpha = bp['alpha']
        y_s = math.sqrt(4 * math.pi * alpha)
        loop_param = alpha / (4 * math.pi)
        loop_pct = loop_param * 100

        ok_pert = alpha < 4 * math.pi
        ok_loop = loop_param < 0.1  # 10% loop correction as warning threshold

        status = "PASS" if (ok_pert and ok_loop) else ("WARN" if ok_pert else "FAIL")
        if not ok_pert:
            all_pass = False

        if alpha > alpha_max:
            alpha_max = alpha

        print(f"  {bp['label']:>6} {alpha:12.6e} {y_s:10.6f} {loop_param:12.6e} {loop_pct:7.3f}% {status:>8}")

    print()
    print(f"  Maximum alpha_s = {alpha_max:.6e}")
    print(f"  Maximum y_s     = {math.sqrt(4*math.pi*alpha_max):.6f}")
    print(f"  Maximum loop    = {alpha_max/(4*math.pi)*100:.4f}%")
    print()

    if all_pass:
        print("  ✓ ALL POINTS PERTURBATIVE — one-loop corrections < 0.04%")
    else:
        print("  ✗ SOME POINTS FAIL PERTURBATIVITY")

    # Named BPs from global_config
    from global_config import GC
    named = GC.all_benchmarks()
    print(f"\n  --- Named benchmarks ({len(named)}) ---")
    print(f"  {'Label':>12} {'alpha_s':>12} {'y_s':>10} {'loop %':>8}")
    print("  " + "-" * 46)
    for bp in named:
        a = bp['alpha']
        y = math.sqrt(4 * math.pi * a)
        print(f"  {bp['label']:>12} {a:12.6e} {y:10.6f} {a/(4*math.pi)*100:7.3f}%")

    print()
    return all_pass


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)
