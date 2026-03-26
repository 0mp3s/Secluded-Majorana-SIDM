#!/usr/bin/env python3
"""
Check 2: Unitarity bounds on all relic-viable BPs at 9 velocities.

Verifies that sigma_T(VPM) <= sigma_T(unitarity) = (l_max+1)^2 / kappa^2
for each BP at each velocity. This is a necessary condition: no partial
wave can have sin^2(delta_l) > 1.
"""
import sys, os, csv, math, time

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

from v22_raw_scan import sigma_T_vpm, C_KM_S
from output_manager import get_latest
from global_config import GC

_PC = GC.physical_constants()
GEV2_TO_CM2 = _PC["GEV2_to_cm2"]
GEV_IN_G = _PC["GeV_in_g"]

BP_CSV = str(get_latest("v31_true_viable_points"))
VELOCITIES = [12, 30, 60, 100, 200, 500, 1000, 2000, 4700]


def unitarity_limit_dimensionless(kappa, lam):
    """Maximum sigma_tilde allowed by unitarity: sum_l (2l+1) * w_l."""
    if kappa < 1e-15:
        return float('inf')
    x_max = 50.0 if kappa < 5 else (80.0 if kappa < 50 else 100.0)
    l_max = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)
    # Majorana weights: w=1 (even), w=3 (odd)
    sigma_max = 0.0
    for l in range(l_max + 1):
        w = 1.0 if l % 2 == 0 else 3.0
        sigma_max += w * (2 * l + 1)
    return sigma_max


def main():
    print("=" * 80)
    print("  CHECK 2: Unitarity Bounds — All BPs × 9 Velocities")
    print("=" * 80)

    bps = []
    with open(BP_CSV, newline='') as f:
        for i, row in enumerate(csv.DictReader(f)):
            bps.append({
                'label': f'BP{i+1}',
                'm_chi_GeV': float(row['m_chi_GeV']),
                'm_phi_MeV': float(row['m_phi_MeV']),
                'alpha': float(row['alpha']),
            })

    print(f"\n  Loaded {len(bps)} BPs, testing {len(VELOCITIES)} velocities each")
    print(f"  Total checks: {len(bps) * len(VELOCITIES)}")

    # JIT warmup
    print("  Warming up Numba JIT...")
    t0 = time.time()
    _ = sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)
    print(f"  JIT warm-up: {time.time()-t0:.1f}s\n")

    n_fail = 0
    n_total = 0
    worst_ratio = 0.0
    worst_info = ""

    for bp in bps:
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] * 1e-3  # GeV
        alpha = bp['alpha']
        lam = alpha * m_chi / m_phi

        for v in VELOCITIES:
            n_total += 1
            mu = m_chi / 2.0
            k = mu * (v / C_KM_S)
            kappa = k / m_phi

            if kappa < 1e-15:
                continue

            # sigma/m in cm^2/g
            sm = sigma_T_vpm(m_chi, m_phi, alpha, float(v))

            # Convert to dimensionless for unitarity comparison
            sigma_phys_cm2 = sm * (m_chi * GEV_IN_G)  # cm^2
            sigma_phys_gev = sigma_phys_cm2 / GEV2_TO_CM2  # GeV^-2
            sigma_tilde = sigma_phys_gev * (k**2 / (2.0 * math.pi))

            sigma_max = unitarity_limit_dimensionless(kappa, lam)
            ratio = sigma_tilde / sigma_max if sigma_max > 0 else 0

            if ratio > worst_ratio:
                worst_ratio = ratio
                worst_info = f"{bp['label']} @ {v} km/s"

            if sigma_tilde > sigma_max * 1.01:  # 1% tolerance
                n_fail += 1
                print(f"  FAIL: {bp['label']} @ {v} km/s: "
                      f"sigma_tilde={sigma_tilde:.4e}, max={sigma_max:.4e}, ratio={ratio:.4f}")

    print(f"\n  Results: {n_total} checks, {n_fail} failures")
    print(f"  Worst ratio (sigma/sigma_max): {worst_ratio:.6f} ({worst_info})")

    if n_fail == 0:
        print("\n  ✓ ALL UNITARITY BOUNDS SATISFIED")
    else:
        print(f"\n  ✗ {n_fail} UNITARITY VIOLATIONS FOUND")

    return n_fail == 0


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)
