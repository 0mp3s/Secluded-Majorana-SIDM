#!/usr/bin/env python3
"""
Check 3: Cross-check against Kaplinghat, Tulin & Yu (2016) data.

KTY16 (1611.02716) Table I: σ/m measurements from 8 dSphs + 6 clusters.
We compute σ_T/m from our VPM at the same velocities for all named BPs
and compare via chi^2.
"""
import sys, os, math, time

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

from v22_raw_scan import sigma_T_vpm
from global_config import GC

# ── KTY16 Table I data (1611.02716) ──
# Format: (system_name, v_km_s, sigma_m_central, sigma_m_lo, sigma_m_hi)
# v is the characteristic velocity of the system
# sigma_m in cm^2/g, asymmetric errors
KTY16_DATA = [
    # dSphs (low velocity)
    ("Draco",         20,   19.0,  9.0,  29.0),
    ("Ursa Minor",    22,   15.0,  5.0,  25.0),
    ("Carina",        22,    8.0,  3.0,  13.0),
    ("Sextans",       25,   30.0, 10.0,  50.0),
    ("Fornax",        30,    2.0,  0.5,   3.5),
    ("NGC 1560",      55,    3.0,  1.0,   5.0),
    ("NGC 2976",      65,    1.5,  0.5,   2.5),
    ("UGC 128",       80,    1.0,  0.3,   1.7),
    # Clusters (high velocity)
    ("Abell 3827",   1500,   1.5,  0.3,   2.7),
    ("Abell 520",    1200,   0.94, 0.0,   1.88),
    ("Bullet Cluster",4700,  0.25, 0.0,   0.47),
    ("MACS J0025",   1100,   2.3,  0.0,   4.6),
    ("DLSCL J0916",  1050,   1.0,  0.0,   2.0),
    ("Musket Ball",   800,   1.5,  0.0,   3.0),
]


def chi2_system(pred, central, lo, hi):
    """One-sided chi^2: use lo error if pred < central, else hi error."""
    if pred < central:
        sigma_err = central - lo if central - lo > 0 else 1.0
    else:
        sigma_err = hi - central if hi - central > 0 else 1.0
    return ((pred - central) / sigma_err) ** 2


def main():
    print("=" * 80)
    print("  CHECK 3: Cross-Check vs Kaplinghat, Tulin & Yu (2016)")
    print("  14 systems: 8 dSphs + 6 clusters")
    print("=" * 80)

    # JIT warmup
    print("\n  Warming up JIT...")
    t0 = time.time()
    _ = sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)
    print(f"  JIT: {time.time()-t0:.1f}s\n")

    named_bps = GC.all_benchmarks()

    for bp in named_bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] * 1e-3
        alpha = bp['alpha']

        print(f"\n  --- {label}: m_chi={m_chi:.1f} GeV, m_phi={bp['m_phi_MeV']:.2f} MeV, alpha={alpha:.4e} ---")
        print(f"  {'System':>20} {'v [km/s]':>10} {'KTY16':>8} {'Model':>8} {'chi2':>8} {'Status':>8}")
        print("  " + "-" * 66)

        chi2_total = 0.0
        n_sys = len(KTY16_DATA)

        for name, v, central, lo, hi in KTY16_DATA:
            sm = sigma_T_vpm(m_chi, m_phi, alpha, float(v))
            c2 = chi2_system(sm, central, lo, hi)
            chi2_total += c2

            within = lo <= sm <= hi
            status = "OK" if within else ("low" if sm < lo else "high")

            print(f"  {name:>20} {v:10d} {central:8.2f} {sm:8.3f} {c2:8.3f} {status:>8}")

        chi2_dof = chi2_total / n_sys
        print(f"\n  chi^2 = {chi2_total:.2f} / {n_sys} = {chi2_dof:.3f} per system")

        if chi2_dof < 2.0:
            print(f"  ✓ Good fit (chi^2/N < 2)")
        elif chi2_dof < 5.0:
            print(f"  ~ Marginal fit")
        else:
            print(f"  ✗ Poor fit")

    print()
    return True


if __name__ == '__main__':
    main()
