#!/usr/bin/env python3
"""Quick comparison: old (symmetric) vs new (one-sided upper limit) chi2 for 17 relic BPs."""
# === path setup (auto-generated) ================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# =================================================================

import sys, os, math, csv
from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm
from output_manager import get_latest

sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)  # JIT warmup

_CFG = load_config(__file__)

OBSERVATIONS = GC.observations_as_tuples()
OBS_VELOCITIES = sorted(set(v for _, v, *_ in OBSERVATIONS))


def chi2_old(sigma_at_v):
    c2 = 0.0
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        theory = sigma_at_v.get(v, 0.0)
        if theory >= central:
            sigma = hi - central if hi > central else 0.5 * central
        else:
            sigma = central - lo if central > lo else 0.5 * central
        if sigma <= 0:
            sigma = 0.5 * max(central, 0.01)
        c2 += ((theory - central) / sigma) ** 2
    return c2


def chi2_new(sigma_at_v):
    c2 = 0.0
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        theory = sigma_at_v.get(v, 0.0)
        if lo == 0.0 and theory <= hi:
            continue
        if theory >= central:
            sigma = hi - central if hi > central else 0.5 * central
        else:
            sigma = central - lo if central > lo else 0.5 * central
        if sigma <= 0:
            sigma = 0.5 * max(central, 0.01)
        c2 += ((theory - central) / sigma) ** 2
    return c2


def pulls_detail(sigma_at_v):
    """Return per-system pulls for the new (one-sided) scheme."""
    pulls = []
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        theory = sigma_at_v.get(v, 0.0)
        if lo == 0.0 and theory <= hi:
            pulls.append((name, v, theory, central, 0.0))
            continue
        if theory >= central:
            sigma = hi - central if hi > central else 0.5 * central
        else:
            sigma = central - lo if central > lo else 0.5 * central
        if sigma <= 0:
            sigma = 0.5 * max(central, 0.01)
        pull = (theory - central) / sigma
        pulls.append((name, v, theory, central, pull))
    return pulls


if __name__ == '__main__':
    relic_csv = _CFG.get("benchmark_csv") or str(get_latest("v31_true_viable_points"))
    _DIR = os.path.dirname(os.path.abspath(__file__))
    if not os.path.isabs(relic_csv):
        relic_csv = os.path.normpath(os.path.join(_DIR, relic_csv))
    relic_points = []
    with open(relic_csv) as f:
        for row in csv.DictReader(f):
            mc = float(row['m_chi_GeV'])
            mp_gev = float(row['m_phi_MeV']) / 1000.0
            al = float(row['alpha'])
            om = float(row['omega_h2'])
            relic_points.append((mc, mp_gev, al, om))

    ndof_old = 13 - 3
    ndof_new = 11 - 3  # 2 upper limits contribute 0 when satisfied

    print("=" * 100)
    print("  One-sided upper limit fix: Bullet Cluster + Harvey+15")
    print("  Old: 13 systems all two-sided symmetric chi2")
    print("  New: 11 measurements + 2 one-sided upper limits (chi2=0 when theory < limit)")
    print("=" * 100)
    print()
    print(f"{'#':>3} {'m_chi':>8} {'m_phi':>8} {'alpha':>12} "
          f"{'chi2_old':>10} {'c2/n_old':>10} {'chi2_new':>10} {'c2/n_new':>10} {'delta':>8}")
    print("-" * 100)

    for i, (mc, mp, al, om) in enumerate(relic_points):
        sigma_at_v = {}
        for v in OBS_VELOCITIES:
            sigma_at_v[v] = sigma_T_vpm(mc, mp, al, float(v))
        c2o = chi2_old(sigma_at_v)
        c2n = chi2_new(sigma_at_v)
        print(f"{i+1:>3} {mc:>8.2f} {mp*1e3:>8.4f} {al:>12.4e} "
              f"{c2o:>10.4f} {c2o/ndof_old:>10.4f} {c2n:>10.4f} {c2n/ndof_new:>10.4f} "
              f"{c2o-c2n:>8.4f}")

    # Detailed pulls for best relic BP
    print()
    print("=" * 80)
    print("  Detailed pulls for best relic BP (one-sided scheme)")
    print("=" * 80)

    # Find best
    best_idx = None
    best_c2 = 1e30
    all_sigma = []
    for i, (mc, mp, al, om) in enumerate(relic_points):
        sigma_at_v = {}
        for v in OBS_VELOCITIES:
            sigma_at_v[v] = sigma_T_vpm(mc, mp, al, float(v))
        all_sigma.append(sigma_at_v)
        c2 = chi2_new(sigma_at_v)
        if c2 < best_c2:
            best_c2 = c2
            best_idx = i

    mc, mp, al, om = relic_points[best_idx]
    print(f"\n  BP{best_idx+1}: m_chi={mc:.2f}, m_phi={mp*1e3:.4f} MeV, alpha={al:.4e}")
    print(f"  chi2_new = {best_c2:.4f}, chi2_new/ndof = {best_c2/ndof_new:.4f}")
    print()
    print(f"  {'System':<25} {'v':>6} {'Theory':>10} {'Central':>10} {'Pull':>8}")
    print(f"  {'-'*65}")
    for name, v, theory, central, pull in pulls_detail(all_sigma[best_idx]):
        print(f"  {name:<25} {v:>6} {theory:>10.4f} {central:>10.2f} {pull:>+8.2f}")


if __name__ == '__main__':
    try:
        import sys as _sys, os as _os
        _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'core'))
        from tg_notify import notify
        notify("\u2705 bullet_onesided done!")
    except Exception:
        pass
