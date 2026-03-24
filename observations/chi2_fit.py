#!/usr/bin/env python3
"""
V10 — v34_chi2_fit.py
======================
χ² fit of our Majorana-scalar SIDM model to REAL astrophysical data.

Full VPM solver on all 13 observational velocities.
Parallel: multiprocessing across all CPU cores.

Uses existing scan results (80,142 viable points from v22 + 17 relic
BPs from v31) and evaluates χ² against 13 observational systems.

Data: Kaplinghat+16, Kamada+17, Randall+08, Harvey+15, Elbert+15.
Output: v34_chi2_fit.png, v34_output.txt, v34_results.csv
"""
# === path setup (auto-generated) ================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# =================================================================

import sys, os, time, math
import numpy as np
import csv
import multiprocessing as mp_lib
from functools import partial
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from config_loader import load_config

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)
os.environ['PYTHONIOENCODING'] = 'utf-8'

_DIR = os.path.dirname(os.path.abspath(__file__))
if _DIR not in sys.path:
    sys.path.insert(0, _DIR)

# ================================================================
#  OBSERVATIONAL DATA  — loaded from config.json if available
# ================================================================
_CFG = load_config(__file__)

_DEFAULT_OBS = [
    ("Draco dSph",          12,   0.6,  0.1,  2.0,  "KTY16"),
    ("Fornax dSph",         12,   0.8,  0.2,  3.0,  "KTY16"),
    ("NGC 2976",            60,   2.0,  0.5,  5.0,  "KTY16"),
    ("NGC 1560",            55,   3.0,  1.0,  8.0,  "KTY16"),
    ("IC 2574",             50,   1.5,  0.3,  5.0,  "KTY16"),
    ("NGC 720 (group)",    250,   0.5,  0.1,  1.5,  "KTY16"),
    ("NGC 1332 (group)",   280,   0.3,  0.05, 1.0,  "KTY16"),
    ("Abell 611",         1200,   0.1,  0.02, 0.3,  "KTY16"),
    ("Abell 2537",        1100,   0.15, 0.03, 0.4,  "KTY16"),
    ("Diverse RC band",     40,   3.0,  0.5,  10.0, "KKPY17"),
    ("Bullet Cluster",    4700,   0.7,  0.0,  1.25, "Randall+08"),
    ("72 cluster mergers", 1000,  0.2,  0.0,  0.47, "Harvey+15"),
    ("TBTF dwarfs",         30,   1.0,  0.5,  5.0,  "Elbert+15"),
]

OBSERVATIONS = [tuple(o) for o in _CFG.get("observations", _DEFAULT_OBS)]

OBS_VELOCITIES = sorted(set(v for _, v, *_ in OBSERVATIONS))


def compute_chi2(sigma_at_v):
    chi2 = 0.0
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        theory = sigma_at_v.get(v, 0.0)
        if theory >= central:
            sigma = hi - central if hi > central else 0.5 * central
        else:
            sigma = central - lo if central > lo else 0.5 * central
        if sigma <= 0:
            sigma = 0.5 * max(central, 0.01)
        chi2 += ((theory - central) / sigma) ** 2
    return chi2


# ================================================================
#  Worker function for multiprocessing
# ================================================================
_solver = None

def _init_worker():
    """Each worker imports and warms up the JIT solver."""
    global _solver
    import sys, os
    d = os.path.dirname(os.path.abspath(__file__))
    if d not in sys.path:
        sys.path.insert(0, d)
    from v22_raw_scan import sigma_T_vpm
    _solver = sigma_T_vpm
    # warm up JIT in this process
    _solver(20.0, 10e-3, 1e-3, 100.0)


def _eval_one(args):
    """Evaluate one (mc, mp_GeV, alpha) → (chi2, mc, mp, al, {v: sigma})."""
    mc, mp, al = args
    sigma_at_v = {}
    for v in OBS_VELOCITIES:
        try:
            s = _solver(mc, mp, al, float(v))
            if math.isnan(s) or math.isinf(s):
                return None
            sigma_at_v[v] = s
        except Exception:
            return None
    c2 = compute_chi2(sigma_at_v)
    return (c2, mc, mp, al, sigma_at_v)


def run():
    t0 = time.time()
    hdr = "=" * 72
    ncpu = os.cpu_count() or 4
    nworkers = _CFG.get("n_workers", max(1, ncpu - 2))

    print(hdr)
    print("  V10 — v34 chi2 Fit to Observational Data")
    print(f"  {len(OBSERVATIONS)} systems x scan points + relic BPs")
    print(f"  Parallel: {nworkers} workers on {ncpu} cores")
    print(hdr)

    # ============================================================
    # LOAD EXISTING DATA  (paths overridable via config.json)
    # ============================================================
    raw_csv_path = _CFG.get("scan_data_csv", os.path.join(DATA_DIR, 'all_viable_raw_v8.csv'))
    if not os.path.isabs(raw_csv_path):
        raw_csv_path = os.path.normpath(os.path.join(_DIR, raw_csv_path))
    raw_points = []
    with open(raw_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mc = float(row['m_chi_GeV'])
            mp_gev = float(row['m_phi_GeV'])   # keep GeV for solver
            al = float(row['alpha'])
            raw_points.append((mc, mp_gev, al))
    print(f"  Loaded {len(raw_points)} raw viable points from scan")

    relic_csv_path = _CFG.get("relic_bp_csv", os.path.join(DATA_DIR, 'v31_true_viable_points.csv'))
    if not os.path.isabs(relic_csv_path):
        relic_csv_path = os.path.normpath(os.path.join(_DIR, relic_csv_path))
    relic_points = []
    with open(relic_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mc = float(row['m_chi_GeV'])
            mp_gev = float(row['m_phi_MeV']) / 1000.0  # MeV -> GeV
            al = float(row['alpha'])
            om = float(row['omega_h2'])
            relic_points.append((mc, mp_gev, al, om))
    print(f"  Loaded {len(relic_points)} relic-viable BPs from v31")

    # ============================================================
    # SAMPLE (every N-th for ~5000 points)
    # ============================================================
    step = max(1, len(raw_points) // 5000)
    sample = raw_points[::step]
    total_pts = len(sample) + len(relic_points)
    print(f"\n  Evaluating {len(sample)} sampled + {len(relic_points)} relic = {total_pts} points")
    print(f"  {len(OBS_VELOCITIES)} unique velocities => {total_pts * len(OBS_VELOCITIES)} VPM calls")
    print(f"  Using {nworkers} parallel workers...\n", flush=True)

    # ============================================================
    # PARALLEL EVALUATION (pool of workers)
    # ============================================================
    all_tasks = list(sample) + [(mc, mp, al) for mc, mp, al, _ in relic_points]

    results = []
    relic_results = []
    n_sample = len(sample)
    done = 0

    with mp_lib.Pool(processes=nworkers, initializer=_init_worker) as pool:
        for i, res in enumerate(pool.imap_unordered(_eval_one, all_tasks, chunksize=32)):
            done += 1
            if done % 500 == 0 or done == total_pts:
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total_pts - done) / rate if rate > 0 else 0
                print(f"    {done}/{total_pts}  [{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining, {rate:.1f} pts/s]", flush=True)
            if res is None:
                continue
            results.append(res)

    elapsed_scan = time.time() - t0
    print(f"\n  Scan done: {len(results)} valid / {total_pts} total in {elapsed_scan:.1f}s")

    # Separate relic BPs (last 17 tasks)
    # Re-evaluate relic BPs from results (they're mixed in with imap_unordered)
    # Easier: just tag them. Let's collect relic separately.
    relic_set = set()
    for mc, mp, al, om in relic_points:
        relic_set.add((round(mc, 6), round(mp, 10), round(al, 10)))

    relic_results = []
    free_results = []
    for r in results:
        c2, mc, mp, al, sv = r
        key = (round(mc, 6), round(mp, 10), round(al, 10))
        if key in relic_set:
            # find omega
            for mc2, mp2, al2, om2 in relic_points:
                if abs(mc2 - mc) < 1e-4 and abs(mp2 - mp) < 1e-10 and abs(al2 - al) < 1e-10:
                    relic_results.append((c2, mc, mp, al, om2, sv))
                    break
        else:
            free_results.append(r)

    free_results.sort(key=lambda x: x[0])
    relic_results.sort(key=lambda x: x[0])

    # ============================================================
    # Fine-scan around top 50 (serial, fast — only 400 points)
    # ============================================================
    if len(free_results) >= 20:
        print(f"  Fine-scanning around top 50...", flush=True)
        fine_tasks = []
        for c2, mc0, mp0, al0, _ in free_results[:50]:
            for af in [0.8, 0.85, 0.9, 0.95, 1.05, 1.1, 1.15, 1.2]:
                fine_tasks.append((mc0, mp0, al0 * af))

        with mp_lib.Pool(processes=nworkers, initializer=_init_worker) as pool:
            for res in pool.imap_unordered(_eval_one, fine_tasks, chunksize=16):
                if res is not None:
                    free_results.append(res)
        free_results.sort(key=lambda x: x[0])
        print(f"  Fine-scan done: {len(fine_tasks)} extra points evaluated")

    ndof = len(OBSERVATIONS) - 3

    # ============================================================
    # WRITE CSV — all results
    # ============================================================
    csv_out = os.path.join(_DIR, 'v34_results.csv')
    vel_cols = [f"sigma_m_{v}" for v in OBS_VELOCITIES]
    with open(csv_out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['chi2', 'chi2_dof', 'm_chi_GeV', 'm_phi_MeV', 'alpha', 'lambda', 'is_relic', 'omega_h2'] + vel_cols)
        for c2, mc, mp, al, sv in free_results:
            lam = al * mc / mp if mp > 0 else 0
            row = [f"{c2:.4f}", f"{c2/ndof:.4f}", f"{mc:.6f}", f"{mp*1e3:.6f}",
                   f"{al:.6e}", f"{lam:.4f}", "0", ""]
            row += [f"{sv.get(v, 0):.6e}" for v in OBS_VELOCITIES]
            w.writerow(row)
        for c2, mc, mp, al, om, sv in relic_results:
            lam = al * mc / mp if mp > 0 else 0
            row = [f"{c2:.4f}", f"{c2/ndof:.4f}", f"{mc:.6f}", f"{mp*1e3:.6f}",
                   f"{al:.6e}", f"{lam:.4f}", "1", f"{om:.6f}"]
            row += [f"{sv.get(v, 0):.6e}" for v in OBS_VELOCITIES]
            w.writerow(row)
    print(f"\n  Saved CSV -> {csv_out}  ({len(free_results) + len(relic_results)} rows)")

    # ============================================================
    # RESULTS: FREE BEST FIT
    # ============================================================
    c2_bf, mc_bf, mp_bf, al_bf, sv_bf = free_results[0]
    lam_bf = al_bf * mc_bf / mp_bf

    print(f"\n  {'='*60}")
    print(f"  FREE BEST FIT (chi2/dof = {c2_bf:.2f}/{ndof} = {c2_bf/ndof:.2f})")
    print(f"  {'='*60}")
    print(f"    m_chi = {mc_bf:.4f} GeV")
    print(f"    m_phi = {mp_bf*1e3:.4f} MeV")
    print(f"    alpha = {al_bf:.6e}")
    print(f"    lambda= {lam_bf:.4f}")

    print(f"\n  {'System':<25s} {'v':>6s} {'Theory':>8s} {'Obs':>8s} {'Pull':>7s}")
    print(f"  {'-'*60}")
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        theory = sv_bf.get(v, 0.0)
        if theory >= central:
            sig = hi - central if hi > central else 0.5 * central
        else:
            sig = central - lo if central > lo else 0.5 * central
        if sig <= 0:
            sig = 0.5 * max(central, 0.01)
        pull = (theory - central) / sig
        print(f"  {name:<25s} {v:>6.0f} {theory:>8.4f} {central:>8.2f} {pull:>+7.2f}")

    print(f"\n  TOP 10 FREE-FIT:")
    print(f"  {'#':>3s} {'chi2':>8s} {'m_chi':>8s} {'m_phi':>8s} {'alpha':>12s} {'lam':>6s} {'s/m(30)':>8s} {'s/m(1k)':>8s}")
    print(f"  {'-'*70}")
    for i, (c2, mc, mp, al, sv) in enumerate(free_results[:10]):
        lam = al * mc / mp
        sm1k = sv.get(1000, sv.get(1100, sv.get(1200, 0)))
        print(f"  {i+1:>3d} {c2:>8.2f} {mc:>8.2f} {mp*1e3:>8.2f} {al:>12.4e} {lam:>6.2f} {sv.get(30,0):>8.4f} {sm1k:>8.4f}")

    # ============================================================
    # RELIC-CONSTRAINED BEST FIT (17 BPs)
    # ============================================================
    if relic_results:
        c2_r, mc_r, mp_r, al_r, om_r, sv_r = relic_results[0]
        lam_r = al_r * mc_r / mp_r

        print(f"\n  {'='*60}")
        print(f"  RELIC-CONSTRAINED BEST FIT (Omega_h2 = {om_r:.4f})")
        print(f"  (chi2/dof = {c2_r:.2f}/{ndof} = {c2_r/ndof:.2f})")
        print(f"  {'='*60}")
        print(f"    m_chi = {mc_r:.4f} GeV")
        print(f"    m_phi = {mp_r*1e3:.4f} MeV")
        print(f"    alpha = {al_r:.6e}")
        print(f"    lambda= {lam_r:.4f}")

        print(f"\n  {'System':<25s} {'v':>6s} {'Theory':>8s} {'Obs':>8s} {'Pull':>7s}")
        print(f"  {'-'*60}")
        for name, v, central, lo, hi, ref in OBSERVATIONS:
            theory = sv_r.get(v, 0.0)
            if theory >= central:
                sig = hi - central if hi > central else 0.5 * central
            else:
                sig = central - lo if central > lo else 0.5 * central
            if sig <= 0:
                sig = 0.5 * max(central, 0.01)
            pull = (theory - central) / sig
            print(f"  {name:<25s} {v:>6.0f} {theory:>8.4f} {central:>8.2f} {pull:>+7.2f}")

        print(f"\n  ALL 17 RELIC BPs (ranked by chi2):")
        print(f"  {'#':>3s} {'chi2':>8s} {'c2/dof':>8s} {'m_chi':>8s} {'m_phi':>8s} {'alpha':>12s} {'Omh2':>8s} {'s/m(30)':>8s}")
        print(f"  {'-'*75}")
        for i, (c2, mc, mp, al, om, sv) in enumerate(relic_results):
            print(f"  {i+1:>3d} {c2:>8.2f} {c2/ndof:>8.2f} {mc:>8.2f} {mp*1e3:>8.2f} {al:>12.4e} {om:>8.4f} {sv.get(30,0):>8.4f}")

    # ============================================================
    # SUMMARY
    # ============================================================
    bp1_chi2 = None
    for c2, mc, mp, al, om, sv in relic_results:
        if abs(mc - 20.69) < 0.1:
            bp1_chi2 = c2
            break

    print(f"\n  {'='*60}")
    print(f"  SUMMARY")
    print(f"  {'='*60}")
    print(f"    Free best-fit:    chi2/dof = {c2_bf/ndof:.2f}")
    if relic_results:
        print(f"    Relic best-fit:   chi2/dof = {c2_r/ndof:.2f}")
    if bp1_chi2 is not None:
        print(f"    BP1:              chi2/dof = {bp1_chi2/ndof:.2f}")
    print(f"    Perfect fit:      chi2/dof = 1.00")
    print(f"    Total time:       {time.time()-t0:.1f}s")

    # ============================================================
    # PLOT
    # ============================================================
    print(f"\n  Generating plot...")
    from v22_raw_scan import sigma_T_vpm
    sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)  # warm up in main

    v_plot = np.logspace(np.log10(8), np.log10(6000), 200)
    curve_bf = np.array([sigma_T_vpm(mc_bf, mp_bf, al_bf, v) for v in v_plot])
    curve_bp1 = np.array([sigma_T_vpm(20.69, 11.34e-3, 1.048e-3, v) for v in v_plot])

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    ax.plot(v_plot, curve_bf, 'b-', lw=2.5,
            label=f'Best fit (free): $m_\\chi$={mc_bf:.1f}, $m_\\phi$={mp_bf*1e3:.1f} MeV, '
                  f'$\\alpha$={al_bf:.2e}\n    $\\chi^2$/dof={c2_bf/ndof:.2f}')

    if relic_results:
        curve_r = np.array([sigma_T_vpm(mc_r, mp_r, al_r, v) for v in v_plot])
        ax.plot(v_plot, curve_r, 'g-.', lw=2.0,
                label=f'Best fit (relic): $m_\\chi$={mc_r:.1f}, $m_\\phi$={mp_r*1e3:.1f} MeV\n'
                      f'    $\\Omega h^2$={om_r:.3f}, $\\chi^2$/dof={c2_r/ndof:.2f}')

    bp1_label = f'BP1: $m_\\chi$=20.7, $m_\\phi$=11.3'
    if bp1_chi2 is not None:
        bp1_label += f'\n    $\\chi^2$/dof={bp1_chi2/ndof:.2f}'
    ax.plot(v_plot, curve_bp1, 'r--', lw=1.8, alpha=0.7, label=bp1_label)

    colors = {'KTY16': 'darkgreen', 'KKPY17': 'purple',
              'Randall+08': 'orange', 'Harvey+15': 'brown', 'Elbert+15': 'red'}
    markers = {'KTY16': 'o', 'KKPY17': 's', 'Randall+08': 'D',
               'Harvey+15': '^', 'Elbert+15': 'v'}
    plotted = set()
    for name, v, central, lo, hi, ref in OBSERVATIONS:
        label = ref if ref not in plotted else None
        plotted.add(ref)
        ax.errorbar(v, central, yerr=[[central - lo], [hi - central]],
                     fmt=markers.get(ref, 'o'), color=colors.get(ref, 'gray'),
                     ms=8, capsize=4, capthick=1.5, elinewidth=1.5, label=label)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$v$ [km/s]', fontsize=14)
    ax.set_ylabel('$\\sigma/m$ [cm$^2$/g]', fontsize=14)
    ax.set_title('$\\chi^2$ Fit to Astrophysical Data (80k scan points)', fontsize=15)
    ax.set_xlim(8, 6000)
    ax.set_ylim(1e-3, 20)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')
    ax.axvspan(10, 60, color='cyan', alpha=0.06)
    ax.axvspan(800, 2000, color='pink', alpha=0.06)
    ax.text(25, 15, 'Dwarfs', fontsize=9, ha='center', color='teal')
    ax.text(1200, 15, 'Clusters', fontsize=9, ha='center', color='crimson')

    plt.tight_layout()
    outpng = os.path.join(_DIR, 'v34_chi2_fit.png')
    fig.savefig(outpng, dpi=150)
    plt.close()
    print(f"  Saved figure -> {outpng}")
    print(f"\n  Total time: {time.time()-t0:.1f}s")
    print(hdr)


if __name__ == '__main__':
    mp_lib.freeze_support()
    run()
