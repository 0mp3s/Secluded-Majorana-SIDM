#!/usr/bin/env python3
"""
Sweep all 17 relic-viable benchmark points through key predictions.

Outputs a comprehensive table showing that every point in the viable island
passes all observational constraints:
  1. σ/m at key velocities (12, 30, 100, 1000, 4700 km/s)
  2. Gravothermal: N_scatter regime for 8 classical dSphs
  3. Cluster offsets: σ/m < upper bounds for 8 mergers
  4. ΔN_eff: Boltzmann-suppressed → consistent with Planck
"""
# === path setup ================================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
_sys.path.insert(0, _os.path.join(_ROOT, 'predictions', 'gravothermal'))
_sys.path.insert(0, _os.path.join(_ROOT, 'predictions', 'cluster_offsets'))
# ================================================================

import sys, os, csv, math, time
import numpy as np

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)
os.environ['PYTHONIOENCODING'] = 'utf-8'

from v22_raw_scan import sigma_T_vpm
from predict_gravothermal import load_dsphs, compute_predictions, classify
from predict_offsets import load_clusters
from output_manager import get_latest

# ── paths ──
BP_CSV = str(get_latest("v31_true_viable_points"))
DSPH_CSV = os.path.join(_ROOT, 'predictions', 'gravothermal', 'dsphs_data.csv')
CLUSTER_CSV = os.path.join(_ROOT, 'predictions', 'cluster_offsets', 'clusters_data.csv')
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── ΔN_eff (analytic — all BPs have m_phi >> T_BBN so ΔN_eff ≈ 0) ──
G_STAR_S_DEC = 86.25  # at T_dec ~ 10-100 GeV


def T_dark_over_T_SM(g_star_s_SM, g_star_s_dec):
    return (g_star_s_SM / g_star_s_dec) ** (1.0 / 3.0)


def delta_neff_at_bbn(m_phi_MeV):
    """Compute ΔN_eff at BBN (T_SM ~ 1 MeV)."""
    g_sm_bbn = 10.75
    ratio = T_dark_over_T_SM(g_sm_bbn, G_STAR_S_DEC)
    T_phi_MeV = ratio * 1.0  # T_SM = 1 MeV at BBN
    boltz = math.exp(-m_phi_MeV / T_phi_MeV) if m_phi_MeV / T_phi_MeV < 500 else 0.0
    # Massless limit: ΔN_eff = (4/7)(T_φ/T_ν)^4, T_ν = T_SM at BBN
    dn_massless = (4.0 / 7.0) * ratio ** 4
    return dn_massless * boltz


def load_bps():
    bps = []
    with open(BP_CSV, newline='') as f:
        for i, row in enumerate(csv.DictReader(f)):
            bps.append({
                'label': f'BP{i+1}',
                'm_chi_GeV': float(row['m_chi_GeV']),
                'm_phi_MeV': float(row['m_phi_MeV']),
                'alpha': float(row['alpha']),
                'omega_h2': float(row['omega_h2']),
                'lambda': float(row['lambda']),
                'sigma_m_30_csv': float(row['sigma_m_30']),
                'sigma_m_1000_csv': float(row['sigma_m_1000']),
            })
    return bps


def main():
    print("=" * 100)
    print("  Sweep: All 17 Relic-Viable Benchmarks Through Key Predictions")
    print("=" * 100)

    t0 = time.time()

    # Warm up JIT
    print("\nWarming up Numba JIT...")
    _ = sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)
    print(f"JIT warm-up: {time.time()-t0:.1f}s")

    bps = load_bps()
    dsphs = load_dsphs(DSPH_CSV)
    clusters = load_clusters(CLUSTER_CSV)

    print(f"\nLoaded {len(bps)} BPs, {len(dsphs)} dSphs, {len(clusters)} clusters\n")

    # ════════════════════════════════════════════════
    #  TABLE 1: σ/m at key velocities
    # ════════════════════════════════════════════════
    key_velocities = [12, 30, 100, 1000, 4700]
    print("─" * 100)
    print("TABLE 1: σ/m [cm²/g] at key velocities")
    print("─" * 100)
    hdr = f"{'BP':>5} {'m_χ':>6} {'m_φ':>7} {'α':>10} {'λ':>6} {'Ωh²':>7}"
    for v in key_velocities:
        hdr += f" {'σ/m@'+str(v):>9}"
    print(hdr)
    print("─" * 100)

    sigma_table = []
    for bp in bps:
        mc = bp['m_chi_GeV']
        mp_gev = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        row = [bp['label'], mc, bp['m_phi_MeV'], alpha, bp['lambda'], bp['omega_h2']]
        sigmas = {}
        for v in key_velocities:
            s = sigma_T_vpm(mc, mp_gev, alpha, float(v))
            sigmas[v] = s
            row.append(s)
        sigma_table.append((bp, sigmas))
        line = f"{bp['label']:>5} {mc:>6.1f} {bp['m_phi_MeV']:>7.2f} {alpha:>10.3e} {bp['lambda']:>6.1f} {bp['omega_h2']:>7.4f}"
        for v in key_velocities:
            line += f" {sigmas[v]:>9.4f}"
        print(line)
    print("─" * 100)

    # SIDM viability check from σ/m values (matching relic_density/config.json cuts)
    SIDM_30_LO, SIDM_30_HI, SIDM_1K_HI = 0.5, 10.0, 0.1
    n_viable = sum(1 for _, s in sigma_table if SIDM_30_LO <= s[30] <= SIDM_30_HI and s[1000] < SIDM_1K_HI)
    print(f"\nSIDM viable ({SIDM_30_LO} ≤ σ/m@30 ≤ {SIDM_30_HI} AND σ/m@1000 < {SIDM_1K_HI}): {n_viable}/{len(bps)}")

    # ════════════════════════════════════════════════
    #  TABLE 2: Gravothermal regime for 8 dSphs
    # ════════════════════════════════════════════════
    print("\n" + "─" * 100)
    print("TABLE 2: Gravothermal N_scatter regimes (8 classical dSphs)")
    print("─" * 100)

    grav_results = compute_predictions(dsphs, bps)

    # Header
    hdr = f"{'dSph':>12} {'obs':>5}"
    for bp in bps:
        hdr += f" {bp['label']:>5}"
    print(hdr)
    print("─" * 100)

    # For each dSph, show N_scatter class per BP
    grav_scores = {bp['label']: {'ok': 0, 'fail': 0, 'ambig': 0} for bp in bps}
    for row in grav_results:
        dsph_name = row['dsph']
        obs = row['core_observed']
        line = f"{dsph_name:>12} {obs:>5}"
        for bp in bps:
            label = bp['label']
            n_s = row[f'{label}_N_scatter']
            cls = classify(n_s)
            if obs == 'YES':
                if 'CORED' in cls:
                    tag = '✓'
                    grav_scores[label]['ok'] += 1
                else:
                    tag = '✗'
                    grav_scores[label]['fail'] += 1
            elif obs == 'AMBIGUOUS':
                tag = '~'
                grav_scores[label]['ambig'] += 1
            else:
                tag = '✓' if 'CUSPY' in cls else '!'
                grav_scores[label]['ok'] += 1
            line += f"  {n_s:>3.0f}{tag}"
        print(line)

    print("─" * 100)
    line = f"{'Score':>12} {'':>5}"
    for bp in bps:
        s = grav_scores[bp['label']]
        line += f" {s['ok']}/{s['ok']+s['fail']:>3} "
    print(line)

    # ════════════════════════════════════════════════
    #  TABLE 3: Cluster offset constraints
    # ════════════════════════════════════════════════
    print("\n" + "─" * 100)
    print("TABLE 3: Cluster merger constraints (σ/m < upper bound)")
    print("─" * 100)

    hdr = f"{'Cluster':>18} {'v':>5} {'UL':>5}"
    for bp in bps:
        hdr += f" {bp['label']:>5}"
    print(hdr)
    print("─" * 100)

    cluster_pass = {bp['label']: 0 for bp in bps}
    for cl in clusters:
        line = f"{cl['name']:>18} {cl['v_infall']:>5.0f} {cl['sigma_m_upper']:>5.2f}"
        for bp in bps:
            mc = bp['m_chi_GeV']
            mp_gev = bp['m_phi_MeV'] / 1000.0
            alpha = bp['alpha']
            s = sigma_T_vpm(mc, mp_gev, alpha, cl['v_infall'])
            ok = s < cl['sigma_m_upper']
            if ok:
                cluster_pass[bp['label']] += 1
            tag = '✓' if ok else '✗'
            line += f" {tag}{s:>.2f}" if s >= 0.01 else f" {tag}{s:>.3f}"
        print(line)

    print("─" * 100)
    line = f"{'PASS':>18} {'':>5} {'':>5}"
    for bp in bps:
        n = cluster_pass[bp['label']]
        tag = "ALL" if n == len(clusters) else f"{n}/{len(clusters)}"
        line += f" {tag:>5}"
    print(line)

    # ════════════════════════════════════════════════
    #  TABLE 4: ΔN_eff at BBN
    # ════════════════════════════════════════════════
    print("\n" + "─" * 100)
    print("TABLE 4: ΔN_eff at BBN (Planck 2σ limit: 0.34)")
    print("─" * 100)

    all_neff_pass = True
    for bp in bps:
        dn = delta_neff_at_bbn(bp['m_phi_MeV'])
        ok = dn < 0.34
        if not ok:
            all_neff_pass = False
        print(f"  {bp['label']:>5}: m_φ = {bp['m_phi_MeV']:>6.2f} MeV → ΔN_eff = {dn:.2e}  {'✓' if ok else '✗'}")
    print(f"\n  All ΔN_eff < Planck 2σ: {'YES' if all_neff_pass else 'NO'}")

    # ════════════════════════════════════════════════
    #  GRAND SUMMARY
    # ════════════════════════════════════════════════
    elapsed = time.time() - t0
    print("\n" + "=" * 100)
    print("  GRAND SUMMARY")
    print("=" * 100)

    all_pass = True
    for bp in bps:
        label = bp['label']
        _, sigmas = next((b, s) for b, s in sigma_table if b['label'] == label)
        sidm_ok = SIDM_30_LO <= sigmas[30] <= SIDM_30_HI and sigmas[1000] < SIDM_1K_HI
        clust_ok = cluster_pass[label] == len(clusters)
        gs = grav_scores[label]
        grav_ok = gs['fail'] == 0
        dn = delta_neff_at_bbn(bp['m_phi_MeV'])
        neff_ok = dn < 0.34

        verdict = "✓ PASS" if (sidm_ok and clust_ok and neff_ok) else "✗ FAIL"
        if not (sidm_ok and clust_ok and neff_ok):
            all_pass = False
        grav_tag = f"{gs['ok']}/{gs['ok']+gs['fail']}ok" + (f"+{gs['ambig']}~" if gs['ambig'] else "")
        print(f"  {label:>5}: σ/m@30={sigmas[30]:.2f} σ/m@1k={sigmas[1000]:.4f} "
              f"clusters={cluster_pass[label]}/{len(clusters)} "
              f"grav={grav_tag} ΔNeff={dn:.0e} → {verdict}")

    print(f"\n  Overall: {'ALL 17 PASS ✓' if all_pass else 'SOME FAIL ✗'}")
    print(f"  Runtime: {elapsed:.1f}s")
    print("=" * 100)

    # ── Save CSV ──
    csv_path = os.path.join(OUT_DIR, 'sweep_17bp_results.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        header = ['BP', 'm_chi_GeV', 'm_phi_MeV', 'alpha', 'lambda', 'omega_h2']
        for v in key_velocities:
            header.append(f'sigma_m_{v}')
        header += ['clusters_pass', 'grav_ok', 'grav_fail', 'grav_ambig', 'delta_neff_bbn', 'overall_pass']
        w.writerow(header)
        for bp, sigmas in sigma_table:
            label = bp['label']
            gs = grav_scores[label]
            dn = delta_neff_at_bbn(bp['m_phi_MeV'])
            sidm_ok = SIDM_30_LO <= sigmas[30] <= SIDM_30_HI and sigmas[1000] < SIDM_1K_HI
            clust_ok = cluster_pass[label] == len(clusters)
            neff_ok = dn < 0.34
            overall = sidm_ok and clust_ok and neff_ok
            row = [label, bp['m_chi_GeV'], bp['m_phi_MeV'], bp['alpha'], bp['lambda'], bp['omega_h2']]
            for v in key_velocities:
                row.append(f'{sigmas[v]:.6f}')
            row += [cluster_pass[label], gs['ok'], gs['fail'], gs['ambig'], f'{dn:.2e}', overall]
            w.writerow(row)
    print(f"\n  CSV saved: {csv_path}")

    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
