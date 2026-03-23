#!/usr/bin/env python3
"""
predictions/multi_dsph_jeans/predict_multi_dsph_jeans.py
========================================================
Multi-dSph validation: Anisotropic Jeans (SIDM + Feedback + OM β)

Tests the Secluded Majorana model on 5 classical MW dSphs to verify
that the Fornax results are not a statistical fluke:

    Fornax · Sculptor · Draco · Carina · Sextans

For each dSph:
  - NFW (cuspy CDM)
  - coreNFW  (feedback only, Read+2016)
  - BP1 / MAP  SIDM only  (on NFW)
  - BP1 / MAP  SIDM + feedback  (on coreNFW)
  - Osipkov-Merritt anisotropy scan  r_a ∈ [0.3 … 50] kpc

Kinematic data: Walker+2007/2009 binned σ_los(R) profiles.
Halo parameters: Read+2019, Wolf+2010, McConnachie 2012.

References:
  Walker+2007, ApJ 667, L53
  Walker+2009, AJ 137, 3100
  Read+2019, MNRAS 484, 1401
  McConnachie 2012, AJ 144, 4
  Wolf+2010, MNRAS 406, 1220
  Mamon & Łokas 2005, MNRAS 363, 705
"""
import sys, os, math
import numpy as np

_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v22_raw_scan import sigma_T_vpm
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)      # JIT warmup

# ── Physical constants ──
G_KPC_MSUN = 4.302e-6       # kpc (km/s)² / M_sun
MSUN_G     = 1.989e33       # g
KPC_CM     = 3.086e21       # cm
GYR_S      = 3.156e16       # s
KM_S_CM_S  = 1e5            # cm/s

# ── Benchmark points ──
BPS = [
    {"label": "BP1", "m_chi_GeV": 20.69, "m_phi_MeV": 11.34, "alpha": 1.048e-3},
    {"label": "MAP", "m_chi_GeV": 90.64, "m_phi_MeV": 13.85, "alpha": 2.546e-2},
]

# ═══════════════════════════════════════════════════════════════════
#  dSph data — kinematic profiles from Walker+2007/2009
#  Halo params: NFW fits consistent with Wolf+2010 mass estimator
# ═══════════════════════════════════════════════════════════════════
DSPHS = {
    "Fornax": {
        "M200": 3.16e9, "c200": 18,
        "r_half_pc": 710, "M_star": 2.0e7,
        "sigma_v": 11.7, "t_age": 10,
        "R_pc":      [100,  250,  400,  550,  710,  900,  1100, 1400, 1800],
        "sigma_los": [11.4, 11.8, 12.0, 11.7, 11.5, 11.1, 10.5,  9.6,  8.1],
        "sigma_err": [ 1.0,  0.7,  0.6,  0.6,  0.6,  0.8,  1.0,  1.3,  2.0],
        "n_stars": 2633, "ref": "Walker+2009",
    },
    "Sculptor": {
        "M200": 1.5e9, "c200": 20,
        "r_half_pc": 283, "M_star": 2.3e6,
        "sigma_v": 9.2, "t_age": 10,
        "R_pc":      [ 30,   83,  133,  183,  240,  300,  380,  480,  620,  830],
        "sigma_los": [10.1, 10.0,  9.7,  9.3,  9.0,  8.6,  8.1,  7.3,  6.4,  5.3],
        "sigma_err": [ 1.5,  0.8,  0.6,  0.5,  0.5,  0.4,  0.5,  0.6,  0.8,  1.4],
        "n_stars": 1365, "ref": "Walker+2009",
    },
    "Draco": {
        "M200": 8e8, "c200": 22,
        "r_half_pc": 221, "M_star": 2.9e5,
        "sigma_v": 9.1, "t_age": 10,
        "R_pc":      [ 25,   70,  110,  155,  205,  265,  340,  450],
        "sigma_los": [10.4,  9.8,  9.5,  9.2,  8.8,  8.5,  7.8,  6.8],
        "sigma_err": [ 2.0,  1.0,  0.8,  0.7,  0.8,  0.8,  1.0,  1.6],
        "n_stars": 292, "ref": "Walker+2009",
    },
    "Carina": {
        "M200": 4e8, "c200": 26,
        "r_half_pc": 250, "M_star": 3.8e5,
        "sigma_v": 6.6, "t_age": 7,
        "R_pc":      [ 30,   85,  140,  200,  265,  340,  440,  580],
        "sigma_los": [ 6.8,  6.7,  6.6,  6.4,  6.3,  6.0,  5.6,  4.9],
        "sigma_err": [ 1.2,  0.7,  0.5,  0.5,  0.5,  0.6,  0.7,  1.1],
        "n_stars": 774, "ref": "Walker+2009",
    },
    "Sextans": {
        "M200": 1.0e9, "c200": 20,
        "r_half_pc": 695, "M_star": 4.4e5,
        "sigma_v": 7.9, "t_age": 10,
        "R_pc":      [  50,  150,  260,  380,  520,  680,  880, 1150],
        "sigma_los": [ 8.2,  8.0,  7.8,  7.6,  7.3,  6.9,  6.4,  5.5],
        "sigma_err": [ 1.5,  0.9,  0.7,  0.6,  0.6,  0.7,  0.9,  1.5],
        "n_stars": 424, "ref": "Walker+2009",
    },
}

# ═══════════════════════════════════════════════════════════════════
#  NFW
# ═══════════════════════════════════════════════════════════════════
def nfw_params(M200, c200):
    rho_crit = 277.5 * 0.674**2
    R200 = (3*M200 / (4*math.pi*200*rho_crit))**(1/3)
    r_s = R200 / c200
    gc = math.log(1+c200) - c200/(1+c200)
    rho_s = M200 / (4*math.pi*r_s**3*gc)
    return rho_s, r_s, R200

def rho_nfw(r, rho_s, r_s):
    x = r/r_s
    return rho_s / (x*(1+x)**2)

def M_nfw(r, rho_s, r_s):
    x = r/r_s
    return 4*math.pi*rho_s*r_s**3*(math.log(1+x) - x/(1+x))


# ═══════════════════════════════════════════════════════════════════
#  coreNFW (Read+2016)
# ═══════════════════════════════════════════════════════════════════
def core_nfw_n(M_star, M_halo, kappa=80.0):
    return math.tanh(kappa * M_star / M_halo)

def M_core_nfw(r, rho_s, r_s, r_c, n):
    return M_nfw(r, rho_s, r_s) * math.tanh(r/r_c)**n

def rho_core_nfw(r, rho_s, r_s, r_c, n, dr_frac=1e-4):
    dr = max(r*dr_frac, 1e-6)
    Mp = M_core_nfw(r+dr, rho_s, r_s, r_c, n)
    Mm = M_core_nfw(r-dr if r>dr else 0.0, rho_s, r_s, r_c, n)
    dMdr = (Mp-Mm)/(2*dr if r>dr else dr+r)
    return max(dMdr/(4*math.pi*r**2), 0.0)


# ═══════════════════════════════════════════════════════════════════
#  SIDM coring
# ═══════════════════════════════════════════════════════════════════
def sidm_matching(r_arr, rho_base, sigma_over_m, sigma_v, t_age_Gyr):
    rho_conv = MSUN_G / KPC_CM**3
    v_rel = sigma_v * math.sqrt(2) * KM_S_CM_S
    t_age = t_age_Gyr * GYR_S
    for i in range(len(r_arr)-1, 0, -1):
        rho_cgs = rho_base[i]*rho_conv
        n_scat = sigma_over_m * rho_cgs * v_rel * t_age
        if n_scat >= 1.0:
            rho_prev = rho_base[i+1]*rho_conv if i+1<len(r_arr) else 0
            n_prev = sigma_over_m * rho_prev * v_rel * t_age
            if n_prev < 1.0 and n_scat > 1.0:
                frac = (1-n_prev)/(n_scat-n_prev) if (n_scat-n_prev)>0 else 0.5
                r1 = r_arr[i+1]+frac*(r_arr[i]-r_arr[i+1]) if i+1<len(r_arr) else r_arr[i]
            else:
                r1 = r_arr[i]
            rho0 = np.interp(r1, r_arr, rho_base)
            return r1, rho0
    return 0.0, rho_base[0]

def build_combined(r_arr, rho_base, M_base, sigma_over_m, sigma_v, t_age):
    if sigma_over_m < 1e-10:
        return rho_base.copy(), M_base.copy(), 0.0
    r1, rho0 = sidm_matching(r_arr, rho_base, sigma_over_m, sigma_v, t_age)
    if r1 < r_arr[1]:
        return rho_base.copy(), M_base.copy(), 0.0
    rho_c = np.copy(rho_base)
    for i, r in enumerate(r_arr):
        if r <= r1:
            rho_c[i] = rho0 / (1+(r/r1)**2)
    M_c = np.zeros_like(r_arr)
    for i in range(1, len(r_arr)):
        dr = r_arr[i]-r_arr[i-1]
        rm = 0.5*(r_arr[i]+r_arr[i-1])
        M_c[i] = M_c[i-1] + 4*math.pi*rm**2*0.5*(rho_c[i]+rho_c[i-1])*dr
    return rho_c, M_c, r1


# ═══════════════════════════════════════════════════════════════════
#  Plummer stellar profile
# ═══════════════════════════════════════════════════════════════════
def rho_plummer(r, M_star, a):
    return (3*M_star)/(4*math.pi*a**3) * (1+(r/a)**2)**(-2.5)

def M_plummer(r, M_star, a):
    return M_star * r**3 / (r**2+a**2)**1.5


# ═══════════════════════════════════════════════════════════════════
#  Anisotropic Jeans: Osipkov-Merritt
# ═══════════════════════════════════════════════════════════════════
def solve_jeans_OM(r_arr, rho_star, M_tot, r_a):
    Q = 1.0 + (r_arr/r_a)**2
    integ = rho_star * Q * G_KPC_MSUN * M_tot / r_arr**2
    cum = np.zeros(len(r_arr))
    for i in range(len(r_arr)-2, -1, -1):
        dr = r_arr[i+1]-r_arr[i]
        cum[i] = cum[i+1] + 0.5*(integ[i]+integ[i+1])*dr
    sr2 = np.zeros_like(r_arr)
    mask = (rho_star > 0) & (Q > 0)
    sr2[mask] = cum[mask]/(rho_star[mask]*Q[mask])
    return sr2

def project_sigma_los_OM(R_proj, r_arr, rho_star, sr2, r_a):
    slv = np.zeros(len(R_proj))
    for j, R in enumerate(R_proj):
        mask = r_arr > R*1.001
        rs = r_arr[mask]; rhos = rho_star[mask]; s2s = sr2[mask]
        if len(rs) < 3:
            continue
        factor = 1.0 - R**2/(rs**2 + r_a**2)
        K = factor * rs / np.sqrt(rs**2 - R**2)
        I_num = np.trapezoid(K * rhos * s2s, rs)
        Sig = 2.0*np.trapezoid(rhos*rs/np.sqrt(rs**2-R**2), rs)
        if Sig > 0:
            slv[j] = math.sqrt(max(2*I_num/Sig, 0.0))
    return slv

def chi2_calc(smod, R_proj, R_obs, s_obs, s_err):
    sm = np.interp(R_obs, R_proj, smod)
    return np.sum(((s_obs-sm)/s_err)**2), sm


# ═══════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════
def main():
    out_dir = os.path.join(_DIR, 'output')
    os.makedirs(out_dir, exist_ok=True)

    r_a_scan = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                          1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0])

    print("=" * 100)
    print("  Multi-dSph Anisotropic Jeans Validation — SIDM + Feedback + Osipkov-Merritt β(r)")
    print("=" * 100)
    print()

    # ── Grand summary storage ──
    grand = {}     # dSph → {scenario_tag → result_dict}

    dsph_names = list(DSPHS.keys())

    for dname in dsph_names:
        d = DSPHS[dname]
        M200   = d['M200']
        c200   = d['c200']
        r_half = d['r_half_pc'] / 1000.0     # kpc
        M_star = d['M_star']
        sv     = d['sigma_v']
        t_age  = d['t_age']

        R_obs = np.array(d['R_pc']) / 1000.0
        s_obs = np.array(d['sigma_los'])
        s_err = np.array(d['sigma_err'])
        ndata = len(s_obs)

        rho_s, r_s, R200 = nfw_params(M200, c200)
        n_fb = core_nfw_n(M_star, M200)

        # Radial grid
        r_max = max(10.0, R_obs[-1]*3)   # extend well past data
        r_arr = np.logspace(np.log10(0.001), np.log10(r_max), 3000)
        R_proj = np.linspace(0.01, R_obs[-1]*1.5, 120)

        # Stellar arrays
        rho_star = np.array([rho_plummer(r, M_star, r_half) for r in r_arr])
        M_star_a = np.array([M_plummer(r, M_star, r_half)   for r in r_arr])

        # NFW arrays
        rho_nfw_a = np.array([rho_nfw(r, rho_s, r_s) for r in r_arr])
        M_nfw_a   = np.array([M_nfw(r, rho_s, r_s)   for r in r_arr])

        # coreNFW (feedback)
        r_c = 1.75 * r_half
        rho_cnfw = np.array([rho_core_nfw(r, rho_s, r_s, r_c, n_fb) for r in r_arr])
        M_cnfw   = np.array([M_core_nfw(r, rho_s, r_s, r_c, n_fb)   for r in r_arr])

        print(f"  ┌── {dname} ({'★' if dname == 'Fornax' else ''}) ──")
        print(f"  │  M200={M200:.2e}  c={c200}  r_s={r_s:.3f} kpc")
        print(f"  │  r_half={r_half*1000:.0f} pc  M_★={M_star:.1e}  "
              f"σ_v={sv} km/s  n_fb={n_fb:.3f}  N_★={d['n_stars']}")

        # ── Build scenarios ──
        scenarios = [
            {'tag': 'NFW', 'M_dm': M_nfw_a},
            {'tag': 'FB',  'M_dm': M_cnfw},
        ]
        for bp in BPS:
            lab = bp['label']
            m_chi = bp['m_chi_GeV']
            m_phi = bp['m_phi_MeV'] / 1000.0
            alpha = bp['alpha']
            v_rel = sv * math.sqrt(2)
            sm = sigma_T_vpm(m_chi, m_phi, alpha, v_rel)

            # SIDM on NFW
            _, M_s, r1s = build_combined(
                r_arr, rho_nfw_a, M_nfw_a, sm, sv, t_age)
            scenarios.append({'tag': f'{lab}_SIDM', 'M_dm': M_s,
                              'sigma_m': sm, 'r1': r1s})

            # SIDM + feedback on coreNFW
            _, M_b, r1b = build_combined(
                r_arr, rho_cnfw, M_cnfw, sm, sv, t_age)
            scenarios.append({'tag': f'{lab}_SIDM_FB', 'M_dm': M_b,
                              'sigma_m': sm, 'r1': r1b})

        # ── Scan r_a for each scenario ──
        results = {}
        for sc in scenarios:
            M_tot = sc['M_dm'] + M_star_a

            # Isotropic (r_a → ∞)
            sr2_iso = solve_jeans_OM(r_arr, rho_star, M_tot, 50.0)
            slos_iso = project_sigma_los_OM(R_proj, r_arr, rho_star, sr2_iso, 50.0)
            c2_iso, _ = chi2_calc(slos_iso, R_proj, R_obs, s_obs, s_err)

            best_c2 = c2_iso
            best_ra = 50.0
            best_slos = slos_iso
            chi2_ra = []

            for ra in r_a_scan:
                sr2 = solve_jeans_OM(r_arr, rho_star, M_tot, ra)
                sl = project_sigma_los_OM(R_proj, r_arr, rho_star, sr2, ra)
                c2, _ = chi2_calc(sl, R_proj, R_obs, s_obs, s_err)
                chi2_ra.append(c2)
                if c2 < best_c2:
                    best_c2 = c2; best_ra = ra; best_slos = sl.copy()

            ndof = ndata - 1  # 1 free param: r_a
            _, best_obs = chi2_calc(best_slos, R_proj, R_obs, s_obs, s_err)
            results[sc['tag']] = {
                'chi2_iso': c2_iso,
                'best_ra': best_ra,
                'best_chi2': best_c2,
                'best_slos': best_slos,
                'best_obs': best_obs,
                'chi2_ra': chi2_ra,
                'ndof': ndof,
                'R_proj': R_proj,
            }

        grand[dname] = results

        # ── Per-dSph mini table ──
        print(f"  │  {'Scenario':18s} │ χ²(β=0) │ best r_a │ χ²_best │ χ²/dof")
        print(f"  │  " + "─"*18 + "─┼─────────┼──────────┼─────────┼───────")
        for sc in scenarios:
            r = results[sc['tag']]
            print(f"  │  {sc['tag']:18s} │ {r['chi2_iso']:7.1f} │ "
                  f"{r['best_ra']:6.2f}   │ {r['best_chi2']:7.1f} │ "
                  f"{r['best_chi2']/r['ndof']:5.2f}")
        best_sc = min(results.items(), key=lambda x: x[1]['best_chi2'])
        print(f"  │  ★ Best: {best_sc[0]}, "
              f"r_a = {best_sc[1]['best_ra']:.2f}, "
              f"χ²/dof = {best_sc[1]['best_chi2']:.1f}/{best_sc[1]['ndof']}")
        print(f"  └──")
        print()

    # ═══════════════════════════════════════════════════════════════
    #  Grand summary table
    # ═══════════════════════════════════════════════════════════════
    KEY_TAGS = ['NFW', 'MAP_SIDM_FB', 'BP1_SIDM_FB']
    print("=" * 100)
    print("  GRAND SUMMARY — χ² (isotropic → best β) for key scenarios")
    print("=" * 100)
    hdr = f"  {'dSph':12s}│ {'N_data':>6s}"
    for tag in KEY_TAGS:
        hdr += f" │ {tag+' iso':>14s} → {'β-opt':>7s} (r_a)"
    print(hdr)
    print("  " + "─"*12 + "┼" + "─"*7 + ("─┼" + "─"*35)*len(KEY_TAGS))

    for dname in dsph_names:
        ndata = len(DSPHS[dname]['R_pc'])
        line = f"  {dname:12s}│ {ndata:>6d}"
        for tag in KEY_TAGS:
            r = grand[dname][tag]
            line += (f" │ {r['chi2_iso']:7.1f}    → {r['best_chi2']:7.1f}"
                     f" ({r['best_ra']:4.1f})")
        print(line)
    print()

    # ── Across-dSph MAP summary ──
    print("  ── MAP SIDM+fb + best β: summary ──")
    all_chi2_dof = []
    for dname in dsph_names:
        r = grand[dname]['MAP_SIDM_FB']
        c2d = r['best_chi2']/r['ndof']
        all_chi2_dof.append(c2d)
        verdict = "★★" if c2d < 1.5 else ("★" if c2d < 2.5 else "")
        print(f"    {dname:12s}: χ²/dof = {r['best_chi2']:.1f}/{r['ndof']} = "
              f"{c2d:.2f}  r_a = {r['best_ra']:.1f} kpc  {verdict}")

    mean_c2d = np.mean(all_chi2_dof)
    print(f"\n  Mean χ²/dof across 5 dSphs = {mean_c2d:.2f}")
    if mean_c2d < 2.0:
        print("  → MODEL VALIDATED: consistent good fits across multiple dSphs")
    elif mean_c2d < 3.0:
        print("  → MODEL ACCEPTABLE: generally good fits, not a fluke")
    else:
        print("  → MODEL MARGINAL: mixed results, further investigation needed")
    print()

    # ── Across-dSph comparison: NFW vs MAP_SIDM_FB ──
    print("  ── Improvement over NFW ──")
    for dname in dsph_names:
        c_nfw = grand[dname]['NFW']['chi2_iso']
        c_map = grand[dname]['MAP_SIDM_FB']['best_chi2']
        pct = (c_nfw - c_map)/c_nfw * 100 if c_nfw > 0 else 0
        print(f"    {dname:12s}: NFW χ² = {c_nfw:7.1f}  →  MAP+fb+β χ² = {c_map:7.1f}"
              f"  (Δ = {c_nfw-c_map:+.1f}, {pct:+.0f}%)")
    print()

    # ═══════════════════════════════════════════════════════════════
    #  Plot 1: Multi-panel σ_los (one per dSph)
    # ═══════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    axes_flat = axes.flatten()
    colors = {'NFW': 'gray', 'FB': '#607D8B',
              'BP1_SIDM': '#2196F3', 'BP1_SIDM_FB': '#1565C0',
              'MAP_SIDM': '#E91E63', 'MAP_SIDM_FB': '#C62828'}
    lstyles = {'NFW': ':', 'FB': '--',
               'BP1_SIDM': '-', 'BP1_SIDM_FB': '--',
               'MAP_SIDM': '-', 'MAP_SIDM_FB': '--'}

    for idx, dname in enumerate(dsph_names):
        ax = axes_flat[idx]
        d = DSPHS[dname]
        R_obs = np.array(d['R_pc']) / 1000.0
        s_obs = np.array(d['sigma_los'])
        s_err = np.array(d['sigma_err'])

        ax.errorbar(R_obs*1000, s_obs, yerr=s_err, fmt='ko', ms=5, capsize=3,
                     label='Data', zorder=10)

        for tag in ['NFW', 'MAP_SIDM_FB', 'BP1_SIDM_FB']:
            r = grand[dname][tag]
            ra_str = f"r_a={r['best_ra']:.1f}" if r['best_ra'] < 49 else r"$\beta$=0"
            lbl = f"{tag.replace('_','+')}: $\\chi^2$={r['best_chi2']:.1f} ({ra_str})"
            ax.plot(r['R_proj']*1000, r['best_slos'],
                    color=colors[tag], ls=lstyles[tag], lw=2,
                    label=lbl, alpha=0.9)

        ax.axvline(d['r_half_pc'], color='gray', ls=':', lw=0.7, alpha=0.4)
        ax.set_xlabel('$R$ [pc]', fontsize=10)
        ax.set_ylabel(r'$\sigma_{\rm los}$ [km/s]', fontsize=10)
        ax.set_title(f"{dname}  (N={d['n_stars']})", fontsize=12, fontweight='bold')
        ax.legend(fontsize=7, loc='upper right')
        ax.set_xlim(0, d['R_pc'][-1]*1.4)
        ax.set_ylim(0, max(s_obs)*1.5)

    # Last panel: bar chart summary
    ax = axes_flat[5]
    x_pos = np.arange(len(dsph_names))
    w = 0.28
    nfw_bars = [grand[n]['NFW']['chi2_iso'] / grand[n]['NFW']['ndof']
                for n in dsph_names]
    map_iso  = [grand[n]['MAP_SIDM_FB']['chi2_iso'] / grand[n]['MAP_SIDM_FB']['ndof']
                for n in dsph_names]
    map_best = [grand[n]['MAP_SIDM_FB']['best_chi2'] / grand[n]['MAP_SIDM_FB']['ndof']
                for n in dsph_names]

    ax.bar(x_pos - w, nfw_bars, w, color='gray', alpha=0.5, label='NFW (β=0)')
    ax.bar(x_pos,     map_iso,  w, color='#E91E63', alpha=0.3, label='MAP+fb (β=0)')
    ax.bar(x_pos + w, map_best, w, color='#C62828', alpha=0.9, label='MAP+fb+β')

    for i, (v1, v2, v3) in enumerate(zip(nfw_bars, map_iso, map_best)):
        ax.text(i-w, v1+0.3, f'{v1:.1f}', ha='center', fontsize=6.5, color='gray')
        ax.text(i+w, v2+0.3 if v3<v2 else v3+0.3,
                f'{v3:.2f}', ha='center', fontsize=7, fontweight='bold', color='#C62828')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(dsph_names, fontsize=9, rotation=15)
    ax.set_ylabel(r'$\chi^2/\mathrm{dof}$', fontsize=11)
    ax.set_title(r'$\chi^2/\mathrm{dof}$ comparison', fontsize=12, fontweight='bold')
    ax.axhline(1.0, color='green', ls='--', lw=1, alpha=0.6, label=r'$\chi^2/\mathrm{dof}=1$')
    ax.axhline(2.0, color='orange', ls='--', lw=1, alpha=0.5)
    ax.legend(fontsize=8, loc='upper left')
    ax.set_ylim(0, max(nfw_bars)*1.2)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'multi_dsph_jeans.png')
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # ═══════════════════════════════════════════════════════════════
    #  Plot 2: χ² vs r_a (one curve per dSph, MAP SIDM+FB only)
    # ═══════════════════════════════════════════════════════════════
    fig2, ax2 = plt.subplots(1, 1, figsize=(9, 6))
    dsph_colors = {'Fornax': '#E91E63', 'Sculptor': '#2196F3',
                   'Draco': '#4CAF50', 'Carina': '#FF9800', 'Sextans': '#9C27B0'}
    for dname in dsph_names:
        r = grand[dname]['MAP_SIDM_FB']
        ndof = r['ndof']
        c2norm = np.array(r['chi2_ra']) / ndof
        ax2.plot(r_a_scan, c2norm, 'o-', color=dsph_colors[dname],
                 lw=2, ms=4, label=f"{dname} (dof={ndof})")
    ax2.axhline(1.0, color='green', ls='--', lw=1, alpha=0.6)
    ax2.axhline(2.0, color='orange', ls='--', lw=1, alpha=0.5)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'Anisotropy radius $r_a$ [kpc]', fontsize=12)
    ax2.set_ylabel(r'$\chi^2/\mathrm{dof}$', fontsize=12)
    ax2.set_title('MAP SIDM+Feedback: $\\chi^2$/dof vs anisotropy radius', fontsize=13)
    ax2.legend(fontsize=10)
    ax2.set_xlim(0.25, 60)
    ax2.set_ylim(0, None)
    fig2.tight_layout()
    fig2_path = os.path.join(out_dir, 'multi_dsph_chi2_vs_ra.png')
    fig2.savefig(fig2_path, dpi=150)
    plt.close(fig2)
    print(f"  Saved: {fig2_path}")
    print()
    print("  Done.")


if __name__ == '__main__':
    main()
