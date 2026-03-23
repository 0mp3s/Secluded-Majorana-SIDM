#!/usr/bin/env python3
"""
model_validations/ufd_crater/predict_ufd.py
============================================
§7.5 — Ultra-Faint Dwarf (UFD) core size predictions + Crater II case study.

Physics:
  SIDM thermalisation produces cores. The core radius r_1 satisfies
      ρ(r_1) × (σ/m) × σ_v × t_age ≈ 1
  where ρ is the NFW density profile (Kaplinghat+2016).

  For UFDs, the key question is whether σ/m at v ~ 3–8 km/s is large enough
  to produce observable cores within a Hubble time.

  Crater II is a special case: very extended (r_half ≈ 1066 pc) but extremely
  low σ_v ≈ 2.7 km/s. In CDM, it should have been tidally disrupted.
  SIDM core formation + tidal stripping may explain its properties.

Targets (Simon+2019, McConnachie+2012, Caldwell+2017):
  - Classical dSphs: Draco, Sculptor, Fornax, Carina
  - UFDs: Tucana II, Segue 1, Reticulum II, Tucana III
  - Special: Crater II

Output: model_validations/ufd_crater/output/
"""
import sys, os, csv, math, json
import numpy as np

# ── path bootstrap ──
_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ── constants ──
G_N = 4.302e-6           # kpc (km/s)² / M_sun
SEC_PER_GYR = 3.156e16   # s / Gyr
KPC_CM = 3.086e21         # cm / kpc
MSUN_G = 1.989e33         # g / M_sun
KM_S_CM_S = 1e5           # cm/s per km/s

# ── output directory ──
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Dwarf galaxy data ──
# Sources: Simon+2019 (review), McConnachie+2012, Read+2019, Caldwell+2017
# M_200 and c_200: from abundance matching (Moster+2013) or Jeans modeling (Read+2019)
# sigma_v: stellar LOS velocity dispersion (proxy for DM velocity)
# Comprehensive dwarf data table
GALAXIES = [
    # Classical dSphs (Walker+2009, Read+2019)
    # name, sigma_v [km/s], r_half [pc], M_200 [M_sun], c_200, core_obs, t_age [Gyr]
    ("Fornax",     11.7,  710, 3.16e9, 18,  "YES",       12.0),
    ("Sculptor",    9.2,  283, 1.5e9,  20,  "YES",       12.0),
    ("Draco",       9.1,  221, 1.2e9,  22,  "AMBIGUOUS", 12.0),
    ("Carina",      6.6,  250, 0.8e9,  20,  "YES",       12.0),
    ("Sextans",     7.9,  695, 1.0e9,  18,  "YES",       12.0),
    ("Leo I",       9.2,  251, 1.5e9,  20,  "YES",        7.0),
    ("Leo II",      6.6,  176, 0.8e9,  22,  "YES",       10.0),
    ("UMi",         9.5,  181, 1.5e9,  22,  "AMBIGUOUS", 12.0),

    # Ultra-Faint Dwarfs (Simon+2019)
    # M_200 from abundance matching: M_200 ~ 10^8–10^9 for UFDs
    ("Tucana II",   3.3,   165, 2.0e8,  25,  "UNKNOWN",  12.0),
    ("Segue 1",     3.7,    29, 3.0e8,  28,  "UNKNOWN",  12.0),
    ("Ret II",      3.3,    32, 2.0e8,  25,  "UNKNOWN",  12.0),
    ("Tucana III",  1.5,    44, 5.0e7,  30,  "UNKNOWN",  12.0),
    ("Carina II",   3.4,    76, 2.5e8,  25,  "UNKNOWN",  12.0),
    ("Grus I",      2.9,    62, 1.5e8,  28,  "UNKNOWN",  12.0),

    # Special: Crater II (Torrealba+2016, Caldwell+2017)
    ("Crater II",   2.7, 1066, 5.0e8,  18,  "YES (inferred)", 10.0),
]

# Benchmark points
BPS = [
    {"label": "BP1", "m_chi": 20.69, "m_phi": 11.34e-3, "alpha": 1.048e-3},
    {"label": "BP9", "m_chi": 37.9,  "m_phi": 16.36e-3, "alpha": 1.840e-3},
    {"label": "MAP", "m_chi": 90.64, "m_phi": 13.85e-3, "alpha": 2.546e-2},
]


def nfw_rho(r, rho_s, r_s):
    """NFW density profile [M_sun/kpc^3]."""
    x = r / r_s
    return rho_s / (x * (1 + x)**2)


def nfw_params(M200, c200):
    """Compute NFW rho_s, r_s from M200 and c200."""
    # R200 from M200: M200 = (4/3)π R200^3 × 200 × ρ_crit
    rho_crit = 277.5  # h^2 M_sun/kpc^3, with h=0.674
    h = 0.674
    rho_200 = 200 * rho_crit * h**2  # M_sun/kpc^3
    R200 = (3 * M200 / (4 * math.pi * rho_200))**(1.0/3)
    r_s = R200 / c200
    g_c = math.log(1 + c200) - c200 / (1 + c200)
    rho_s = M200 / (4 * math.pi * r_s**3 * g_c)
    return rho_s, r_s, R200


def find_core_radius(rho_s, r_s, sigma_m, sigma_v_km, t_age_Gyr):
    """
    Find SIDM core radius r_1 where ρ(r_1) × (σ/m) × v_rel × t_age = 1.

    Units: ρ in g/cm³, σ/m in cm²/g, v in cm/s, t in s.
    """
    if sigma_m <= 0:
        return 0.0

    v_rel_cm = sigma_v_km * math.sqrt(2) * KM_S_CM_S  # relative velocity
    t_s = t_age_Gyr * SEC_PER_GYR

    # Target: rho_cgs × sigma_m × v_rel × t = 1
    # rho_cgs = rho_Msun_kpc3 × MSUN_G / KPC_CM^3
    target_rho_cgs = 1.0 / (sigma_m * v_rel_cm * t_s)
    target_rho_Msun = target_rho_cgs * KPC_CM**3 / MSUN_G

    # Binary search for r where NFW ρ(r) = target
    if nfw_rho(1e-4, rho_s, r_s) < target_rho_Msun:
        return 0.0  # sigma/m too low; not even central density is enough

    r_lo, r_hi = 1e-4, 10 * r_s  # kpc
    for _ in range(100):
        r_mid = (r_lo + r_hi) / 2
        if nfw_rho(r_mid, rho_s, r_s) > target_rho_Msun:
            r_lo = r_mid
        else:
            r_hi = r_mid
    return r_mid  # kpc


def n_scatter(rho_s_Msun, sigma_m, sigma_v_km, t_age_Gyr):
    """Number of scatterings: Γ × t = (σ/m) × ρ × v × t."""
    rho_cgs = rho_s_Msun * MSUN_G / KPC_CM**3  # g/cm³
    v_rel_cm = sigma_v_km * math.sqrt(2) * KM_S_CM_S
    t_s = t_age_Gyr * SEC_PER_GYR
    return sigma_m * rho_cgs * v_rel_cm * t_s


def classify(n_scat):
    """Classify gravothermal state."""
    if n_scat < 1:
        return "CUSPY"
    elif n_scat < 100:
        return "CORED"
    else:
        return "COLLAPSE"


def crater_ii_analysis(bp, rho_s, r_s):
    """
    Crater II case study: test whether SIDM core + tidal stripping
    can explain r_half ≈ 1066 pc with σ_v ≈ 2.7 km/s.

    CDM problem: M(<r_half) ≈ 580 × 2.7² × 1066 ≈ 4.5e6 M_sun
    → mean density ≈ 0.01 M_sun/pc³ → too diffuse for CDM NFW
    → should be tidally disrupted.

    SIDM resolution: core formation lowers central density without
    reducing total mass → larger r_half at same σ_v.
    Tidal radius estimate: r_t ≈ r_MW × (M_sat/3M_MW)^{1/3}
    """
    sigma_v = 2.7  # km/s
    v_rel = sigma_v * math.sqrt(2)

    sigma_m = sigma_T_vpm(bp["m_chi"], bp["m_phi"], bp["alpha"], v_rel)

    # Core radius
    r_core = find_core_radius(rho_s, r_s, sigma_m, sigma_v, 10.0)

    # Walker+2009 mass estimator
    r_half_kpc = 1.066  # kpc
    M_half = 580 * sigma_v**2 * 1066  # M_sun

    # SIDM core density (isothermal inside r_core)
    if r_core > 0:
        rho_core = nfw_rho(r_core, rho_s, r_s)  # constant inside core
    else:
        rho_core = nfw_rho(1e-3, rho_s, r_s)

    # NFW central density (at 10 pc)
    rho_nfw_center = nfw_rho(0.01, rho_s, r_s)

    # Tidal stripping estimate:
    # Crater II is at D ≈ 117 kpc from MW center (Torrealba+2016)
    # Pericenter distance ≈ 50 kpc (estimated from orbit)
    # Tidal radius: r_t = D_peri × (M_sat / 3 M_MW(<D_peri))^{1/3}
    D_peri = 50.0  # kpc (estimated pericenter)
    M_MW_50 = 4.0e11  # M_sun within 50 kpc (Bland-Hawthorn & Gerhard 2016)
    M_sat = 5.0e8  # M_sun (Crater II halo)
    r_tidal = D_peri * (M_sat / (3 * M_MW_50))**(1.0/3)

    return {
        "sigma_m": sigma_m,
        "r_core_pc": r_core * 1000,
        "rho_core": rho_core,
        "rho_nfw_center": rho_nfw_center,
        "density_reduction": rho_nfw_center / max(rho_core, 1),
        "M_half": M_half,
        "r_tidal_kpc": r_tidal,
        "core_within_tidal": r_core < r_tidal,
    }


def main():
    print("=" * 100)
    print("§7.5 — UFD CORE SIZE PREDICTIONS + CRATER II CASE STUDY")
    print("=" * 100)

    # ── Main prediction table ──
    print(f"\n{'Galaxy':>12s}  {'σ_v':>5s}  {'r_h [pc]':>8s}  {'M200':>10s}  {'c':>4s}  ", end="")
    for bp in BPS:
        print(f"{'σ/m(' + bp['label'] + ')':>10s}  {'r_c(' + bp['label'] + ')':>10s}  {'N(' + bp['label'] + ')':>8s}  {'State':>8s}  ", end="")
    print()
    print("-" * 200)

    all_results = []

    for (name, sigma_v, r_half_pc, M200, c200, core_obs, t_age) in GALAXIES:
        rho_s, r_s, R200 = nfw_params(M200, c200)
        v_rel = sigma_v * math.sqrt(2)

        row = {"name": name, "sigma_v": sigma_v, "r_half_pc": r_half_pc,
               "M200": M200, "c200": c200, "core_obs": core_obs}

        line = f"{name:>12s}  {sigma_v:5.1f}  {r_half_pc:8d}  {M200:10.2e}  {c200:4d}  "

        for bp in BPS:
            sigma_m = sigma_T_vpm(bp["m_chi"], bp["m_phi"], bp["alpha"], v_rel)
            r_core = find_core_radius(rho_s, r_s, sigma_m, sigma_v, t_age)
            r_core_pc = r_core * 1000
            n_scat = n_scatter(rho_s, sigma_m, sigma_v, t_age)
            state = classify(n_scat)

            row[f"{bp['label']}_sigma_m"] = sigma_m
            row[f"{bp['label']}_r_core_pc"] = r_core_pc
            row[f"{bp['label']}_N_scatter"] = n_scat
            row[f"{bp['label']}_state"] = state

            line += f"{sigma_m:10.4f}  {r_core_pc:10.1f}  {n_scat:8.2f}  {state:>8s}  "

        all_results.append(row)
        print(line)

    # ── Summary statistics ──
    print("\n" + "=" * 100)
    print("SUMMARY — Core size predictions")
    print("=" * 100)

    for bp in BPS:
        label = bp["label"]
        print(f"\n  {label}:")

        # Classical dSphs
        classical = [r for r in all_results if r["M200"] >= 5e8
                     and r["name"] not in ("Crater II",)]
        ufd = [r for r in all_results if r["M200"] < 5e8]
        crater = [r for r in all_results if r["name"] == "Crater II"]

        if classical:
            cores = [r[f"{label}_r_core_pc"] for r in classical]
            sigmas = [r[f"{label}_sigma_m"] for r in classical]
            print(f"    Classical dSphs: r_core = {min(cores):.0f} — {max(cores):.0f} pc, "
                  f"σ/m = {min(sigmas):.3f} — {max(sigmas):.3f} cm²/g")

        if ufd:
            cores = [r[f"{label}_r_core_pc"] for r in ufd]
            sigmas = [r[f"{label}_sigma_m"] for r in ufd]
            n_zero = sum(1 for c in cores if c < 1)
            print(f"    UFDs: r_core = {min(cores):.0f} — {max(cores):.0f} pc, "
                  f"σ/m = {min(sigmas):.3f} — {max(sigmas):.3f} cm²/g"
                  f" ({n_zero} with negligible cores)")

        if crater:
            c = crater[0]
            print(f"    Crater II: r_core = {c[f'{label}_r_core_pc']:.0f} pc, "
                  f"σ/m = {c[f'{label}_sigma_m']:.4f} cm²/g, "
                  f"N_scatter = {c[f'{label}_N_scatter']:.2f}")

    # ── Key predictions for paper ──
    print("\n" + "=" * 100)
    print("★ KEY PREDICTIONS FOR §7.5:")
    print("=" * 100)

    # Check universality: do all UFDs get similar r_core/r_s?
    for bp in BPS:
        label = bp["label"]
        ufd_data = [(r["name"], r[f"{label}_r_core_pc"],
                      r["r_half_pc"], r[f"{label}_sigma_m"])
                     for r in all_results if r["M200"] < 5e8]
        if ufd_data and all(d[1] > 0 for d in ufd_data):
            ratios = [d[1] / d[2] for d in ufd_data]
            print(f"\n  {label} — r_core/r_half for UFDs:")
            for d in ufd_data:
                r_ratio = d[1] / d[2] if d[2] > 0 else 0
                print(f"    {d[0]:>12s}: r_core = {d[1]:6.0f} pc, "
                      f"r_half = {d[2]:5d} pc, r_core/r_half = {r_ratio:.2f}, "
                      f"σ/m = {d[3]:.4f}")
            mean_ratio = np.mean(ratios)
            std_ratio = np.std(ratios)
            print(f"    → Mean r_core/r_half = {mean_ratio:.2f} ± {std_ratio:.2f}")
            if std_ratio / max(mean_ratio, 0.01) < 0.3:
                print(f"    → ★ UNIVERSAL: scatter < 30% → testable prediction!")
            else:
                print(f"    → NOT universal: scatter = {100*std_ratio/max(mean_ratio,0.01):.0f}%")
        else:
            print(f"\n  {label} — some UFDs have negligible cores (σ/m too low at v ~ 3 km/s)")

    # ── Crater II Special Analysis ──
    print("\n" + "=" * 100)
    print("★ CRATER II CASE STUDY")
    print("=" * 100)

    crater_M200, crater_c200 = 5.0e8, 18
    rho_s_c, r_s_c, R200_c = nfw_params(crater_M200, crater_c200)
    print(f"  Crater II halo: M200 = {crater_M200:.1e} M_sun, c = {crater_c200}, "
          f"R200 = {R200_c:.1f} kpc")
    print(f"  NFW: rho_s = {rho_s_c:.3e} M_sun/kpc³, r_s = {r_s_c:.3f} kpc")
    print(f"  Observed: r_half = 1066 pc, σ_v = 2.7 km/s")

    for bp in BPS:
        result = crater_ii_analysis(bp, rho_s_c, r_s_c)
        print(f"\n  {bp['label']}:")
        print(f"    σ/m(v_rel = {2.7*math.sqrt(2):.1f} km/s) = {result['sigma_m']:.4f} cm²/g")
        print(f"    r_core = {result['r_core_pc']:.0f} pc")
        print(f"    ρ_NFW(10 pc) / ρ_core = {result['density_reduction']:.1f}×")
        print(f"    M(<r_half) Walker = {result['M_half']:.2e} M_sun")
        print(f"    Tidal radius = {result['r_tidal_kpc']:.1f} kpc")
        print(f"    Core within tidal radius? {result['core_within_tidal']}")

        if result['r_core_pc'] > 500:
            print(f"    ★ SIDM core ({result['r_core_pc']:.0f} pc) comparable to r_half (1066 pc)")
            print(f"      → Explains why Crater II survives: core lowers ρ_center by "
                  f"{result['density_reduction']:.0f}×")
            print(f"      → Tidal stripping removes outer material but core persists")
        elif result['r_core_pc'] > 100:
            print(f"    ⚠ SIDM core ({result['r_core_pc']:.0f} pc) smaller than r_half (1066 pc)")
            print(f"      → Partial explanation only")
        else:
            print(f"    ✗ Negligible core — σ/m too low at v ≈ 3 km/s")

    # ── Plot ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: r_core vs sigma_v for all dwarfs
    ax = axes[0]
    colors = {"BP1": "#FF5722", "BP9": "#2196F3", "MAP": "#4CAF50"}
    markers = {"BP1": "o", "BP9": "s", "MAP": "D"}

    for bp in BPS:
        label = bp["label"]
        sigma_vs = [r["sigma_v"] for r in all_results]
        r_cores = [r[f"{label}_r_core_pc"] for r in all_results]
        ax.scatter(sigma_vs, r_cores, c=colors[label], marker=markers[label],
                   s=50, label=label, alpha=0.8, edgecolors='k', linewidths=0.5)

    # Add r_half as reference
    sigma_vs = [r["sigma_v"] for r in all_results]
    r_halfs = [r["r_half_pc"] for r in all_results]
    ax.scatter(sigma_vs, r_halfs, c='gray', marker='x', s=30, alpha=0.5, label=r'$r_{1/2}$ (observed)')

    ax.set_xlabel(r'$\sigma_v$ [km/s]', fontsize=12)
    ax.set_ylabel(r'$r_{\rm core}$ [pc]', fontsize=12)
    ax.set_title('SIDM Core Radius vs Velocity Dispersion', fontsize=13)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_xlim(1, 20)

    # Panel 2: N_scatter classification
    ax = axes[1]
    galaxy_names = [r["name"] for r in all_results]
    x_pos = np.arange(len(galaxy_names))
    width = 0.25

    for i, bp in enumerate(BPS):
        label = bp["label"]
        n_scats = [r[f"{label}_N_scatter"] for r in all_results]
        ax.bar(x_pos + i * width, n_scats, width, label=label,
               color=colors[label], alpha=0.8, edgecolor='k', linewidth=0.3)

    ax.axhline(y=1.0, color='green', linestyle='--', alpha=0.7, label='Core threshold')
    ax.axhline(y=100.0, color='red', linestyle='--', alpha=0.7, label='Collapse threshold')
    ax.set_yscale('log')
    ax.set_ylabel(r'$N_{\rm scatter}$', fontsize=12)
    ax.set_title('Gravothermal Classification', fontsize=13)
    ax.set_xticks(x_pos + width)
    ax.set_xticklabels(galaxy_names, rotation=45, ha='right', fontsize=7)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.2, axis='y')

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, 'ufd_predictions.png')
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\n  Plot saved: {plot_path}")

    # ── Save CSV ──
    csv_path = os.path.join(OUT_DIR, 'ufd_predictions.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        header = ['galaxy', 'sigma_v_km_s', 'r_half_pc', 'M200', 'c200', 'core_obs']
        for bp in BPS:
            label = bp["label"]
            header += [f'{label}_sigma_m', f'{label}_r_core_pc',
                       f'{label}_N_scatter', f'{label}_state']
        w.writerow(header)
        for r in all_results:
            row = [r["name"], r["sigma_v"], r["r_half_pc"], f"{r['M200']:.2e}",
                   r["c200"], r["core_obs"]]
            for bp in BPS:
                label = bp["label"]
                row += [f"{r[f'{label}_sigma_m']:.6f}",
                        f"{r[f'{label}_r_core_pc']:.1f}",
                        f"{r[f'{label}_N_scatter']:.4f}",
                        r[f"{label}_state"]]
            w.writerow(row)
    print(f"  CSV saved: {csv_path}")

    print("\n✓ §7.5 UFD predictions DONE.")


if __name__ == '__main__':
    main()
