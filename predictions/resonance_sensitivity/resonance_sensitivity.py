#!/usr/bin/env python3
"""
predictions/resonance_sensitivity/resonance_sensitivity.py
==========================================================
Resonance-enhanced concentration sensitivity: the Diversity Mechanism.

Near a Ramsauer-Townsend resonance (λ ≈ λ_c), σ/m(v) is extremely
sensitive to small changes in velocity.  Since halo concentration c
controls both ρ_s (→ r_1) and σ_v (→ v_rel → σ/m), the natural
~0.16 dex cosmological scatter in c translates into large scatter
in r_core — a diversity mechanism unique to near-resonant couplings.

This script computes:
  1. For each galaxy × BP × c_factor:
     - NFW params (ρ_s, r_s) from M_200, c
     - σ_v from c-dependent virial relation
     - σ/m(v_rel) from VPM solver
     - r_1 (thermalization radius) from Kaplinghat+2016 criterion
  2. Diversity metric: Δr_core / r_core for each galaxy × BP
  3. Comparison: BP1 (λ=1.91, near resonance) vs MAP (λ=48.6, deep resonant)

Physics:
  BP1 sits at λ = 1.91, just past the first critical coupling λ_c = 1.68.
  The s-wave phase shift δ₀ has crossed π/2, placing BP1 in a
  Ramsauer-Townsend minimum where σ_att/σ_rep = 0.017 at v=30 km/s.
  Small changes in v (from Δc) cause large Δσ/m → large Δr_core.

  MAP at λ = 48.6 is deep in the resonant regime where many partial
  waves contribute.  The cross section varies slowly → small Δr_core.

Output:
  - Table: galaxy × BP × c_factor → r_core, σ/m, diversity metric
  - 2-panel figure: r_core vs c_factor for BP1 (steep) vs MAP (flat)
  - CSV with all results
"""
import sys, os, math, csv, time
import numpy as np

# ---------- path bootstrap ----------
_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
# ------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config_loader import load_config
from global_config import GC
from v22_raw_scan import sigma_T_vpm

# JIT warmup
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

# ---- constants ----
_RHO_CRIT    = GC.cosmological_constants()["rho_crit_Msun_kpc3"]  # M_sun/kpc³ (h=0.674)
G_N          = 4.302e-6       # kpc (km/s)² / M_sun
KPC_CM       = 3.086e21       # cm
MSUN_G       = 1.989e33       # g
SEC_PER_GYR  = 3.156e16       # s
KM_S_TO_CM_S = 1e5

# =========================================================================
#  Galaxy data: 8 classical dSphs + 6 UFDs
#  (name, σ_v [km/s], r_half [pc], M_200 [M_sun], c_200_fid, t_age [Gyr])
# =========================================================================
GALAXIES = [
    # Classical dSphs (Walker+2009, Read+2019)
    ("Fornax",      11.7,  710, 3.16e9, 18, 12.0),
    ("Sculptor",     9.2,  283, 1.50e9, 20, 12.0),
    ("Draco",        9.1,  221, 1.20e9, 22, 12.0),
    ("Carina",       6.6,  250, 0.80e9, 20, 12.0),
    ("Sextans",      7.9,  695, 1.00e9, 18, 12.0),
    ("Leo I",        9.2,  251, 1.50e9, 20,  7.0),
    ("Leo II",       6.6,  176, 0.80e9, 22, 10.0),
    ("UMi",          9.5,  181, 1.50e9, 22, 12.0),
    # UFDs (Simon+2019)
    ("Tucana II",    3.3,  165, 2.0e8,  25, 12.0),
    ("Segue 1",      3.7,   29, 3.0e8,  28, 12.0),
    ("Ret II",       3.3,   32, 2.0e8,  25, 12.0),
    ("Carina II",    3.4,   76, 2.5e8,  25, 12.0),
    ("Grus I",       2.9,   62, 1.5e8,  28, 12.0),
    ("Crater II",    2.7, 1066, 5.0e8,  18, 10.0),
]

# =========================================================================
#  NFW helpers
# =========================================================================

def nfw_rho(r_kpc, rho_s, r_s):
    x = max(r_kpc / r_s, 1e-8)
    return rho_s / (x * (1.0 + x)**2)


def nfw_mass(r_kpc, rho_s, r_s):
    x = r_kpc / r_s
    return 4.0 * math.pi * rho_s * r_s**3 * (math.log(1 + x) - x / (1 + x))


def nfw_params_from_M200_c(M_200, c):
    """Return (ρ_s [M_sun/kpc³], r_s [kpc]) from M_200 and c."""
    rho_crit = _RHO_CRIT
    R_200 = (3.0 * M_200 / (4.0 * math.pi * 200.0 * rho_crit))**(1.0/3.0)
    r_s = R_200 / c
    f_c = math.log(1 + c) - c / (1 + c)
    delta_c = (200.0 / 3.0) * c**3 / f_c
    rho_s = rho_crit * delta_c
    return rho_s, r_s


def nfw_sigma_v(M_200, c):
    """Approximate 1D velocity dispersion [km/s] from M_200, c.

    Uses V_max / sqrt(3) where V_max = sqrt(G M_200 f(x_max) / (r_s x_max))
    and x_max ≈ 2.163.
    """
    rho_s, r_s = nfw_params_from_M200_c(M_200, c)
    x_max = 2.163
    f_xmax = math.log(1 + x_max) - x_max / (1 + x_max)
    V_max = math.sqrt(G_N * 4.0 * math.pi * rho_s * r_s**2 * f_xmax / x_max)
    return V_max / math.sqrt(3.0)


# =========================================================================
#  SIDM core finder
# =========================================================================

def find_r1(rho_s, r_s, sigma_over_m, sigma_v_km_s, t_age_s):
    """Thermalization radius r_1 [kpc]: ρ(r_1) × (σ/m) × v × t = 1."""
    v_cm = sigma_v_km_s * KM_S_TO_CM_S
    target_rho_cgs = 1.0 / (sigma_over_m * v_cm * t_age_s)
    target_rho = target_rho_cgs * KPC_CM**3 / MSUN_G

    r_lo, r_hi = 1e-4, 10.0 * r_s
    if nfw_rho(r_lo, rho_s, r_s) < target_rho:
        return r_lo
    if nfw_rho(r_hi, rho_s, r_s) > target_rho:
        return r_hi
    for _ in range(100):
        mid = 0.5 * (r_lo + r_hi)
        if nfw_rho(mid, rho_s, r_s) > target_rho:
            r_lo = mid
        else:
            r_hi = mid
    return 0.5 * (r_lo + r_hi)


# =========================================================================
#  Main computation
# =========================================================================

def run():
    t0 = time.time()
    cfg = load_config(__file__)

    bps = GC.benchmarks_from_labels(cfg.get('benchmark_labels', []))
    c_factors = cfg.get("c_factors", [0.8, 0.9, 1.0, 1.1, 1.2])
    t_age_gyr = cfg.get("halo_age_Gyr", 10.0)
    out_dir = os.path.join(_DIR, cfg.get("output_dir", "output"))
    os.makedirs(out_dir, exist_ok=True)

    t_age_s = t_age_gyr * SEC_PER_GYR

    hdr = "=" * 76
    print(hdr)
    print("  Resonance-Enhanced Concentration Sensitivity")
    print(f"  {len(GALAXIES)} galaxies × {len(bps)} BPs × {len(c_factors)} c-factors")
    print(f"  c_factors = {c_factors}")
    print(hdr)

    # =====================================================================
    #  MODE A: Fixed σ_v — only ρ_s changes with c (density effect)
    #  MODE B: NFW-derived σ_v — both ρ_s AND v_rel change (full effect)
    # =====================================================================
    # Storage: results[mode][bp_label][gal_name] = [(c_factor, c, sigma_v, sigma_m, r1)]
    results = {"fixed": {}, "nfw": {}}
    all_rows = []

    for mode in ["fixed", "nfw"]:
        mode_label = "FIXED σ_v (density effect only)" if mode == "fixed" \
            else "NFW-DERIVED σ_v (density + velocity effects)"
        print(f"\n{'='*76}")
        print(f"  MODE: {mode_label}")
        print(f"{'='*76}")

        for bp in bps:
            label = bp["label"]
            m_chi = bp["m_chi_GeV"]
            m_phi = bp["m_phi_MeV"] / 1000.0  # MeV → GeV for solver
            alpha = bp["alpha"]
            lam = alpha * m_chi / m_phi

            results[mode][label] = {}
            print(f"\n  {label}: m_χ={m_chi} GeV, m_φ={bp['m_phi_MeV']} MeV, "
                  f"α={alpha:.3e}, λ={lam:.2f}")
            print(f"  {'Galaxy':<12} {'c_fac':>5} {'c':>5} {'σ_v':>7} {'v_rel':>7} "
                  f"{'σ/m':>8} {'r_1':>7} {'r_1/r_h':>7}")
            print("  " + "-" * 68)

            for gal in GALAXIES:
                name, sigma_v_obs, r_half_pc, M_200, c_fid, t_age_g = gal
                r_half_kpc = r_half_pc / 1000.0
                t_s = t_age_g * SEC_PER_GYR
                results[mode][label][name] = []

                for cf in c_factors:
                    c = c_fid * cf
                    rho_s, r_s = nfw_params_from_M200_c(M_200, c)

                    if mode == "fixed":
                        sigma_v = sigma_v_obs
                    else:
                        # σ_v from NFW V_max, scaled so c=c_fid gives σ_v_obs
                        sigma_v_nfw = nfw_sigma_v(M_200, c)
                        sigma_v_fid = nfw_sigma_v(M_200, c_fid)
                        sigma_v = sigma_v_obs * (sigma_v_nfw / sigma_v_fid)

                    v_rel = math.sqrt(2.0) * sigma_v

                    # σ/m at this velocity
                    sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_rel)

                    # Thermalization radius
                    r_1 = find_r1(rho_s, r_s, sigma_m, sigma_v, t_s)

                    results[mode][label][name].append(
                        (cf, c, sigma_v, sigma_m, r_1))

                    ratio = r_1 / r_half_kpc if r_half_kpc > 0 else 0

                    print(f"  {name:<12} {cf:>5.2f} {c:>5.1f} {sigma_v:>7.1f} "
                          f"{v_rel:>7.1f} {sigma_m:>8.3f} {r_1:>7.3f} {ratio:>7.2f}")

                    all_rows.append({
                        "mode": mode, "bp": label, "galaxy": name,
                        "c_factor": cf, "c": c,
                        "sigma_v": sigma_v, "v_rel": v_rel,
                        "sigma_m": sigma_m, "r_1_kpc": r_1,
                        "r_half_kpc": r_half_kpc,
                        "r_1_over_r_half": ratio,
                        "lambda": lam,
                    })

    # =====================================================================
    #  Diversity metric for BOTH modes
    # =====================================================================
    diversity = {}  # diversity[mode][label] = [(name, r1_min, r1_med, r1_max, delta)]

    for mode in ["fixed", "nfw"]:
        mode_label = "FIXED σ_v" if mode == "fixed" else "NFW-DERIVED σ_v"
        diversity[mode] = {}

        print(f"\n{'='*76}")
        print(f"  DIVERSITY METRIC ({mode_label}): "
              f"Δr_1/r_1(median) for ±20% concentration scatter")
        print(f"  {'Galaxy':<12} ", end="")
        for bp in bps:
            print(f"  {'r1_min':>7} {'r1_med':>7} {'r1_max':>7} {'Δr1/r1':>7}", end="")
        print()
        print("  " + "-" * (12 + len(bps) * 32))

        for gal in GALAXIES:
            name = gal[0]
            print(f"  {name:<12} ", end="")
            for bp in bps:
                label = bp["label"]
                entries = results[mode][label][name]
                r1_vals = [e[4] for e in entries]
                r1_min, r1_max = min(r1_vals), max(r1_vals)
                r1_med = [e[4] for e in entries if abs(e[0] - 1.0) < 0.01][0]
                delta = (r1_max - r1_min) / r1_med if r1_med > 0 else 0
                if label not in diversity[mode]:
                    diversity[mode][label] = []
                diversity[mode][label].append((name, r1_min, r1_med, r1_max, delta))
                print(f"  {r1_min:>7.3f} {r1_med:>7.3f} {r1_max:>7.3f} "
                      f"{delta:>7.1%}", end="")
            print()

    # Comparison summary
    print(f"\n{'='*76}")
    print("  COMPARISON: Fixed σ_v vs NFW-derived σ_v")
    print(f"{'='*76}")
    print(f"  {'BP':<6} {'Mode':<12} {'Classical':>12} {'UFDs':>12} {'All':>12}")
    print("  " + "-" * 54)
    for bp in bps:
        label = bp["label"]
        lam = bp["alpha"] * bp["m_chi_GeV"] / (bp["m_phi_MeV"] / 1000)
        for mode in ["fixed", "nfw"]:
            classical = [d[4] for d in diversity[mode][label][:8]]
            ufds = [d[4] for d in diversity[mode][label][8:]]
            all_d = classical + ufds
            mode_tag = "fixed" if mode == "fixed" else "nfw"
            print(f"  {label:<6} {mode_tag:<12} {np.mean(classical):>11.1%} "
                  f"{np.mean(ufds):>11.1%} {np.mean(all_d):>11.1%}")

    # Velocity effect isolation
    print(f"\n  VELOCITY EFFECT (NFW/fixed - 1):")
    print(f"  {'BP':<6} {'λ':>6} {'Classical':>12} {'UFDs':>12} {'All':>12}")
    print("  " + "-" * 48)
    for bp in bps:
        label = bp["label"]
        lam = bp["alpha"] * bp["m_chi_GeV"] / (bp["m_phi_MeV"] / 1000)
        c_fixed = [d[4] for d in diversity["fixed"][label][:8]]
        c_nfw = [d[4] for d in diversity["nfw"][label][:8]]
        u_fixed = [d[4] for d in diversity["fixed"][label][8:]]
        u_nfw = [d[4] for d in diversity["nfw"][label][8:]]
        a_fixed = c_fixed + u_fixed
        a_nfw = c_nfw + u_nfw
        boost_c = np.mean(c_nfw) / np.mean(c_fixed) - 1
        boost_u = np.mean(u_nfw) / np.mean(u_fixed) - 1
        boost_a = np.mean(a_nfw) / np.mean(a_fixed) - 1
        print(f"  {label:<6} {lam:>6.1f} {boost_c:>+11.1%} "
              f"{boost_u:>+11.1%} {boost_a:>+11.1%}")

    # =====================================================================
    #  Write CSV
    # =====================================================================
    csv_path = os.path.join(out_dir, "resonance_sensitivity.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "mode", "bp", "galaxy", "c_factor", "c", "sigma_v", "v_rel",
            "sigma_m", "r_1_kpc", "r_half_kpc", "r_1_over_r_half", "lambda",
        ])
        w.writeheader()
        w.writerows(all_rows)
    print(f"\n  CSV: {csv_path}")

    # =====================================================================
    #  Figures
    # =====================================================================

    # --- Figure 1: r_core vs c_factor — fixed σ_v vs NFW σ_v ---
    n_bp = len(bps)
    fig, axes = plt.subplots(2, n_bp, figsize=(7 * n_bp, 12), squeeze=False)
    classical_names = [g[0] for g in GALAXIES[:8]]
    colors = plt.cm.tab10(np.linspace(0, 1, 8))

    for row, mode in enumerate(["fixed", "nfw"]):
        mode_label = "Fixed $\\sigma_v$" if mode == "fixed" \
            else "NFW-derived $\\sigma_v$"
        for col, bp in enumerate(bps):
            ax = axes[row][col]
            label = bp["label"]
            lam = bp["alpha"] * bp["m_chi_GeV"] / (bp["m_phi_MeV"] / 1000)
            for i, name in enumerate(classical_names):
                entries = results[mode][label][name]
                cfs = [e[0] for e in entries]
                r1s = [e[4] for e in entries]
                ax.plot(cfs, r1s, 'o-', color=colors[i], label=name,
                        markersize=5)
            ax.set_xlabel("$c/c_{\\rm fid}$", fontsize=11)
            ax.set_ylabel("$r_1$ [kpc]", fontsize=11)
            ax.set_title(f"{label} ($\\lambda={lam:.1f}$) — {mode_label}",
                         fontsize=11)
            ax.legend(fontsize=7, ncol=2, loc="best")
            ax.grid(True, alpha=0.3)

    fig.suptitle("Resonance-Enhanced Concentration Sensitivity: "
                 "$r_{\\rm core}$ vs Halo Concentration",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig1_path = os.path.join(out_dir, "resonance_sensitivity_classical.png")
    fig.savefig(fig1_path, dpi=150, bbox_inches='tight')
    fig.savefig(fig1_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Fig1: {fig1_path}")

    # --- Figure 2: Diversity bar chart — fixed vs NFW side by side ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    all_names = [g[0] for g in GALAXIES]
    x = np.arange(len(all_names))
    width = 0.8 / n_bp
    bp_colors = ['#2196F3', '#FF5722', '#4CAF50', '#9C27B0'][:n_bp]

    for ax_idx, mode in enumerate(["fixed", "nfw"]):
        ax = axes[ax_idx]
        mode_title = "Density effect only (fixed $\\sigma_v$)" if mode == "fixed" \
            else "Full effect (NFW-derived $\\sigma_v$)"
        for i, bp in enumerate(bps):
            label = bp["label"]
            lam = bp['alpha'] * bp['m_chi_GeV'] / (bp['m_phi_MeV'] / 1000)
            deltas = [d[4] * 100 for d in diversity[mode][label]]
            ax.bar(x + (i - (n_bp - 1) / 2) * width, deltas, width,
                   label=f"{label} ($\\lambda={lam:.1f}$)",
                   color=bp_colors[i], alpha=0.8)
        ax.set_xlabel("Galaxy", fontsize=11)
        if ax_idx == 0:
            ax.set_ylabel("$\\Delta r_1 / r_1$ [%]  ($\\Delta c/c = \\pm 20\\%$)",
                          fontsize=11)
        ax.set_title(mode_title, fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(all_names, rotation=45, ha='right', fontsize=8)
        ax.axvline(x=7.5, color='gray', ls='--', alpha=0.5,
                   label='Classical | UFD')
        ax.legend(fontsize=9)
        ax.grid(True, axis='y', alpha=0.3)

    fig.suptitle("Core Size Diversity from Cosmological Concentration Scatter",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig2_path = os.path.join(out_dir, "resonance_diversity_metric.png")
    fig.savefig(fig2_path, dpi=150, bbox_inches='tight')
    fig.savefig(fig2_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Fig2: {fig2_path}")

    # --- Figure 3: σ/m(v) with shaded band from c-scatter ---
    fig, axes = plt.subplots(1, n_bp, figsize=(7 * n_bp, 6), squeeze=False)
    v_range = np.logspace(np.log10(1), np.log10(300), 80)

    for ax_idx, bp in enumerate(bps):
        ax = axes[0][ax_idx]
        label = bp["label"]
        m_chi = bp["m_chi_GeV"]
        m_phi = bp["m_phi_MeV"] / 1000.0
        alpha = bp["alpha"]
        lam = alpha * m_chi / m_phi

        sigma_m_curve = [sigma_T_vpm(m_chi, m_phi, alpha, float(v))
                         for v in v_range]
        ax.plot(v_range, sigma_m_curve, 'k-', lw=2, label="$\\sigma/m(v)$")

        # Mark Fornax and typical UFD velocity bands
        fornax_v = 11.7 * math.sqrt(2)
        ufd_v = 3.5 * math.sqrt(2)
        ax.axvspan(ufd_v * 0.7, ufd_v * 1.3, alpha=0.15, color='purple',
                   label=f"UFD $v_{{\\rm rel}}$ band")
        ax.axvspan(fornax_v * 0.7, fornax_v * 1.3, alpha=0.15, color='green',
                   label=f"dSph $v_{{\\rm rel}}$ band")

        # Derivative: d(ln σ/m)/d(ln v) ≈ sensitivity
        log_v = np.log(v_range)
        log_s = np.log(np.maximum(sigma_m_curve, 1e-30))
        dlogS_dlogV = np.gradient(log_s, log_v)

        ax2 = ax.twinx()
        ax2.plot(v_range, np.abs(dlogS_dlogV), 'r--', alpha=0.6, lw=1)
        ax2.set_ylabel("$|d\\ln(\\sigma/m) / d\\ln v|$", color='red',
                       fontsize=10)
        ax2.tick_params(axis='y', labelcolor='red')
        ax2.set_ylim(0, 5)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("$v_{\\rm rel}$ [km/s]", fontsize=12)
        ax.set_ylabel("$\\sigma/m$ [cm²/g]", fontsize=12)
        ax.set_title(f"{label} ($\\lambda = {lam:.1f}$)", fontsize=13)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_xlim(1, 300)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Velocity Sensitivity: $\\sigma/m(v)$ and Log-Slope",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig3_path = os.path.join(out_dir, "sigma_v_sensitivity.png")
    fig.savefig(fig3_path, dpi=150, bbox_inches='tight')
    fig.savefig(fig3_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  Fig3: {fig3_path}")

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f} s")
    print(hdr)


if __name__ == "__main__":
    run()


if __name__ == '__main__':
    try:
        import sys as _sys, os as _os
        _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'core'))
        from tg_notify import notify
        notify("\u2705 resonance_sensitivity done!")
    except Exception:
        pass
