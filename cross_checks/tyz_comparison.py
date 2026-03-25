#!/usr/bin/env python3
"""
V10 — v35_tyz_comparison.py
============================
Compare our VPM solver against:
  1. Exact Born transfer cross section (numerical integral, Majorana)
  2. Classical geometric estimate

This validates VPM in the Born regime (λ ≪ 1) where analytics are exact,
and shows the resonant regime (λ ~ 1-30) where full partial-wave analysis
(VPM) is essential and fitting formulas break down.

Output: table + PNG figure.

References:
  - Born limit: standard first-order scattering amplitude
  - Tulin & Yu, Phys. Reports 730, 1 (2018), §III
  - Tulin, Yu, Zurek, PRD 87, 115007 (2013)
"""
# === path setup (auto-generated) ================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# =================================================================

import sys, os, math
import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(__file__))
from v22_raw_scan import sigma_T_vpm, C_KM_S
from config_loader import load_config

# Load BP1 from global config
_CFG = load_config(__file__)
from global_config import GC

# ── Constants (sourced from global_config.json) ──
_PC = GC.physical_constants()
GEV2_TO_CM2 = _PC["GEV2_to_cm2"]
GEV_IN_G    = _PC["GeV_in_g"]
_BP1 = GC.benchmark("BP1")
BP1_M_CHI = _BP1["m_chi_GeV"]
BP1_M_PHI = _BP1["m_phi_MeV"] * 1e-3   # → GeV
BP1_ALPHA = _BP1["alpha"]


def sigma_T_born_majorana(m_chi_GeV, m_phi_GeV, alpha, v_km_s):
    """
    Exact Born transfer cross section for identical Majorana fermions.
    Numerically integrates the symmetrized Born amplitude squared
    with (1 - cosθ) weight.

    For Majorana ψ: spin-1/2 identical fermions, so
      dσ/dΩ = ¼[3|f(θ)-f(π-θ)|² + |f(θ)+f(π-θ)|²]
    where f(θ) = -2μα/(q² + m_φ²) is the Yukawa Born amplitude,
    q = 2k sin(θ/2), and the 3/1 are triplet/singlet spin weights.
    """
    v = v_km_s / 299792.458  # v/c, dimensionless
    mu = m_chi_GeV / 2.0
    k = mu * v
    if k < 1e-30:
        return 0.0

    def integrand(cos_theta):
        t_fwd = 2.0 * k**2 * (1.0 - cos_theta) + m_phi_GeV**2
        t_bwd = 2.0 * k**2 * (1.0 + cos_theta) + m_phi_GeV**2
        f_fwd = 2.0 * mu * alpha / t_fwd
        f_bwd = 2.0 * mu * alpha / t_bwd
        # Majorana: triplet (antisymmetric spatial, weight 3) + singlet (weight 1)
        dsigma = 0.25 * (3.0 * (f_fwd - f_bwd)**2 + (f_fwd + f_bwd)**2)
        return (1.0 - cos_theta) * dsigma * 2.0 * math.pi

    result, _ = quad(integrand, -1.0, 1.0, limit=200)
    # Factor 1/2 for identical particles (avoid double-counting pairs)
    sigma_cm2 = 0.5 * result * GEV2_TO_CM2
    return sigma_cm2 / (m_chi_GeV * GEV_IN_G)


# ──────────────────────────────────────────────────
#  Test points
# ──────────────────────────────────────────────────

TEST_POINTS = [
    # Born regime: λ = α m_χ / m_φ < 1
    ("Born-1",   50.0,  0.100,  1.0e-4,  100.0),
    ("Born-2",   20.0,  0.050,  5.0e-5,   50.0),
    ("Born-3",   10.0,  0.030,  1.0e-4,   30.0),
    ("Born-4",  100.0,  0.200,  5.0e-5,  200.0),
    ("Born-5",   30.0,  0.080,  2.0e-4,  100.0),
    # Transition: λ ~ 1
    ("Trans-1",  20.0,  0.020,  5.0e-4,   50.0),
    ("Trans-2",  50.0,  0.050,  8.0e-4,  100.0),
    # Resonant regime: λ ~ 2-10
    ("Res-1",    BP1_M_CHI, BP1_M_PHI, BP1_ALPHA, 30.0),  # BP1 relic
    ("Res-2",    50.0,  0.015,  1.5e-3,   50.0),
    ("Res-3",    30.0,  0.012,  2.0e-3,  100.0),
    # Deep resonant / classical: λ > 10
    ("Deep-1",  100.0,  0.00915,2.84e-3,  30.0),  # Free best-fit (λ=31)
    ("Deep-2",  100.0,  0.010,  5.0e-3,   30.0),
]


def run():
    print("=" * 85)
    print("  V10 — v35: VPM vs Born Analytic (Majorana)")
    print("  Born formula: exact numerical integral of symmetrized amplitude")
    print("=" * 85)

    header = (f"  {'Label':10s}  {'λ':>6s}  {'κ':>8s}  "
              f"{'VPM σ/m':>12s}  {'Born σ/m':>12s}  {'VPM/Born':>9s}")
    print(header)
    print("  " + "-" * 80)

    labels, vpm_vals, born_vals, lam_vals = [], [], [], []

    for label, mc, mp, al, v in TEST_POINTS:
        lam = al * mc / mp
        vc = v / 299792.458
        kappa = (mc / 2.0) * vc / mp

        sig_vpm  = sigma_T_vpm(mc, mp, al, v)
        sig_born = sigma_T_born_majorana(mc, mp, al, v)

        ratio = sig_vpm / sig_born if sig_born > 0 else float('inf')

        labels.append(label)
        vpm_vals.append(sig_vpm)
        born_vals.append(sig_born)
        lam_vals.append(lam)

        regime_tag = ""
        if lam < 1:
            regime_tag = " ← should agree"
        elif lam > 10:
            regime_tag = " ← resonant, Born invalid"

        print(f"  {label:10s}  {lam:6.2f}  {kappa:8.4f}  "
              f"{sig_vpm:12.4e}  {sig_born:12.4e}  {ratio:9.3f}{regime_tag}")

    # Summary by regime
    born_ratios = [vpm_vals[i]/born_vals[i] for i in range(len(labels))
                   if lam_vals[i] < 1 and born_vals[i] > 0]
    trans_ratios = [vpm_vals[i]/born_vals[i] for i in range(len(labels))
                    if 0.5 <= lam_vals[i] <= 2 and born_vals[i] > 0]

    print()
    print("  Summary:")
    if born_ratios:
        print(f"    Born regime (λ<1):   VPM/Born = {np.mean(born_ratios):.3f} "
              f"± {np.std(born_ratios):.3f}  (expect ≈1)")
    if trans_ratios:
        print(f"    Transition (λ~1):    VPM/Born = {np.mean(trans_ratios):.3f} "
              f"± {np.std(trans_ratios):.3f}  (deviations expected)")
    print()

    # ── Velocity scan at BP1 ──
    bp1_lam = BP1_ALPHA * BP1_M_CHI / BP1_M_PHI
    print(f"  Velocity scan at BP1 (m_χ={BP1_M_CHI} GeV, m_φ={BP1_M_PHI*1e3:.2f} MeV, α={BP1_ALPHA:.4e}, λ={bp1_lam:.2f}):")
    print(f"  {'v [km/s]':>10s}  {'VPM':>12s}  {'Born':>12s}  {'VPM/Born':>9s}")
    print("  " + "-" * 50)
    v_scan = [10, 20, 30, 50, 100, 200, 500, 1000, 2000]
    bp1_vpm, bp1_born, bp1_v = [], [], []
    for v in v_scan:
        sv = sigma_T_vpm(BP1_M_CHI, BP1_M_PHI, BP1_ALPHA, float(v))
        sb = sigma_T_born_majorana(BP1_M_CHI, BP1_M_PHI, BP1_ALPHA, float(v))
        r = sv / sb if sb > 0 else float('inf')
        bp1_vpm.append(sv)
        bp1_born.append(sb)
        bp1_v.append(v)
        print(f"  {v:10d}  {sv:12.4e}  {sb:12.4e}  {r:9.3f}")

    # ── Plot ──
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

        # Panel 1: VPM vs Born scatter
        colors = ['tab:blue' if l < 1 else 'tab:orange' if l < 5 else 'tab:red'
                  for l in lam_vals]
        for i in range(len(labels)):
            c = colors[i]
            lab = None
            if lam_vals[i] < 1 and not any(lam_vals[j] < 1 for j in range(i)):
                lab = r'Born ($\lambda<1$)'
            elif 1 <= lam_vals[i] < 5 and not any(1 <= lam_vals[j] < 5 for j in range(i)):
                lab = r'Transition ($\lambda\sim 1$–$5$)'
            elif lam_vals[i] >= 5 and not any(lam_vals[j] >= 5 for j in range(i)):
                lab = r'Resonant ($\lambda>5$)'
            ax1.scatter(born_vals[i], vpm_vals[i], c=c, s=80, zorder=3, label=lab)
            ax1.annotate(labels[i], (born_vals[i], vpm_vals[i]), fontsize=7,
                        xytext=(5, 5), textcoords='offset points')

        all_v = [x for x in vpm_vals + born_vals if x > 0]
        lo, hi = min(all_v) * 0.2, max(all_v) * 5
        ax1.plot([lo, hi], [lo, hi], 'k--', alpha=0.5, label='Perfect agreement')
        ax1.set_xscale('log'); ax1.set_yscale('log')
        ax1.set_xlim(lo, hi); ax1.set_ylim(lo, hi)
        ax1.set_xlabel(r'Born analytic $\sigma/m$ [cm$^2$/g]', fontsize=11)
        ax1.set_ylabel(r'VPM $\sigma/m$ [cm$^2$/g]', fontsize=11)
        ax1.set_title('VPM vs Born: multi-point', fontsize=12)
        ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')

        # Panel 2: BP1 velocity scan
        ax2.loglog(bp1_v, bp1_vpm, 'o-', color='tab:blue', label='VPM (full partial-wave)')
        ax2.loglog(bp1_v, bp1_born, 's--', color='tab:orange', label='Born analytic')
        ax2.set_xlabel('$v$ [km/s]', fontsize=11)
        ax2.set_ylabel(r'$\sigma/m$ [cm$^2$/g]', fontsize=11)
        ax2.set_title(r'BP1: VPM vs Born ($\lambda=2.19$)', fontsize=12)
        ax2.legend(fontsize=10); ax2.grid(True, alpha=0.3)

        out = os.path.join(os.path.dirname(__file__), "v35_tyz_comparison.png")
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        print(f"\n  Saved figure -> {out}")
        plt.close(fig)
    except ImportError:
        print("\n  matplotlib not available, skipping plot.")

    print("=" * 85)


if __name__ == "__main__":
    run()
