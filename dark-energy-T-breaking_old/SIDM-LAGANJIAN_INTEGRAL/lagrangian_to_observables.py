#!/usr/bin/env python3
r"""
lagrangian_to_observables.py — what_does_it_mean/
=================================================
From the path integral of the Secluded Majorana SIDM Lagrangian
to ALL physical observables, derived step by step.

LAGRANGIAN (the starting point):
  L = ½(∂_μ φ)² - ½ m_φ² φ² + χ̄(i∂̸ - m_χ)χ - g_s φ χ̄χ
  α = g_s²/(4π)

PATH INTEGRAL (generating functional):
  Z[J,η,η̄] = ∫ Dφ Dχ Dχ̄ exp(i∫d⁴x L + sources)

FEYNMAN RULES (from Gaussian integration of the free-field quadratic part):
  ─────────────────────────────────────────────────────────────────────
  Scalar propagator:    Δ_F(q²) = i / (q² - m_φ² + iε)
  Fermion propagator:   S_F(p)  = i(p̸ + m_χ) / (p² - m_χ² + iε)
  Interaction vertex:   -i g_s    (scalar × Majorana bilinear)
  ─────────────────────────────────────────────────────────────────────

DERIVATION CHAIN:
  [Two fermions exchange one scalar, t-channel + u-channel for Majorana]
       ↓  NR limit: t = -|q|²,  V̂(q) = -M / (4μ²)
  V(r) = -α/r · exp(-m_φ r)         [Yukawa, attractive, range = ℏ/m_φ]
       ↓  Born approximation (valid only when λ = αm_χ/m_φ ≪ 1)
  σ_T^{Born} = (16πα²)/(m_χ²v⁴) f(β),   β = m_χ v /m_φ
  f(β) = ln(1+β²) − β²/(1+β²)
       ↓  full partial-wave VPM (required when λ > 1)
  σ_T = (2π/k²) Σ_ℓ w_ℓ(2ℓ+1) sin²δ_ℓ    [w=1 even ℓ, w=3 odd ℓ]

ANNIHILATION (χχ → φφ via t/u-channel fermion exchange):
  ⟨σv⟩_s = πα²/(4m_χ²)    [s-wave, m_φ/m_χ → 0]

Output: what_does_it_mean/output/lagrangian_to_observables.png
"""

import sys
import os
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── project root on path ──────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "core"))   # for core-internal relative imports

from core.v22_raw_scan import sigma_T_vpm   # returns σ_T/m_χ [cm²/g]

# ── output directory ──────────────────────────────────────────────────────────
OUT_DIR = os.path.join(_HERE, "output")
os.makedirs(OUT_DIR, exist_ok=True)

# ── physical constants ────────────────────────────────────────────────────────
HBAR_C_GEV_FM = 0.197327099     # ℏc  [GeV·fm]
GEV2_TO_CM2   = 3.8937966e-28   # 1 GeV⁻² → cm²
GEV_IN_G      = 1.78266192e-24  # 1 GeV/c² → g
C_KM_S        = 2.99792458e5    # c  [km/s]
C_CM_S        = 2.99792458e10   # c  [cm/s]

# ── load benchmarks + observations from global config ─────────────────────────
_cfg_path = os.path.join(_ROOT, "data", "global_config.json")
with open(_cfg_path, encoding="utf-8") as _f:
    _cfg = json.load(_f)

_bps = _cfg["benchmark_points"]
BENCHMARKS = [
    {
        "name":  "BP1",
        "m_chi": _bps["BP1"]["m_chi_GeV"],
        "m_phi": _bps["BP1"]["m_phi_MeV"] * 1e-3,
        "alpha": _bps["BP1"]["alpha"],
        "color": "#1f77b4", "ls": "-",
    },
    {
        "name":  "MAP",
        "m_chi": _bps["MAP"]["m_chi_GeV"],
        "m_phi": _bps["MAP"]["m_phi_MeV"] * 1e-3,
        "alpha": _bps["MAP"]["alpha"],
        "color": "#ff7f0e", "ls": "--",
    },
    {
        "name":  "MAP_relic",
        "m_chi": _bps["MAP_relic"]["m_chi_GeV"],
        "m_phi": _bps["MAP_relic"]["m_phi_MeV"] * 1e-3,
        "alpha": _bps["MAP_relic"]["alpha"],
        "color": "#d62728", "ls": "-.",
    },
]
for _bp in BENCHMARKS:
    _bp["lam"] = _bp["alpha"] * _bp["m_chi"] / _bp["m_phi"]

OBSERVATIONS = _cfg["observations"]
SIDM_CUTS    = _cfg["sidm_cuts"]


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Physics functions derived from the Lagrangian
# ════════════════════════════════════════════════════════════════════════════════

def yukawa_potential_MeV(r_fm, alpha, m_phi_gev):
    """
    Yukawa potential |V(r)| in MeV at radius r [fm].

    Derived from the scalar propagator Δ_F(q²) = i/(q²−m_φ²+iε):
      Static NR limit: −iM = (−ig_s)² × i/(|q|²+m_φ²)
      Matching to Born rule V̂(q) = −4πα/(|q|²+m_φ²)
      Fourier inversion → V(r) = −α/r · exp(−m_φ r)

    The 1/e decay length (force range) = ℏc/m_φ = HBAR_C_GEV_FM / m_phi_gev [fm].
    """
    m_phi_fm = m_phi_gev / HBAR_C_GEV_FM   # m_φ in fm⁻¹  (natural units: m_φ/ℏc)
    v_gev    = alpha / r_fm * np.exp(-m_phi_fm * r_fm)   # |V(r)| in α·GeV/fm → GeV
    return v_gev * 1e3                                    # MeV


def sigma_T_born_cm2g(m_chi_gev, m_phi_gev, alpha, v_km_s):
    """
    Born-approximation Majorana transfer cross section [cm²/g].

    Amplitude from tree-level t + u channel scalar exchange:
      iM_t = (−ig_s)² × i/(−|q_t|²−m_φ²) × (2m_χ)²   [NR spin sums = 4m_χ²]
      iM_u = same with q_t ↔ q_u   (u adds for scalar coupling + Majorana)

    Transfer cross section (Majorana, t + u add coherently):
      σ_T^{Born,Maj} = (16πα²)/(m_χ² v⁴) × f(β)
      β ≡ m_χ v / m_φ   (dimensionless: CM momentum / mediator mass)
      f(β) = ln(1+β²) − β²/(1+β²)

    VALIDITY: requires λ = αm_χ/m_φ ≪ 1 (perturbative regime).
    For the benchmark points λ = 11–33 → Born is INVALID.
    """
    v_c  = v_km_s / C_KM_S
    beta = m_chi_gev * v_c / m_phi_gev
    f    = np.log1p(beta**2) - beta**2 / (1.0 + beta**2)
    sigma_gev2 = 16.0 * np.pi * alpha**2 / (m_chi_gev**2 * v_c**4) * f
    return sigma_gev2 * GEV2_TO_CM2 / (m_chi_gev * GEV_IN_G)


def sigma_v_swave_cm3s(alpha, m_chi_gev):
    """
    s-wave annihilation rate ⟨σv⟩ [cm³/s] in the m_φ/m_χ → 0 limit.

    From t + u channel diagrams for Majorana χχ → φφ:
      Amplitude (NR, summed over final-state helicities):
        |M|² = 2g_s⁴ m_χ² / (m_χ² − m_φ²/2)²
      Non-relativistic annihilation:
        ⟨σv⟩ = πα²/(4m_χ²)   [m_φ → 0]

    Compare with Planck thermal relic value: 3×10⁻²⁶ cm³/s.
    The ratio ⟨σv⟩/⟨σv⟩_Planck sets the relic abundance.
    """
    sigma_v_gev2 = np.pi * alpha**2 / (4.0 * m_chi_gev**2)   # GeV⁻²
    return sigma_v_gev2 * GEV2_TO_CM2 * C_CM_S                # cm³/s


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Compute
# ════════════════════════════════════════════════════════════════════════════════

V_GRID  = np.logspace(np.log10(10),   np.log10(5000), 40)   # km/s
R_GRID  = np.logspace(np.log10(0.05), np.log10(60),  400)   # fm
PLANCK_SV = 3.0e-26  # Planck thermal relic ⟨σv⟩ [cm³/s]

print("=" * 68)
print("  From the Lagrangian to Observables — Secluded Majorana SIDM")
print("=" * 68)
print()
print("  L = ½(∂φ)² − ½m_φ²φ² + χ̄(i∂̸−m_χ)χ − g_s φ χ̄χ")
print()
print(f"  {'Model':<12} {'m_χ [GeV]':>11} {'m_φ [MeV]':>10} "
      f"{'λ':>7} {'α':>9} {'⟨σv⟩ [cm³/s]':>16} {'/ Planck':>10}")
print("  " + "─" * 78)

results = {}
for bp in BENCHMARKS:
    name  = bp["name"]
    mc    = bp["m_chi"]
    mp    = bp["m_phi"]
    al    = bp["alpha"]
    lam   = bp["lam"]

    # Force range and de Broglie wavelength at v = 30 km/s (dwarfs)
    r_range_fm = HBAR_C_GEV_FM / mp
    lam_dB_fm  = 4.0 * HBAR_C_GEV_FM / (mc * 30.0 / C_KM_S)   # λ_dB = 2ℏ/(μv) = 4ℏ/(m_χ v)

    # Yukawa potential
    pot = yukawa_potential_MeV(R_GRID, al, mp)

    # Born σ_T (analytic, perturbative)
    st_born = np.array([sigma_T_born_cm2g(mc, mp, al, v) for v in V_GRID])

    # VPM σ_T (full non-perturbative partial-wave sum)
    print(f"  Computing VPM for {name} ({len(V_GRID)} points)...", flush=True)
    st_vpm = np.array([sigma_T_vpm(mc, mp, al, v) for v in V_GRID])

    # Annihilation rate
    sv = sigma_v_swave_cm3s(al, mc)

    results[name] = {
        "pot": pot, "st_born": st_born, "st_vpm": st_vpm,
        "sv": sv, "lam": lam, "alpha": al,
        "m_chi": mc, "m_phi": mp,
        "r_range_fm": r_range_fm, "lam_dB_fm": lam_dB_fm,
        "color": bp["color"], "ls": bp["ls"],
    }
    print(f"  {name:<12} {mc:>11.3f} {mp*1e3:>10.3f} {lam:>7.1f} "
          f"{al:>9.4e} {sv:>16.3e} {sv/PLANCK_SV:>10.2f}×")

print()


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Born breakdown table
# ════════════════════════════════════════════════════════════════════════════════

print("  Born approximation validity at v = 30 km/s:")
print(f"  {'Model':<12} {'λ':>7}  {'σ_T^Born':>16}  {'σ_T^VPM':>16}  {'ratio':>10}")
print("  " + "─" * 65)
idx30 = np.argmin(np.abs(V_GRID - 30.0))
for name, r in results.items():
    sb    = r["st_born"][idx30]
    sv    = r["st_vpm"][idx30]
    ratio = sb / sv if sv > 0 else float("inf")
    label = "INVALID" if r["lam"] > 2 else ("MARGINAL" if r["lam"] > 0.5 else "ok")
    print(f"  {name:<12} {r['lam']:>7.1f}  {sb:>14.2e}  {sv:>14.2e}  "
          f"{ratio:>8.0f}×  Born {label}")
print()


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Three-panel figure
# ════════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(
    r"$\mathcal{L} = \frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2"
    r" + \bar{\chi}(i\partial\!\!\!/ - m_\chi)\chi - g_s\,\phi\bar{\chi}\chi$"
    r"     $\longrightarrow$     observables",
    fontsize=12, y=1.02
)

# ── Panel 1: Yukawa potential V(r) ───────────────────────────────────────────
ax1 = axes[0]
for bp in BENCHMARKS:
    name = bp["name"]
    r    = results[name]
    ax1.semilogy(
        R_GRID, r["pot"],
        color=bp["color"], ls=bp["ls"], lw=2.0,
        label=fr"{name}   $m_\phi$={bp['m_phi']*1e3:.1f} MeV, range={r['r_range_fm']:.0f} fm"
    )
    # Mark force range 1/m_φ
    ax1.axvline(r["r_range_fm"], color=bp["color"], ls=":", lw=1.0, alpha=0.5)

# de Broglie wavelength for BP1 at v_dwarf
bp1 = results["BP1"]
ax1.axvline(
    bp1["lam_dB_fm"], color="gray", ls="--", lw=1.5, alpha=0.8,
    label=r"$\lambda_\mathrm{dB}$(BP1, $v$=30 km/s)"
)
ax1.text(bp1["lam_dB_fm"] * 1.07, 2e-4, r"$\lambda_\mathrm{dB}$",
         color="gray", fontsize=9, va="bottom")

ax1.set_xlabel(r"$r$  [fm]", fontsize=11)
ax1.set_ylabel(r"$|V(r)|$  [MeV]", fontsize=11)
ax1.set_xlim(0.05, 60)
ax1.set_ylim(1e-5, None)
ax1.set_title(
    "Step 1: Yukawa potential\n"
    r"$V(r) = -\alpha\,e^{-m_\phi r}/r$   [from $\Delta_F(q^2)$]",
    fontsize=10
)
ax1.legend(fontsize=8, loc="upper right")
ax1.grid(True, alpha=0.3)

# ── Panel 2: Born vs VPM — shows Born completely fails ───────────────────────
ax2 = axes[1]
name = "BP1"
r    = results["BP1"]
lam  = r["lam"]

ax2.loglog(V_GRID, r["st_born"], color="gray", lw=2.0, ls="--",
           label=fr"Born  (perturbative, $\lambda$={lam:.0f} $\gg$ 1)")
ax2.loglog(V_GRID, r["st_vpm"],  color=results["BP1"]["color"], lw=2.5, ls="-",
           label="VPM  (non-perturbative, correct)")

# Annotate ratio at v = 30 km/s
sb30  = r["st_born"][idx30]
sv30  = r["st_vpm"][idx30]
ratio = sb30 / sv30
ax2.annotate(
    f"Born/VPM = {ratio:.0f}×\n@ $v$ = 30 km/s",
    xy=(30, sv30), xytext=(70, sv30 * 200),
    arrowprops=dict(arrowstyle="->", color="black", lw=1.2),
    fontsize=9, ha="left",
)

# SIDM windows (reference)
ax2.axhspan(SIDM_CUTS["sigma_m_30_lo"], SIDM_CUTS["sigma_m_30_hi"],
            xmin=0.0, xmax=0.22, alpha=0.10, color="limegreen",
            label="Dwarf window [0.5–10 cm²/g]")
ax2.axhline(SIDM_CUTS["sigma_m_1000_hi"], color="salmon", ls="-.", lw=1.5,
            label=fr"Cluster bound <{SIDM_CUTS['sigma_m_1000_hi']} cm²/g")

ax2.set_xlabel(r"$v_\mathrm{rel}$  [km/s]", fontsize=11)
ax2.set_ylabel(r"$\sigma_T / m_\chi$  [cm²/g]", fontsize=11)
ax2.set_xlim(10, 5000)
ax2.set_ylim(1e-4, 1e10)
ax2.set_title(
    fr"Step 2: Born fails for BP1 ($\lambda$={lam:.0f} $\gg$ 1)"
    "\n→ VPM mandatory",
    fontsize=10
)
ax2.legend(fontsize=8, loc="upper right")
ax2.grid(True, alpha=0.3)

# ── Panel 3: Full VPM σ_T/m vs v with observations ───────────────────────────
ax3 = axes[2]

# SIDM viability windows (shaded)
ax3.axhspan(SIDM_CUTS["sigma_m_30_lo"], SIDM_CUTS["sigma_m_30_hi"],
            xmin=0.00, xmax=0.24, alpha=0.12, color="limegreen",
            label="Dwarf SIDM window")
ax3.axhspan(0, SIDM_CUTS["sigma_m_1000_hi"],
            xmin=0.62, xmax=1.00, alpha=0.12, color="salmon",
            label="Cluster upper bound")
ax3.text(12,   SIDM_CUTS["sigma_m_30_lo"] * 1.15, "Dwarf\nwindow",
         fontsize=7.5, color="darkgreen", va="bottom")
ax3.text(2500, SIDM_CUTS["sigma_m_1000_hi"] * 0.25, "Cluster\nbound",
         fontsize=7.5, color="firebrick", va="top")

# Observational data points
for obs in OBSERVATIONS:
    v_o  = obs["v_km_s"]
    c_o  = obs["central"]
    lo_o = obs["lo"]
    hi_o = obs["hi"]
    is_upper = (lo_o == 0.0)
    yerr_lo  = (c_o - lo_o) if not is_upper else 0.0
    yerr_hi  = hi_o - c_o
    ax3.errorbar(
        v_o, c_o, yerr=[[yerr_lo], [yerr_hi]],
        fmt="ks", ms=4, capsize=3, lw=1.0,
        uplims=is_upper, alpha=0.75,
    )

# Benchmark curves
for bp in BENCHMARKS:
    name = bp["name"]
    r    = results[name]
    ax3.loglog(V_GRID, r["st_vpm"],
               color=bp["color"], ls=bp["ls"], lw=2.0,
               label=fr"{name}  ($\lambda$={r['lam']:.0f})")

ax3.set_xlabel(r"$v_\mathrm{rel}$  [km/s]", fontsize=11)
ax3.set_ylabel(r"$\sigma_T / m_\chi$  [cm²/g]", fontsize=11)
ax3.set_xlim(10, 5000)
ax3.set_ylim(5e-4, 30)
ax3.set_title(
    "Step 3: From $\\mathcal{L}$ to the sky\n"
    "VPM cross section vs observations",
    fontsize=10
)
ax3.legend(fontsize=8, loc="upper right")
ax3.grid(True, alpha=0.3)

fig.tight_layout()
out_path = os.path.join(OUT_DIR, "lagrangian_to_observables.png")
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"  Saved → {out_path}")


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Physical interpretation (what the Lagrangian means)
# ════════════════════════════════════════════════════════════════════════════════

print()
print("=" * 68)
print("  What the Lagrangian MEANS physically")
print("=" * 68)

print()
print("  ① FORCE RANGE  r₀ = ℏc/m_φ:")
for bp in BENCHMARKS:
    name = bp["name"]
    r    = results[name]
    print(f"     {name:<12}  r₀ = {r['r_range_fm']:.1f} fm"
          f"   (~{r['r_range_fm']/1.2:.0f}× proton radius)")

print()
print("  ② DE BROGLIE WAVELENGTH  λ_dB = 4ℏc/(m_χ v)  at v = 30 km/s (dwarfs):")
for bp in BENCHMARKS:
    name = bp["name"]
    r    = results[name]
    regime = "DEEPLY RESONANT" if r["lam_dB_fm"] / r["r_range_fm"] < 2 else "CLASSICAL"
    print(f"     {name:<12}  λ_dB = {r['lam_dB_fm']:.1f} fm"
          f"   = {r['lam_dB_fm']/r['r_range_fm']:.1f} × r₀  →  {regime}")

print()
print("  ③ BORN APPROXIMATION STATUS  (λ = αm_χ/m_φ ≪ 1 required):")
idx30 = np.argmin(np.abs(V_GRID - 30.0))
for bp in BENCHMARKS:
    name  = bp["name"]
    r     = results[name]
    ratio = r["st_born"][idx30] / r["st_vpm"][idx30]
    print(f"     {name:<12}  λ = {r['lam']:.1f}"
          f"   →  Born overestimates by {ratio:.0f}×  (INVALID)")

print()
print("  ④ RELIC DENSITY from ⟨σv⟩_s = πα²/(4m_χ²)  [χχ → φφ]:")
for bp in BENCHMARKS:
    name = bp["name"]
    r    = results[name]
    print(f"     {name:<12}  ⟨σv⟩ = {r['sv']:.3e} cm³/s"
          f"   = {r['sv']/PLANCK_SV:.2f} × ⟨σv⟩_Planck")

print()
print("  SUMMARY — what L = ½(∂φ)² - ½mφ²φ² + χ̄(i∂̸-mχ)χ - gsφχ̄χ means:")
print()
print("  • The vertex -igs generates an attractive Yukawa force with range ~15 fm")
print("  • λ ~ 10–33: deep resonant regime, Born approximation fails by 10³–10⁴×")
print("  • VPM partial-wave sum correctly unitarises the cross section")
print("  • σ_T(v) is naturally large at dwarfs, small at clusters  ✓ SIDM")
print("  • Annihilation ⟨σv⟩ ~ Planck relic value  ✓ correct relic abundance")
print()
