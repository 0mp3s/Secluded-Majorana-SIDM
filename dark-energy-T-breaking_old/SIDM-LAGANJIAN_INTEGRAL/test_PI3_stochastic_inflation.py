"""
Test PI-3: Stochastic Inflation — P(θ_i) Distribution
====================================================================
Core question:
  Is θ_i ~ 2 rad (required by PI-7 → H₀ = 67.4) a NATURAL outcome of
  de Sitter quantum diffusion during inflation, or is it a free parameter?

Physics:
  During inflation, every scalar field σ undergoes quantum diffusion:
    ⟨(δθ)²⟩ = (H_inf / 2πf)²   per e-fold   [Starobinsky 1986]

  Fokker-Planck equation in e-fold time N:
    ∂P/∂N = D ∂²P/∂θ² + ∂/∂θ[ V'(θ)/(3H_inf² f²) · P ]
    D = (H_inf / 2πf)²

  Since V(θ) = Λ_d⁴(1−cosθ) with Λ_d ~ 2 meV  and  H_inf ~ 10¹³ GeV:
    drift / diffusion  ~  Λ_d⁴ / (3 H_inf² f²)  ≪  D
  → PURE DIFFUSION regime: drift term negligible

  After N_e e-folds, starting from θ = 0 (shift symmetry):
    σ_rms = (H_inf / 2πf) × √N_e
    P(θ_i) ≈ Gaussian(0, σ_rms)   if σ_rms ≪ π
    P(θ_i) ≈ Uniform[0, π]         if σ_rms ≫ π  (ergodic)

Three results:
  1. At CMB-allowed H_inf ~ 10¹³ GeV: σ_rms ~ 4×10⁻⁶ rad → θ_i is NOT set
     by inflation. Extreme fine-tuning required.
  2. Natural θ_i needs H_inf ~ 4.7×10¹⁸ GeV — 5 orders above CMB bound.
  3. V(θ) drift term during inflation is negligible by ~40 orders of magnitude.

Implication for PI-7:
  H₀ = 67.4 is a CONSISTENCY CHECK given θ_i, not a parameter-free prediction.
  → PI-1 (thermal attractor in V_eff(T)) is the next key test.

MAP PARAMETERS (fixed):
  m_χ = 94.07 MeV,  m_φ = 11.10 MeV,  α_D = 5.734×10⁻³
  Λ_d = 2.054 meV,  f = 0.240 M_Pl = 2.93×10¹⁸ GeV
"""

import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")

print("=" * 68)
print("TEST PI-3: STOCHASTIC INFLATION — P(θ_i)")
print("=" * 68)

# ─── Physical constants (same convention as PI-7) ────────────────────────────
M_Pl_GeV   = 1.221e19      # Planck mass [GeV]
hbar_GeV_s = 6.582e-25     # ℏ [GeV·s]

# Model parameters (from PI-7)
f_frac     = 0.240         # f / M_Pl
f_GeV      = f_frac * M_Pl_GeV   # 2.930e18 GeV
Lambda_d_eV  = 2.0535e-3   # dark confinement scale [eV]
Lambda_d_GeV = Lambda_d_eV * 1e-9
m_sigma_GeV  = Lambda_d_GeV**2 / f_GeV   # GMOR: m_σ = Λ_d²/f ≈ H₀
theta_A4   = np.arcsin(1/3)               # A₄ angle = 0.3398 rad
theta_i_target = 2.0                      # rad — required for H₀ = 67.4 (PI-7)

# Inflation parameters
N_efolds            = 60          # standard inflation e-folds
H_inf_standard      = 1e13        # GeV — typical GUT-scale inflation
H_inf_Planck_bound  = 6.0e13      # GeV — Planck 2018 upper bound (r < 0.036)
H_inf_BICEP_bound   = 2.5e13      # GeV — BICEP/Keck 2021

print(f"\n  Model parameters:")
print(f"    f              = {f_frac:.3f} × M_Pl = {f_GeV:.3e} GeV")
print(f"    Λ_d            = {Lambda_d_eV:.4f} meV = {Lambda_d_GeV:.3e} GeV")
print(f"    m_σ            = {m_sigma_GeV:.3e} GeV  (~H₀ ✓)")
print(f"    θ_A₄           = {theta_A4:.4f} rad  (A₄ geometry, frozen by PI-2)")
print(f"    θ_i,target     = {theta_i_target:.1f} rad  (needed for H₀ = 67.4 in PI-7)")
print(f"    N_e            = {N_efolds}")

# ─── Helper functions ─────────────────────────────────────────────────────────
def sigma_rms(H_inf, f, N):
    """RMS misalignment angle after N e-folds of de Sitter diffusion."""
    return (H_inf / (2.0 * np.pi * f)) * np.sqrt(N)

def P_gaussian(theta, sigma):
    """Gaussian P(θ) — valid for σ_rms ≪ π."""
    return np.exp(-0.5 * (theta / sigma)**2) / (np.sqrt(2.0 * np.pi) * sigma)

def log10_P_gaussian(theta, sigma):
    """log₁₀ P(θ) to avoid overflow."""
    return -0.5 * (theta / sigma)**2 / np.log(10) - np.log10(np.sqrt(2.0 * np.pi) * sigma)

def P_circle(theta_arr, D, N, n_modes=2000):
    """
    Exact P(θ, N) on the circle [0, 2π] starting from θ=0:
      P = 1/(2π) + (1/π) Σ_{n≥1} exp(−n²DN) cos(nθ)
    Converges to uniform for large DN.
    """
    result = np.ones_like(theta_arr, dtype=float) / (2.0 * np.pi)
    for n in range(1, n_modes + 1):
        coeff = np.exp(-n * n * D * N)
        if coeff < 1e-300:
            break
        result += coeff * np.cos(n * theta_arr) / np.pi
    return result

# ─── Section 1: Standard inflation ───────────────────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 1 — Standard Inflation (H_inf = 10¹³ GeV)")
print(f"{'─'*68}")

D_std     = (H_inf_standard / (2.0 * np.pi * f_GeV))**2
sig_std   = sigma_rms(H_inf_standard, f_GeV, N_efolds)
ratio_std = theta_i_target / sig_std        # how many σ away
log10_P_std = log10_P_gaussian(theta_i_target, sig_std)

print(f"  D = (H_inf/2πf)²     = {D_std:.3e}  [per e-fold]")
print(f"  σ_rms after {N_efolds} e-folds = {sig_std:.3e} rad")
print(f"  θ_i,target / σ_rms    = {ratio_std:.3e}  (how many σ away)")
print(f"  P(θ_i = 2 rad)       ~ 10^{log10_P_std:.2e}")
print(f"  Fine-tuning level    ~ 1 : 10^{abs(log10_P_std):.1e}")
print()

if sig_std < 0.01:
    print("  REGIME: Tiny diffusion (σ_rms ≪ 1 rad)")
    print("          θ_i is NOT set by de Sitter quantum diffusion at this scale.")
    print("          The misalignment angle is effectively the initial classical value.")
elif sig_std > np.pi:
    print("  REGIME: Ergodic (σ_rms ≫ π)")
    print("          θ_i is uniform over [0, π] — P(θ=2) ~ 0.32  NATURAL ✅")
else:
    print("  REGIME: Intermediate diffusion")

# ─── Section 2: Required H_inf for naturalness ───────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 2 — What H_inf gives natural θ_i = 2 rad?")
print(f"{'─'*68}")

# σ_rms = θ_i,target → H_inf* = θ_i × 2πf / √N_e
H_inf_1sigma   = theta_i_target * 2.0 * np.pi * f_GeV / np.sqrt(N_efolds)
# σ_rms = π → ergodic crossover
H_inf_ergodic  = np.pi * 2.0 * np.pi * f_GeV / np.sqrt(N_efolds)

print(f"  For σ_rms = θ_i = {theta_i_target:.1f} rad (1σ natural): H_inf* = {H_inf_1sigma:.3e} GeV")
print(f"  For ergodic crossover  (σ_rms = π):    H_inf  = {H_inf_ergodic:.3e} GeV")
print()
print(f"  Planck 2018 upper bound (r < 0.036):   H_inf < {H_inf_Planck_bound:.2e} GeV")
print(f"  BICEP/Keck 2021 bound  (r < 0.036):   H_inf < {H_inf_BICEP_bound:.2e} GeV")
print()

gap_1sigma  = H_inf_1sigma / H_inf_Planck_bound
gap_ergodic = H_inf_ergodic / H_inf_Planck_bound

print(f"  Gap to natural (1σ → Planck bound):   {H_inf_1sigma:.2e} / {H_inf_Planck_bound:.2e}")
print(f"                                        = {gap_1sigma:.1e}  ({np.log10(gap_1sigma):.1f} orders of magnitude)")
print()

if H_inf_1sigma > H_inf_Planck_bound:
    print("  VERDICT: θ_i = 2 rad is NOT naturally produced by de Sitter fluctuations")
    print(f"           at CMB-allowed inflation scales.")
    print(f"           The gap is {np.log10(gap_1sigma):.0f} orders of magnitude.")
    print()
    print("  Possible resolutions:")
    print("   (A) Accept θ_i as classical initial condition (misalignment scenario)")
    print("       H₀ = 67.4 is then a CONSISTENCY CHECK, not a parameter-free prediction")
    print("   (B) Thermal attractor: V_eff(θ,T) in early universe pulls θ → θ_i ~ 2")
    print("       → this is what PI-1 will test")
    print("   (C) Anthropic/landscape: patches with θ_i ~ 2 give right ρ_Λ → Boltzmann")
else:
    print("  VERDICT: ✅ θ_i = 2 rad naturally produced at observable inflation scales.")

# ─── Section 3: H_inf scan table ─────────────────────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 3 — σ_rms and P(θ=2) vs H_inf")
print(f"{'─'*68}")

header = f"  {'H_inf [GeV]':>14}  {'σ_rms [rad]':>12}  {'θ/σ':>12}  {'log₁₀ P(θ=2)':>16}  Regime"
print(header)
print(f"  {'─'*14}  {'─'*12}  {'─'*12}  {'─'*16}  {'─'*20}")

scan_cases = [
    1e10, 1e12, H_inf_BICEP_bound, H_inf_Planck_bound,
    1e14, 1e15, 1e16, H_inf_1sigma, H_inf_ergodic, 1e19
]

for H in scan_cases:
    s = sigma_rms(H, f_GeV, N_efolds)
    if s < np.pi:
        lp = log10_P_gaussian(theta_i_target, s)
        ratio = theta_i_target / s
        if s < 0.01:
            regime = "tiny diffusion"
        elif s < 0.5:
            regime = "small diffusion"
        else:
            regime = "intermediate"
    else:
        lp = np.log10(1.0 / np.pi)   # uniform → P(θ=2) ~ 1/π
        ratio = 0.0
        regime = "ergodic (uniform)"
        lp_str = f"{lp:>16.2f}"

    marker = " ← Planck bound" if abs(H - H_inf_Planck_bound) < 1 else (
             " ← 1σ natural"     if abs(H - H_inf_1sigma) < H_inf_1sigma * 0.01 else (
             " ← ergodic"        if abs(H - H_inf_ergodic) < H_inf_ergodic * 0.01 else ""))
    if s < np.pi:
        print(f"  {H:>14.2e}  {s:>12.3e}  {ratio:>12.2e}  {lp:>16.2e}  {regime}{marker}")
    else:
        print(f"  {H:>14.2e}  {s:>12.3e}  {'—':>12}  {lp:>16.4f}  {regime}{marker}")

# ─── Section 4: Fokker-Planck drift vs diffusion ─────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 4 — Fokker-Planck: drift vs diffusion at H_inf = 10¹³ GeV")
print(f"{'─'*68}")

# V(θ) = Λ_d⁴(1−cosθ),  V'(θ) ~ Λ_d⁴ sinθ → max at θ=π/2: V'_max = Λ_d⁴
# Drift coefficient μ = V'/(3H²_inf f²) ~ Λ_d⁴/(3 H²_inf f²)
# Diffusion coefficient D = (H_inf/2πf)²
Lambda_d4     = Lambda_d_GeV**4
drift_max     = Lambda_d4 / (3.0 * H_inf_standard**2 * f_GeV**2)   # μ_max
ratio_drift   = drift_max / D_std

print(f"  D  (diffusion)            = (H_inf/2πf)²          = {D_std:.3e}  per e-fold")
print(f"  μ  (drift, θ=π/2)        = Λ_d⁴/(3H²f²)         = {drift_max:.3e}  per e-fold")
print(f"  μ/D                       = Λ_d⁴/(3H⁴)×(2π)²    = {ratio_drift:.3e}")
print()
print(f"  Λ_d⁴ = {Lambda_d4:.3e} GeV⁴")
print(f"  3H_inf²f² = {3*H_inf_standard**2*f_GeV**2:.3e} GeV⁴")
print()

if ratio_drift < 1e-6:
    print(f"  DRIFT/DIFFUSION ratio ~ 10^{np.log10(ratio_drift):.0f}")
    print(f"  ✅ Pure diffusion — V(θ) potential is completely negligible during inflation.")
    print(f"     Including the dark QCD potential changes P(θ_i) by less than 1 part in 10^{abs(int(np.log10(ratio_drift)))}")
else:
    print(f"  Drift contributes at level {ratio_drift:.2e}")

# ─── Section 5: P(θ_i) plot ──────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    theta_arr = np.linspace(0, np.pi, 2000)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    cases = [
        (H_inf_BICEP_bound,  f"H = {H_inf_BICEP_bound:.1e} GeV (BICEP bound)", "blue",   "-"),
        (H_inf_Planck_bound, f"H = {H_inf_Planck_bound:.1e} GeV (Planck bound)", "royalblue", "--"),
        (1e15,               "H = 10¹⁵ GeV",            "green",  "-"),
        (H_inf_1sigma,       f"H = {H_inf_1sigma:.1e} GeV (1σ natural)", "orange", "-"),
        (H_inf_ergodic,      f"H = {H_inf_ergodic:.1e} GeV (ergodic)", "red",    "-"),
    ]

    ax1, ax2 = axes
    for H, label, color, ls in cases:
        s  = sigma_rms(H, f_GeV, N_efolds)
        D_ = (H / (2.0 * np.pi * f_GeV))**2
        if s < np.pi:
            P  = P_gaussian(theta_arr, s)
            P  = np.clip(P, 0.0, None)
        else:
            P  = P_circle(theta_arr, D_, N_efolds, n_modes=500)
            P  = np.clip(P, 0.0, None)
        norm = np.trapezoid(P, theta_arr) if hasattr(np, 'trapezoid') else np.trapz(P, theta_arr)
        if norm > 0:
            P /= norm
        ax1.semilogy(theta_arr, P + 1e-300, label=label, color=color, ls=ls)
        ax2.plot(theta_arr, P, label=label, color=color, ls=ls)

    for ax in [ax1, ax2]:
        ax.axvline(theta_i_target, color='k',      lw=1.5, ls='--', label=f'θ_i = {theta_i_target} rad (target)')
        ax.axvline(theta_A4,       color='purple',  lw=1.0, ls=':',  label=f'θ_A4 = {theta_A4:.2f} rad')
        ax.set_xlabel('θ_i  [rad]', fontsize=11)
        ax.set_xlim(0, np.pi)

    ax1.set_ylabel('P(θ_i)  [log scale]', fontsize=11)
    ax1.set_title('Stochastic Inflation: P(θ_i) — log scale', fontsize=11)
    ax1.set_ylim(1e-5, None)
    ax1.legend(fontsize=7.5, loc='upper right')

    ax2.set_ylabel('P(θ_i)  [linear]', fontsize=11)
    ax2.set_title('Stochastic Inflation: P(θ_i) — linear scale', fontsize=11)
    ax2.legend(fontsize=7.5, loc='upper right')

    plt.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'PI3_stochastic_inflation.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Plot saved: {out_path}")
except Exception as e:
    print(f"\n  (Plot skipped: {e})")

# ─── Section 6: Implication for H₀ prediction ────────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 5 — Implication for H₀ = 67.4 (PI-7 result)")
print(f"{'─'*68}")

print()
print("  PI-7 showed:  V_eff(σ₀) = ρ_Λ  →  H₀ = 67.4  IF  θ_i ~ 2 rad")
print("  PI-3 shows:   θ_i is NOT naturally produced by de Sitter diffusion")
print()
print("  Status of H₀ as a model prediction:")
print()
print("   STRONG prediction:  θ_i fixed by inflation → no free parameter")
print("   WEAK prediction:    θ_i given as initial condition → H₀ = consistency check")
print()
print("  Current status: WEAK — θ_i is a free parameter (misalignment angle)")
print()
print("  What would make it STRONG:")
print("   → PI-1: Does V_eff(θ,T) at freeze-out create an attractor toward θ_i ~ 2?")
print("           If V_eff(T_fo) has a minimum near θ = 2 rad → θ_i is dynamically set")
print("           → H₀ = 67.4 becomes a TRUE parameter-free prediction")

# ─── SUMMARY ─────────────────────────────────────────────────────────────────
print()
print("=" * 68)
print("PI-3 SUMMARY")
print("=" * 68)
print()
print(f"  f = {f_frac:.3f} M_Pl,  θ_i,target = {theta_i_target:.1f} rad,  N_e = {N_efolds}")
print()

rows = [
    ("H_inf [standard]",      f"{H_inf_standard:.0e} GeV"),
    ("σ_rms [standard]",      f"{sig_std:.2e} rad"),
    ("Fine-tuning level",     f"1 : 10^({abs(log10_P_std):.2e})"),
    ("H_inf for σ_rms=θ_i",  f"{H_inf_1sigma:.2e} GeV"),
    ("CMB upper bound",       f"{H_inf_Planck_bound:.1e} GeV (r<0.036)"),
    ("Gap (orders)",          f"{np.log10(gap_1sigma):.0f} orders of magnitude"),
    ("Drift/Diffusion ratio", f"{ratio_drift:.1e}  (pure diffusion ✅)"),
]

col1 = max(len(r[0]) for r in rows) + 2
for label, val in rows:
    print(f"  {label:<{col1}} {val}")

print()
print("  ┌──────────────────────────────────────────────────────────────────┐")
print("  │  VERDICT:                                                        │")
print("  │  ✗ θ_i = 2 rad is NOT naturally produced by stochastic          │")
print("  │    inflation at CMB-allowed energies.                            │")
print("  │    Gap: ~5 orders of magnitude below the needed H_inf.          │")
print("  │                                                                  │")
print("  │  θ_i is a FREE PARAMETER (classical misalignment angle).        │")
print("  │                                                                  │")
print("  │  H₀ = 67.4 (PI-7) is a CONSISTENCY CHECK given θ_i,            │")
print("  │  not a parameter-free prediction — unless PI-1 shows a          │")
print("  │  thermal attractor that dynamically fixes θ_i ~ 2 rad.         │")
print("  └──────────────────────────────────────────────────────────────────┘")
print()
print("  CHAIN STATUS:")
print(f"   ✅ PI-2:  θ_A₄ = arcsin(1/3) stable — S_E ~ 10¹²¹ (Hubble-radius bubble)")
print(f"   ✗  PI-3:  θ_i origin unclear — free parameter at standard H_inf")
print(f"   ✅ PI-7:  H₀ = 67.4 DERIVABLE from V_eff(σ₀) given θ_i")
print()
print(f"  NEXT: PI-1 — Does V_eff(θ,T) create a thermal attractor at θ ~ 2 rad?")
print(f"        If YES → H₀ = 67.4 becomes a genuine zero-free-parameter prediction.")
print()
print("=" * 68)
print("Test PI-3 COMPLETE")
print("=" * 68)
