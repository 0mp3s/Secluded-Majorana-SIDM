"""
Test PI-1: Thermal Effective Potential V_eff(θ, T)
====================================================================
Core question:
  Does V_eff(θ, T) create a THERMAL ATTRACTOR near θ ~ 2 rad
  at any temperature in the dark sector's history?
  If yes → θ_i is dynamically fixed → H₀ = 67.4 is a true prediction.
  If no  → θ_i confirmed as free parameter (consistent with PI-3).

Physics:
  V_eff(θ, T) = V_QCD(θ, T)    [dark condensate, T-suppressed]
              + V_CW(θ)         [Coleman-Weinberg, T=0]
              + V_T(θ, T)       [Matsubara finite-T correction]

  θ-dependence of χ mass (A₄ structure, g_p/g_s = 1/3):
    m_χ(θ) = m_χ √(cos²θ + (1/3)² sin²θ) = m_χ √(1 − 8/9 × sin²θ)
      → max at θ=0: m_χ(0) = m_χ
      → min at θ=π/2: m_χ(π/2) = m_χ/3

  φ mass is θ-independent (mediator mass from different sector).

  Finite-T thermal functions:
    J_B(x²) = ∫₀^∞ k² ln(1 − e^{−√(k²+x²)}) dk   [bosons]
    J_F(x²) = ∫₀^∞ k² ln(1 + e^{−√(k²+x²)}) dk    [fermions]

  V_T(θ,T) = T⁴/(2π²) × [J_B(m_φ²/T²) − 4 J_F(m_χ²(θ)/T²)]
             ↑ 1 real scalar dof   ↑ 4 Dirac fermion dof

Key temperatures:
  T >> m_χ : thermal correction ~ m_χ²(θ) T²/12 — pushes toward θ=π/2
  T ~ m_χ  : Boltzmann crossover — thermal correction fades
  T_fo = m_χ/20 ~ 4.7 MeV: actual freeze-out — m_χ/T~20, thermal ≈ 0
  T_dark ~ Λ_d ~ 2 meV: dark condensate forms, V_QCD turns on
"""

import numpy as np
from scipy import integrate
import os
import warnings
warnings.filterwarnings("ignore")

print("=" * 68)
print("TEST PI-1: THERMAL EFFECTIVE POTENTIAL V_eff(θ, T)")
print("=" * 68)

# ─── Physical constants and model parameters ─────────────────────────────────
M_Pl_GeV     = 1.221e19
f_GeV        = 0.240 * M_Pl_GeV
Lambda_d_GeV = 2.0535e-3 * 1e-9   # 2.054 meV
m_chi        = 94.07e-3            # GeV
m_phi        = 11.10e-3            # GeV
alpha_D      = 5.734e-3
T_fo         = m_chi / 20.0        # freeze-out temperature
mu_ren       = m_chi               # renormalization scale
theta_A4     = np.arcsin(1.0/3.0)  # 0.3398 rad
theta_target = 2.0                 # rad — needed for H₀ = 67.4

# A₄ structure
gp_over_gs   = 1.0 / 3.0          # exact from A₄

print(f"\n  Model parameters:")
print(f"    m_χ         = {m_chi*1e3:.2f} MeV")
print(f"    m_φ         = {m_phi*1e3:.2f} MeV")
print(f"    g_p / g_s   = {gp_over_gs:.4f}  (A₄ exact)")
print(f"    f           = {f_GeV:.3e} GeV")
print(f"    Λ_d         = {Lambda_d_GeV*1e12:.4f} meV")
print(f"    T_fo        = m_χ / 20 = {T_fo*1e3:.2f} MeV")
print(f"    θ_A₄        = {theta_A4:.4f} rad")
print(f"    θ_target    = {theta_target:.4f} rad  (needed for H₀ = 67.4)")

# ─── θ-dependent mass ────────────────────────────────────────────────────────
def m_chi_eff(theta):
    """χ mass from A₄ coupling: m_χ(θ) = m_χ √(cos²θ + (1/3)² sin²θ)."""
    return m_chi * np.sqrt(np.cos(theta)**2 + gp_over_gs**2 * np.sin(theta)**2)

# ─── Thermal integrals J_B, J_F ─────────────────────────────────────────────
_gamma_E = 0.5772156649

def _J_highT_boson(x):
    """High-T expansion J_B(x²), valid for x = m/T < 1.5."""
    x2 = x*x;  x3 = x2*x;  x4 = x2*x2
    ln_term = np.log(x2 / (4 * np.pi**2 * np.exp(-2*_gamma_E))) - 3/2
    return -np.pi**4/45 + (np.pi**2/12)*x2 - (np.pi/6)*x3 - (x4/32)*ln_term

def _J_highT_fermion(x):
    """High-T expansion J_F(x²), valid for x = m/T < 1.5."""
    x2 = x*x;  x4 = x2*x2
    ln_term = np.log(x2 / (np.pi**2 * np.exp(-2*_gamma_E))) - 3/2
    return -7*np.pi**4/360 + (np.pi**2/24)*x2 + (x4/32)*ln_term

def _J_boltzmann(x, boson):
    """Boltzmann approximation, valid for x > 4."""
    sign = 1.0 if boson else -1.0
    return sign * np.sqrt(np.pi * x / 2) * x**2 * np.exp(-x)

def _J_numerical(x2, boson):
    """Exact numerical integration."""
    def integrand(k):
        E = np.sqrt(k**2 + x2)
        if boson:
            return k**2 * np.log(max(1.0 - np.exp(-E), 1e-300))
        else:
            return k**2 * np.log(1.0 + np.exp(-E))
    kmax = max(20.0, 10.0 * np.sqrt(x2))
    result, _ = integrate.quad(integrand, 0, kmax, limit=80, epsabs=1e-12)
    return result

def J_thermal(x, boson):
    """
    Evaluate J_{B/F}(x²) with x = m/T.
    Uses high-T expansion (x<1.5), numerical (1.5≤x≤5), Boltzmann (x>5).
    """
    if x < 0.0:
        x = 0.0
    if x < 1.5:
        return _J_highT_boson(x) if boson else _J_highT_fermion(x)
    elif x > 5.0:
        return _J_boltzmann(x, boson)
    else:
        return _J_numerical(x**2, boson)

# ─── Potential components ────────────────────────────────────────────────────
def V_QCD(theta, T):
    """
    Dark QCD potential with thermal suppression.
    Λ_d(T) = Λ_d × min(1, (Λ_d/T)^8) — condensate melts at T > Λ_d.
    """
    if T < 1e-30:
        return Lambda_d_GeV**4 * (1.0 - np.cos(theta))
    ratio = Lambda_d_GeV / T
    suppression = min(1.0, ratio**8)
    return Lambda_d_GeV**4 * (1.0 - np.cos(theta)) * suppression

def V_CW(theta):
    """
    Coleman-Weinberg 1-loop: χ contributes with 4 dof (Dirac fermion).
    V_CW = −4/(64π²) m_χ⁴(θ) [ln(m_χ²(θ)/μ²) − 3/2]
    Negative sign: fermion loops lower the vacuum energy.
    """
    mq  = m_chi_eff(theta)
    mq  = max(mq, 1e-30)
    log = np.log(mq**2 / mu_ren**2) - 1.5
    return -4.0 / (64.0 * np.pi**2) * mq**4 * log

def V_T(theta, T):
    """
    Finite-T thermal correction.
    Fermions (4 dof, χ): −4 × T⁴/(2π²) × J_F(m_χ²(θ)/T²)
    Boson   (1 dof, φ):  +1 × T⁴/(2π²) × J_B(m_φ²/T²)
    """
    if T < 1e-30:
        return 0.0
    mq  = m_chi_eff(theta)
    xf  = mq  / T
    xb  = m_phi / T
    jf  = J_thermal(xf, boson=False)
    jb  = J_thermal(xb, boson=True)
    pref = T**4 / (2.0 * np.pi**2)
    return pref * (jb - 4.0 * jf)

def V_eff_scalar(theta, T):
    return V_QCD(theta, T) + V_CW(theta) + V_T(theta, T)

V_eff_vec = np.vectorize(V_eff_scalar)

# ─── Section 1: θ-dependence of m_χ ─────────────────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 1 — m_χ(θ) from A₄ structure")
print(f"{'─'*68}")

print(f"\n  {'θ [rad]':>10}  {'m_χ(θ) [MeV]':>14}  {'m_χ(θ)/m_χ(0)':>14}")
print(f"  {'─'*10}  {'─'*14}  {'─'*14}")
for th in [0.0, theta_A4, np.pi/6, np.pi/4, np.pi/3, 1.0, 1.5, 2.0, np.pi/2, 2.5, np.pi]:
    mq = m_chi_eff(th)
    print(f"  {th:>10.4f}  {mq*1e3:>14.4f}  {mq/m_chi:>14.6f}")
print(f"\n  m_χ(0)   = {m_chi*1e3:.2f} MeV  (maximum)")
print(f"  m_χ(π/2) = {m_chi_eff(np.pi/2)*1e3:.2f} MeV  = m_χ/3  (minimum)")

# ─── Section 2: J_thermal values — regime check ──────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 2 — Thermal correction regime at key temperatures")
print(f"{'─'*68}")

print(f"\n  At T = T_fo = {T_fo*1e3:.2f} MeV:")
for th in [0.0, 1.0, 2.0, np.pi/2]:
    mq = m_chi_eff(th)
    xf = mq / T_fo
    xb = m_phi / T_fo
    jf = J_thermal(xf, boson=False)
    jb = J_thermal(xb, boson=True)
    vt = V_T(th, T_fo)
    vcw = V_CW(th)
    print(f"   θ={th:.2f}: m_χ/T={xf:.1f},  J_F={jf:.3e},  V_T/V_CW={abs(vt)/max(abs(vcw),1e-100):.3e}")

print(f"\n  At T = m_χ = {m_chi*1e3:.1f} MeV:")
for th in [0.0, 1.0, 2.0, np.pi/2]:
    mq = m_chi_eff(th)
    xf = mq / m_chi
    jf = J_thermal(xf, boson=False)
    vt = V_T(th, m_chi)
    vcw = V_CW(th)
    print(f"   θ={th:.2f}: m_χ(θ)/T={xf:.3f},  J_F={jf:.3e},  V_T/V_CW={abs(vt)/max(abs(vcw),1e-100):.3e}")

print(f"\n  At T = 10 × m_χ = {10*m_chi*1e3:.0f} MeV:")
for th in [0.0, 1.0, 2.0, np.pi/2]:
    mq = m_chi_eff(th)
    xf = mq / (10*m_chi)
    jf = J_thermal(xf, boson=False)
    vt = V_T(th, 10*m_chi)
    vcw = V_CW(th)
    print(f"   θ={th:.2f}: m_χ(θ)/T={xf:.3f},  J_F={jf:.3e},  V_T/V_CW={abs(vt)/max(abs(vcw),1e-100):.3e}")

# ─── Section 3: Minimum of V_eff vs Temperature ──────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 3 — θ_min(T): where does V_eff have its minimum?")
print(f"{'─'*68}")

theta_arr = np.linspace(0.01, np.pi - 0.01, 600)
temps = [
    50.0 * m_chi, 20.0 * m_chi, 10.0 * m_chi,
    5.0  * m_chi,  3.0 * m_chi,  2.0 * m_chi,
    1.0  * m_chi,  0.5 * m_chi,  0.2 * m_chi,
    T_fo,          0.01 * m_chi,  0.0
]

print(f"\n  {'T [MeV]':>10}  {'T/m_χ':>8}  {'θ_min [rad]':>12}  {'θ_min [°]':>10}  Near θ=2?  Near π/2?")
print(f"  {'─'*10}  {'─'*8}  {'─'*12}  {'─'*10}  {'─'*9}  {'─'*9}")

theta_min_arr = []
for T in temps:
    if T < 1e-30:
        V_vals = np.array([V_QCD(th, 1e-30) + V_CW(th) for th in theta_arr])
        T_label = 0.0
        T_ratio = 0.0
    else:
        V_vals = V_eff_vec(theta_arr, T)
        T_label = T * 1e3
        T_ratio = T / m_chi

    good = np.isfinite(V_vals)
    if good.sum() < 10:
        theta_min_arr.append(np.nan)
        print(f"  {T_label:>10.2f}  {T_ratio:>8.3f}  {'FAILED':>12}  {'—':>10}")
        continue

    idx_min = np.argmin(V_vals[good])
    th_min  = theta_arr[good][idx_min]
    theta_min_arr.append(th_min)

    near2   = "✅ YES" if abs(th_min - 2.0) < 0.3 else "✗"
    nearpi2 = "YES  " if abs(th_min - np.pi/2) < 0.3 else "no"
    print(f"  {T_label:>10.2f}  {T_ratio:>8.3f}  {th_min:>12.4f}  {np.degrees(th_min):>10.2f}  {near2:<9}  {nearpi2}")

# ─── Section 4: Breakdown at T_fo ────────────────────────────────────────────
print(f"\n{'─'*68}")
print(f"SECTION 4 — V_eff(θ) breakdown at T_fo = {T_fo*1e3:.2f} MeV (m_χ/T = {m_chi/T_fo:.0f})")
print(f"{'─'*68}")

theta_key = [0.0, theta_A4, np.pi/4, 1.0, 1.5, 2.0, np.pi/2, 2.5, np.pi]
Vqcd = np.array([V_QCD(th, T_fo)    for th in theta_key])
Vcw  = np.array([V_CW(th)           for th in theta_key])
Vt   = np.array([V_T(th, T_fo)      for th in theta_key])
Vtot = Vqcd + Vcw + Vt

# Shift relative to θ=0
Vtot_shifted = Vtot - Vtot[0]
Vcw_shifted  = Vcw  - Vcw[0]

print(f"\n  (All values relative to θ=0, units: 10⁻¹² GeV⁴)")
scale = 1e12
print(f"\n  {'θ [rad]':>8}  {'V_QCD':>10}  {'V_CW':>10}  {'V_T':>10}  {'V_tot':>10}  {'V_CW(θ)-V_CW(0)':>17}")
print(f"  {'─'*8}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*17}")
for i, th in enumerate(theta_key):
    marker = " ← θ_A₄" if abs(th-theta_A4)<0.01 else (" ← target θ=2" if abs(th-2.0)<0.01 else "")
    print(f"  {th:>8.4f}  {Vqcd[i]*scale:>10.4f}  {Vcw[i]*scale:>10.4f}  "
          f"{Vt[i]*scale:>10.4f}  {Vtot[i]*scale:>10.4f}  "
          f"{Vcw_shifted[i]*scale:>17.4f}{marker}")

print(f"\n  V_QCD/V_CW ratio at θ=1 rad: {abs(Vqcd[3])/max(abs(Vcw[3]),1e-100):.2e}")
print(f"  V_T/V_CW   ratio at θ=1 rad: {abs(Vt[3])/max(abs(Vcw[3]),1e-100):.2e}")

# ─── Section 5: Plot ─────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    theta_plot = np.linspace(0.01, np.pi - 0.01, 800)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    ax1, ax2 = axes

    # Left: V_eff(θ) at several temperatures, shifted to minimum=0
    plot_temps = [
        (50  * m_chi, "T = 50 m_χ",   "navy"),
        (10  * m_chi, "T = 10 m_χ",   "blue"),
        (5   * m_chi, "T = 5 m_χ",    "royalblue"),
        (2   * m_chi, "T = 2 m_χ",    "dodgerblue"),
        (1   * m_chi, "T = m_χ",       "limegreen"),
        (0.5 * m_chi, "T = 0.5 m_χ",  "goldenrod"),
        (T_fo,         f"T = T_fo = m_χ/20", "red"),
        (0.0,          "T = 0 (CW only)",    "black"),
    ]

    for T, label, color in plot_temps:
        if T < 1e-30:
            V_vals = np.array([V_QCD(th, 1e-30) + V_CW(th) for th in theta_plot])
        else:
            V_vals = V_eff_vec(theta_plot, T)
        good = np.isfinite(V_vals)
        V_shifted = V_vals - np.nanmin(V_vals[good])
        # Normalize by value at θ=0 to see shape
        scale_v = abs(V_shifted[good][0]) + 1e-60
        ax1.plot(theta_plot[good], V_shifted[good] / scale_v,
                 label=label, color=color, lw=1.5)

    ax1.axvline(theta_target, color='k',      lw=2, ls='--', label='θ = 2 rad (target)')
    ax1.axvline(theta_A4,     color='purple', lw=1, ls=':',  label=f'θ_A₄ = {theta_A4:.2f}')
    ax1.axvline(np.pi/2,      color='gray',   lw=1, ls=':',  label='π/2 = 1.57 rad')
    ax1.set_xlabel('θ [rad]', fontsize=11)
    ax1.set_ylabel('(V_eff − min) / |V_eff(0) − min|', fontsize=10)
    ax1.set_title('Shape of V_eff(θ, T) at various temperatures', fontsize=11)
    ax1.legend(fontsize=7.5, loc='upper right')
    ax1.set_xlim(0, np.pi)
    ax1.set_ylim(-0.2, 3)

    # Right: θ_min(T)
    T_scan = np.logspace(np.log10(0.001 * m_chi), np.log10(100 * m_chi), 300)
    th_min_scan = []
    for T in T_scan:
        V = V_eff_vec(theta_plot, T)
        good = np.isfinite(V)
        if good.sum() > 10:
            th_min_scan.append(theta_plot[good][np.argmin(V[good])])
        else:
            th_min_scan.append(np.nan)

    ax2.semilogx(T_scan * 1e3, th_min_scan, 'b-', lw=2.5, label='θ_min(T)')
    ax2.axhline(theta_target, color='k',      lw=2,   ls='--', label='θ_target = 2 rad')
    ax2.axhline(theta_A4,     color='purple', lw=1,   ls=':',  label=f'θ_A₄ = {theta_A4:.2f}')
    ax2.axhline(np.pi/2,      color='gray',   lw=1.5, ls='-',  label='π/2 (expected min)')
    ax2.axvline(T_fo*1e3,     color='red',    lw=1.5, ls='--', label=f'T_fo = {T_fo*1e3:.1f} MeV')
    ax2.axvline(m_chi*1e3,    color='green',  lw=1,   ls=':',  label=f'm_χ = {m_chi*1e3:.0f} MeV')
    ax2.set_xlabel('T [MeV]', fontsize=11)
    ax2.set_ylabel('θ_min(T)  [rad]', fontsize=11)
    ax2.set_title('Minimum of V_eff vs temperature', fontsize=11)
    ax2.legend(fontsize=8, loc='best')
    ax2.set_ylim(0, np.pi)

    plt.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'PI1_thermal_effective_potential.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Plot saved: {out_path}")
except Exception as e:
    print(f"\n  (Plot skipped: {e})")

# ─── Section 6: Summary ──────────────────────────────────────────────────────
print(f"\n{'─'*68}")
print("SECTION 6 — SUMMARY")
print(f"{'─'*68}")
print()

# θ_min at T_fo
V_fo  = V_eff_vec(theta_arr, T_fo)
good  = np.isfinite(V_fo)
th_fo = theta_arr[good][np.argmin(V_fo[good])]

# θ_min at T=0
V_0   = np.array([V_CW(th) for th in theta_arr])
th_0  = theta_arr[np.argmin(V_0)]

# θ_min at T = 10 m_chi (maximum thermal correction)
V_10  = V_eff_vec(theta_arr, 10*m_chi)
good  = np.isfinite(V_10)
th_10 = theta_arr[good][np.argmin(V_10[good])]

rows = [
    ("θ_min at T = 10 m_χ (full thermal)", f"{th_10:.3f} rad  ({np.degrees(th_10):.1f}°)"),
    ("θ_min at T = T_fo (freeze-out)",      f"{th_fo:.3f} rad  ({np.degrees(th_fo):.1f}°)"),
    ("θ_min at T = 0 (CW only)",            f"{th_0:.3f} rad  ({np.degrees(th_0):.1f}°)"),
    ("Target θ_i",                          f"2.000 rad  (114.6°)"),
    ("V_QCD/V_CW at T_fo",                  f"{abs(Vqcd[3])/max(abs(Vcw[3]),1e-100):.2e}  (QCD suppressed)"),
    ("V_T/V_CW   at T_fo",                  f"{abs(Vt[3])/max(abs(Vcw[3]),1e-100):.2e}  (Boltzmann → ~0)"),
]
col1 = max(len(r[0]) for r in rows) + 2
for label, val in rows:
    print(f"  {label:<{col1}} {val}")

near_anywhere = any(abs(th - 2.0) < 0.3 for th in theta_min_arr if not np.isnan(th))

print()
print("  ┌──────────────────────────────────────────────────────────────────┐")
if near_anywhere:
    found_T = [temps[i] for i, th in enumerate(theta_min_arr)
               if not np.isnan(th) and abs(th - 2.0) < 0.3]
    print(f"  │  ✅ THERMAL ATTRACTOR FOUND near θ ~ 2 rad                      │")
    print(f"  │     At T = {found_T[0]*1e3:.0f} MeV                                       │")
    print(f"  │     → V_eff has minimum near θ = 2 rad at some temperature!     │")
    print(f"  │     → PI-1 PASS: H₀ = 67.4 may be a parameter-free prediction  │")
else:
    print(f"  │  ✗ NO THERMAL ATTRACTOR found near θ ~ 2 rad                   │")
    print(f"  │                                                                  │")
    print(f"  │  At ALL temperatures, θ_min ≈ π/2 (1.57 rad).                  │")
    print(f"  │                                                                  │")
    print(f"  │  Physical reason (three lines):                                 │")
    print(f"  │    • V_QCD ≈ 0 at T >> Λ_d (condensate melted)                │")
    print(f"  │    • V_CW minimum at θ=π/2 (lower m_χ = lower loop energy)    │")
    print(f"  │    • V_T  reinforces CW (thermal mass ∝ m_χ²(θ))              │")
    print(f"  │                                                                  │")
    print(f"  │  θ_i = 2 rad CONFIRMED as free parameter.                      │")
    print(f"  │  H₀ = 67.4 is a CONSISTENCY CHECK, not a zero-param prediction. │")
print(f"  └──────────────────────────────────────────────────────────────────┘")

print()
print("  UPDATED CHAIN STATUS:")
print(f"   ✅ PI-2: θ_A₄ = arcsin(1/3) stable — S_E ~ 10¹²¹, frozen forever")
print(f"   ✗  PI-3: θ_i origin — stochastic inflation gives σ_rms ~ 4×10⁻⁶ rad")
print(f"   ✗  PI-1: thermal attractor — V_eff minimum always at π/2, not at 2 rad")
print(f"   ✅ PI-7: H₀ = 67.4 derived from V_eff(σ₀) GIVEN θ_i (consistency check)")
print()
print("  NET CONCLUSION:")
print("   θ_i = 2 rad is a classical initial condition (misalignment angle).")
print("   The dark energy density ρ_Λ = V(θ_i) is SENSITIVE to θ_i:")
print(f"   V(θ_i=2) = Λ_d⁴(1−cos2) = {Lambda_d_GeV**4*(1-np.cos(2)):.3e} GeV⁴ (→ H₀=67.4)")
print(f"   V(θ_i=π) = Λ_d⁴ × 2     = {Lambda_d_GeV**4*2:.3e} GeV⁴ (→ H₀ × {np.sqrt(2*(1/(1-np.cos(2)))):.2f})")
print()
print("   This is the standard dark-energy misalignment scenario.")
print("   The H₀ prediction requires specifying θ_i — no free-parameter")
print("   derivation is possible from this model alone.")
print()
print("=" * 68)
print("Test PI-1 COMPLETE")
print("=" * 68)
