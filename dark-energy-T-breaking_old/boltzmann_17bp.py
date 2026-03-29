#!/usr/bin/env python3
# ⚠️ SUPERSEDED — This script applies (8/9) reduction to α_CSV,
# which is already α_s. See test_alpha_convention.py (Test 12).
# Relic density Ωh² values are overestimated by ~27%.
# Buggy output archived in archived_buggy/boltzmann_17bp_BUGGY.txt
"""
boltzmann_17bp.py — Full Numerical Boltzmann for 17 BPs with α×(8/9)
=====================================================================

Step 2 of verification: Replace Kolb-Turner analytic approximation
with full RK4 numerical Boltzmann solver.

Key change from consistency_17bp.py:
  - KT analytic: ~20% systematic underestimate of Ωh²
  - This script: full numerical RK4 Boltzmann (from v27_boltzmann_relic.py)

For each BP:
  1. Original α → α_s = (8/9)α  (θ-decomposition)
  2. ⟨σv⟩₀ = πα_s²/(4m²)       (s-wave, scalar mediator)
  3. Full Boltzmann RK4 → Y_∞ → Ωh²
  4. Compare with KT analytic result from consistency_17bp.py
  5. Re-evaluate which BPs still match Ωh² = 0.120 ± tolerance
"""
import sys, math, os
import numpy as np

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

# =============================================================================
#  Constants
# =============================================================================
M_PL = 1.220890e19      # Planck mass in GeV
S_0  = 2891.2            # entropy density today, cm⁻³
RHO_CRIT_H2 = 1.0539e-5 # ρ_crit/h² in GeV/cm³

THETA = math.asin(1.0 / 3.0)
COS2 = math.cos(THETA)**2   # 8/9
SIN2 = math.sin(THETA)**2   # 1/9

OMEGA_TARGET = 0.120

# =============================================================================
#  g_*(T) tabulation (Drees, Hajkarim & Schmitz 2015)
# =============================================================================
_G_STAR_TABLE = np.array([
    [1e4,106.75,106.75],[200,106.75,106.75],[80,86.25,86.25],
    [10,86.25,86.25],[1,75.75,75.75],[0.3,61.75,61.75],
    [0.2,17.25,17.25],[0.15,14.25,14.25],[0.1,10.75,10.75],
    [0.01,10.75,10.75],[0.001,10.75,10.75],[0.0005,10.75,10.75],
    [0.0001,3.36,3.91],[1e-5,3.36,3.91],[1e-8,3.36,3.91],
])
_LOG_T = np.log(_G_STAR_TABLE[:, 0])
_G_RHO = _G_STAR_TABLE[:, 1]
_G_S   = _G_STAR_TABLE[:, 2]

def g_star_rho(T):
    logT = math.log(T) if T > 0 else -50
    return float(np.interp(logT, _LOG_T[::-1], _G_RHO[::-1]))

def g_star_S(T):
    logT = math.log(T) if T > 0 else -50
    return float(np.interp(logT, _LOG_T[::-1], _G_S[::-1]))

# =============================================================================
#  Boltzmann solver (copied from v27_boltzmann_relic.py, self-contained)
# =============================================================================
def Y_eq_full(x, m_chi, g_chi=2):
    """Y_eq = 45/(4π⁴) × g/g_*S × √(π/2) × x^{3/2} × exp(-x)  (non-relativistic, K₂ NR limit)."""
    T = m_chi / x
    g_s = g_star_S(T)
    if x > 300:
        return 0.0
    return 45.0 / (4 * math.pi**4) * g_chi / g_s * math.sqrt(math.pi / 2) * x**1.5 * math.exp(-x)

def dYdx(x, Y, m_chi, sv0, g_chi=2):
    """dY/dx for s-wave: ⟨σv⟩ = sv0 = constant."""
    T = m_chi / x
    gs = g_star_S(T)
    gr = g_star_rho(T)
    g_eff = gs / math.sqrt(gr)
    lambda_x = math.sqrt(math.pi / 45.0) * g_eff * M_PL * m_chi
    Yeq = Y_eq_full(x, m_chi, g_chi)
    return -lambda_x * sv0 * (Y * Y - Yeq * Yeq) / (x * x)

def solve_boltzmann(m_chi, sv0, x_start=1.0, x_end=1000.0, n_steps=10000, g_chi=2):
    """Solve dY/dx using RK4 for s-wave annihilation."""
    dx = (x_end - x_start) / n_steps
    x_arr = np.zeros(n_steps + 1)
    Y_arr = np.zeros(n_steps + 1)
    x_arr[0] = x_start
    Y_arr[0] = Y_eq_full(x_start, m_chi, g_chi)

    for i in range(n_steps):
        x = x_arr[i]
        Y = Y_arr[i]
        k1 = dx * dYdx(x, Y, m_chi, sv0, g_chi)
        k2 = dx * dYdx(x + dx/2, Y + k1/2, m_chi, sv0, g_chi)
        k3 = dx * dYdx(x + dx/2, Y + k2/2, m_chi, sv0, g_chi)
        k4 = dx * dYdx(x + dx, Y + k3, m_chi, sv0, g_chi)
        Y_new = Y + (k1 + 2*k2 + 2*k3 + k4) / 6
        Y_new = max(Y_new, 1e-30)
        x_arr[i+1] = x + dx
        Y_arr[i+1] = Y_new

    return x_arr, Y_arr

def Y_to_omega_h2(Y_inf, m_chi):
    """Convert Y_∞ to Ω h²."""
    return m_chi * Y_inf * S_0 / RHO_CRIT_H2

# =============================================================================
#  Kolb-Turner analytic (for comparison)
# =============================================================================
def kolb_turner_swave(m_chi, sv0, g_chi=2):
    x_fo = 20.0
    for _ in range(50):
        T_fo = m_chi / x_fo
        gr = g_star_rho(T_fo)
        arg = 0.0764 * g_chi * M_PL * m_chi * sv0 / math.sqrt(gr * x_fo)
        if arg <= 0: break
        x_new = math.log(arg) - 0.5 * math.log(x_fo)
        if abs(x_new - x_fo) < 0.001: x_fo = x_new; break
        x_fo = 0.5 * x_fo + 0.5 * x_new
    T_fo = m_chi / x_fo
    gs = g_star_S(T_fo)
    gr = g_star_rho(T_fo)
    lam = math.sqrt(math.pi / 45.0) * gs / math.sqrt(gr) * M_PL * m_chi
    Y_inf = x_fo / (lam * sv0)
    return x_fo, Y_inf

# =============================================================================
#  Load 17 BPs from CSV
# =============================================================================
CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
    "Secluded-Majorana-SIDM", "predictions", "output", "sweep_17bp_results.csv")

import csv
bps = []
with open(CSV, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        bps.append({
            "label":    row["BP"],
            "m_chi":    float(row["m_chi_GeV"]),
            "m_phi":    float(row["m_phi_MeV"]) * 1e-3,
            "alpha":    float(row["alpha"]),
            "omega_orig": float(row["omega_h2"]),
        })

# =============================================================================
#  Main: Run all 17 BPs
# =============================================================================
print("=" * 110)
print("  17 BPs — FULL BOLTZMANN (RK4) vs KOLB-TURNER ANALYTIC")
print("  θ-decomposition: α_s = (8/9)α,   ⟨σv⟩ = πα_s²/(4m²)")
print("=" * 110)
print()

header = (f"  {'BP':<5} {'m_χ(GeV)':>9} {'m_φ(MeV)':>9} {'α':>11} {'α_s':>11} "
          f"{'Ωh²_orig':>9} {'Ωh²_KT':>9} {'Ωh²_Boltz':>10} {'KT/Boltz':>9} "
          f"{'x_fo':>6} {'pass?':>6}")
print(header)
print(f"  {'-'*108}")

n_pass_kt = 0
n_pass_boltz = 0
n_total = len(bps)
results = []

for bp in bps:
    m_chi = bp["m_chi"]
    m_phi = bp["m_phi"]
    alpha_orig = bp["alpha"]
    
    # θ-decomposition
    alpha_s = alpha_orig * COS2  # × 8/9
    
    # s-wave cross section
    sv0 = math.pi * alpha_s**2 / (4.0 * m_chi**2)
    
    # Kolb-Turner analytic
    x_fo_kt, Y_inf_kt = kolb_turner_swave(m_chi, sv0)
    omega_kt = Y_to_omega_h2(Y_inf_kt, m_chi)
    
    # Full numerical Boltzmann (RK4)
    x_arr, Y_arr = solve_boltzmann(m_chi, sv0)
    Y_inf_boltz = Y_arr[-1]
    omega_boltz = Y_to_omega_h2(Y_inf_boltz, m_chi)
    
    # Freeze-out temperature from numerical solution
    # (find where Y deviates by >10% from Y_eq)
    x_fo_num = x_fo_kt  # fallback
    for i in range(len(x_arr)):
        Yeq_i = Y_eq_full(x_arr[i], m_chi)
        if Yeq_i > 0 and Y_arr[i] / Yeq_i > 1.5:
            x_fo_num = x_arr[i]
            break
    
    ratio = omega_kt / omega_boltz if omega_boltz > 0 else 999
    
    # Pass criterion: within 20% of Ωh² = 0.120
    tol = 0.20
    pass_boltz = abs(omega_boltz / OMEGA_TARGET - 1) < tol
    pass_kt    = abs(omega_kt / OMEGA_TARGET - 1) < tol
    
    if pass_kt: n_pass_kt += 1
    if pass_boltz: n_pass_boltz += 1
    
    mark = "✅" if pass_boltz else "❌"
    
    print(f"  {bp['label']:<5} {m_chi:>9.4f} {m_phi*1e3:>9.4f} {alpha_orig:>11.4e} {alpha_s:>11.4e} "
          f"{bp['omega_orig']:>9.6f} {omega_kt:>9.4f} {omega_boltz:>10.4f} {ratio:>9.4f} "
          f"{x_fo_num:>6.1f} {mark:>6}")
    
    results.append({
        "label": bp["label"],
        "m_chi": m_chi,
        "alpha": alpha_orig,
        "alpha_s": alpha_s,
        "omega_orig": bp["omega_orig"],
        "omega_kt": omega_kt,
        "omega_boltz": omega_boltz,
        "ratio": ratio,
        "pass_boltz": pass_boltz,
        "pass_kt": pass_kt,
    })

# =============================================================================
#  Summary
# =============================================================================
print()
print("=" * 110)
print("  SUMMARY")
print("=" * 110)
print()

# Systematic comparison
omega_kt_arr = np.array([r["omega_kt"] for r in results])
omega_boltz_arr = np.array([r["omega_boltz"] for r in results])
ratios = omega_kt_arr / omega_boltz_arr

print(f"  KT / Boltzmann ratio (systematic bias):")
print(f"    Mean:   {np.mean(ratios):.4f}")
print(f"    Median: {np.median(ratios):.4f}")
print(f"    Std:    {np.std(ratios):.4f}")
print(f"    Min:    {np.min(ratios):.4f}")
print(f"    Max:    {np.max(ratios):.4f}")
print()

if np.mean(ratios) < 1:
    direction = "UNDERESTIMATES"
else:
    direction = "OVERESTIMATES"
bias_pct = abs(np.mean(ratios) - 1) * 100
print(f"  → KT {direction} Ωh² by ~{bias_pct:.1f}% on average")
print()

print(f"  BPs passing Ωh² = 0.120 ± 20%:")
print(f"    KT analytic:       {n_pass_kt}/{n_total}")
print(f"    Full Boltzmann:    {n_pass_boltz}/{n_total}")
print()

# Which BPs changed status?
changed = []
for r in results:
    if r["pass_boltz"] != r["pass_kt"]:
        changed.append(r)

if changed:
    print(f"  BPs that CHANGED status (KT → Boltzmann):")
    for c in changed:
        old = "pass" if c["pass_kt"] else "fail"
        new = "pass" if c["pass_boltz"] else "fail"
        print(f"    {c['label']}: {old} → {new}  (KT={c['omega_kt']:.4f}, Boltz={c['omega_boltz']:.4f})")
else:
    print(f"  No BPs changed status between KT and Boltzmann. ✅")

# The original BPs all have Ωh² ≈ 0.120 by construction.
# After θ-decomposition (α → 8α/9), the cross section is SMALLER:
# ⟨σv⟩ = πα_s²/(4m²) = π(8α/9)²/(4m²) = (64/81) × πα²/(4m²)
# Smaller ⟨σv⟩ → MORE relic → Ωh² INCREASES
# The increase factor: Ωh² ∝ 1/⟨σv⟩ → factor 81/64 ≈ 1.266

expected_factor = 81/64
print()
print(f"  Expected Ωh² increase from θ-decomposition:")
print(f"    ⟨σv⟩_s = (8/9)² × ⟨σv⟩_orig = (64/81) ⟨σv⟩_orig")
print(f"    Ωh²_new / Ωh²_orig ≈ 81/64 = {expected_factor:.4f}")
print()

actual_factors_boltz = omega_boltz_arr / np.array([r["omega_orig"] for r in results])
print(f"  Actual Ωh² increase (Boltzmann):")
print(f"    Mean factor: {np.mean(actual_factors_boltz):.4f}  (expected: {expected_factor:.4f})")
print(f"    Range: [{np.min(actual_factors_boltz):.4f}, {np.max(actual_factors_boltz):.4f}]")
print()

# Interpret
print(f"  INTERPRETATION:")
print(f"  The θ-decomposition increases Ωh² by ~{(np.mean(actual_factors_boltz)-1)*100:.0f}%.")
print(f"  From Ωh² ≈ 0.120, this gives Ωh² ≈ {0.120*np.mean(actual_factors_boltz):.3f}.")
print(f"  This means the ORIGINAL α was tuned to give Ωh²=0.120,")
print(f"  but after θ-decomposition the relic is OVERPRODUCED.")
print()
print(f"  To restore Ωh²=0.120 with θ-decomposition:")
print(f"  Need to re-tune α (or equivalently, find new BPs where")
print(f"  the ORIGINAL scan used α_s = (8/9)α in the Boltzmann equation).")
print()
print(f"  The key question: do any of these 17 BPs STILL pass both")
print(f"  SIDM (with α_s) AND relic (with the FULL α in annihilation)?")
print()

# Check with FULL α in annihilation (the physical scenario):
# ⟨σv⟩ = 2π α_s α_p / m² = 2π (8/9)(1/9) α² / m² = (16π/81) α²/m²
# vs original: πα²/(4m²)
# Ratio: (16π/81) / (π/4) = 64/81 — SAME as before
# So the relic density calculation is the same regardless of interpretation:
# the annihilation cross section is reduced by 64/81.

# Alternative interpretation: α in original scan IS α_total,
# and the annihilation process χχ→φφ uses the FULL α (not just α_s).
# Then: ⟨σv⟩ = πα²/(4m²) UNCHANGED → Ωh² unchanged → BPs still work!
# But: SIDM uses only α_s = (8/9)α → cross section drops.

print(f"  ALTERNATIVE: If χχ→φφ uses FULL α (not decomposed):")
print(f"    Relic: ⟨σv⟩ = πα²/(4m²) → Ωh² = 0.120 (unchanged)")
print(f"    SIDM:  α_s = (8/9)α       → σ/m drops (from consistency_17bp.py)")
print(f"    This scenario: relic OK, SIDM reduced but 10/17 still pass")
