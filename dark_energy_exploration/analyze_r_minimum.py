#!/usr/bin/env python3
"""
Focused analysis: Why is the ΔV minimum at r≈1.17 and not r=1 for BP1?

Physics question from Sonnet:
  "מה קורה ב-r≈1.17 לעומת r=1? ההפרש קטן — אבל האם יש פרשנות פיזיקלית
   למה המינימום נמצא שם ולא ב-r=1?"

Strategy:
  1. Extract all BP1 rows, group by r, find best |ΔV/ρΛ| per r
  2. Check: is r=1.17 an artifact of the discrete grid, or a true minimum?
  3. Compute the analytic dependence: α_s + α_p = rα + α/(8r)
     has derivative d/dr = α - α/(8r²), zero at r = 1/√8 ≈ 0.354
     → the TOTAL coupling sum is monotonically increasing for r > 0.354
     → so the CW loop GROWS with r for r > 0.354
     → ΔV should get LARGER (more negative) with increasing r
     → but the tree-level also depends on μ₃, λ₄
  4. For fixed (μ₃, λ₄), how does ΔV depend on r?
  5. Fine-grained r scan around [0.8, 1.5] for BP1
"""
import csv
import math
import sys
import numpy as np
from pathlib import Path
from collections import defaultdict

# Also run a fine-grained scan
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "core"))
from global_config import GC

# Import physics from vev_scan
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vev_scan import build_V_eff, find_minimum

CSV_PATH = Path(__file__).resolve().parent / "vev_scan_results.csv"

sep = "=" * 80

# ═══════════════════════════════════════════════════════════
#  Part 1: Analyze existing grid data for BP1
# ═══════════════════════════════════════════════════════════

rows = []
with open(CSV_PATH, 'r') as f:
    for r in csv.DictReader(f):
        rows.append(r)

bp1_rows = [r for r in rows if r['BP'] == 'BP1']
print(f"BP1 rows: {len(bp1_rows)}")

# Group by r, find minimum |log10(|ΔV|/ρΛ)| per r (nontrivial+stable)
r_best = defaultdict(lambda: {'best_lr': 999, 'row': None, 'all_lr': []})
for r in bp1_rows:
    rv = float(r['r_alpha_s_over_alpha'])
    if r['nontrivial'] == 'True' and r['stable'] == 'True':
        lr = abs(float(r['log10_ratio_to_Lambda']))
        r_best[rv]['all_lr'].append(lr)
        if lr < r_best[rv]['best_lr']:
            r_best[rv]['best_lr'] = lr
            r_best[rv]['row'] = r

print()
print(sep)
print("  BP1: Best |log₁₀(|ΔV|/ρ_Λ)| per r value")
print(sep)
print()
print(f"  {'r':>8s}  {'α_s+α_p':>10s}  {'λ_VPM':>8s}  {'regime':>12s}  "
      f"{'best_log':>10s}  {'median_log':>10s}  {'n_stable':>8s}")
print("  " + "-" * 80)

for rv in sorted(r_best.keys()):
    info = r_best[rv]
    if info['row'] is None:
        continue
    row = info['row']
    alpha = float(row['alpha'])
    a_s = float(row['alpha_s'])
    a_p = float(row['alpha_p'])
    sum_a = a_s + a_p
    sum_ratio = sum_a / alpha  # (α_s + α_p)/α
    lam_vpm = float(row['lambda_VPM'])
    regime = row['sidm_regime']
    all_lr = sorted(info['all_lr'])
    median_lr = all_lr[len(all_lr)//2] if all_lr else 999

    print(f"  {rv:8.4f}  {sum_ratio:10.4f}α  {lam_vpm:8.3f}  {regime:>12s}  "
          f"{info['best_lr']:10.1f}  {median_lr:10.1f}  {len(all_lr):>8d}")

# ═══════════════════════════════════════════════════════════
#  Part 2: Analytic understanding
# ═══════════════════════════════════════════════════════════

print()
print(sep)
print("  ANALYTIC: Why r≈1.17?")
print(sep)
print()
print("  The CW loop depends on y_s² + y_p² = 4π(α_s + α_p) = 4πα(r + 1/(8r))")
print("  The sum f(r) = r + 1/(8r) has:")
print("    f'(r) = 1 - 1/(8r²) = 0  →  r* = 1/√8 ≈ 0.354  (minimum)")
print("    f(0.354) = 0.707  (CP-symmetric point)")
print("    f(1.0)   = 1.125  (natural point)")
print("    f(1.17)  = 1.277")
print()

# But ΔV depends on the INTERPLAY between tree-level and CW
# Let's check: for a fixed μ₃/m_φ and λ₄, what r minimizes |ΔV|?

print("  The tree-level V₀ = ½m²_φ φ² + (μ₃/6)φ³ + (λ₄/24)φ⁴")
print("  does NOT depend on r at all!")
print()
print("  The CW correction δV depends on y_s and y_p separately:")
print("    M²_eff(φ) = (m_χ + y_s φ/2)² + (y_p φ/2)²")
print("  So it depends on BOTH y_s and y_p, not just their sum.")
print()
print("  At r=1:   y_s = √(4πα),       y_p = √(πα/2)       → y_s/y_p = 2√2 ≈ 2.83")
print("  At r=1.17: y_s = √(4π·1.17α),  y_p = √(πα/(2·1.17)) → y_s/y_p = 2√2·√1.17 ≈ 3.06")
print()

# ═══════════════════════════════════════════════════════════
#  Part 3: Fine-grained r scan for BP1 at the best (μ₃, λ₄)
# ═══════════════════════════════════════════════════════════

print(sep)
print("  FINE SCAN: BP1, r ∈ [0.5, 2.0], 100 points")
print("  at (μ₃/m_φ, λ₄) that gave the minimum in the coarse scan")
print(sep)
print()

# Get BP1 params
bp1 = [b for b in GC.all_benchmarks() if b['label'] == 'BP1'][0]
m_chi = bp1['m_chi_GeV']
m_phi = bp1['m_phi_MeV'] * 1e-3
alpha = bp1['alpha']

RHO_LAMBDA = 2.58e-47

# Best row from coarse scan
best_coarse = r_best[sorted(r_best.keys(), key=lambda x: r_best[x]['best_lr'])[0]]['row']
mu3r_best = float(best_coarse['mu3_over_mphi'])
lam4_best = float(best_coarse['lambda4'])
print(f"  Best coarse: r={float(best_coarse['r_alpha_s_over_alpha']):.4f}, "
      f"μ₃/m_φ={mu3r_best:.2f}, λ₄={lam4_best:.2e}")
print()

mu3 = mu3r_best * m_phi

# Fine scan
r_fine = np.linspace(0.5, 2.0, 200)
print(f"  {'r':>8s}  {'α_s/α':>8s}  {'α_p/α':>8s}  {'(αs+αp)/α':>10s}  "
      f"{'y_s':>10s}  {'y_p':>10s}  {'ΔV (GeV⁴)':>14s}  {'log₁₀|ΔV/ρΛ|':>14s}  "
      f"{'⟨φ⟩ (GeV)':>12s}")
print("  " + "-" * 110)

best_fine_lr = 999
best_fine_r = 0
results_fine = []

for r in r_fine:
    a_s = r * alpha
    a_p = alpha**2 / (8.0 * a_s)
    y_s = math.sqrt(4.0 * math.pi * a_s)
    y_p = math.sqrt(4.0 * math.pi * a_p)

    V_vec, V_scl = build_V_eff(m_chi, m_phi, y_s, y_p, mu3, lam4_best)
    phi_range = max(100.0, 5.0 * m_chi)
    phi_min, V_min, stable = find_minimum(V_vec, V_scl, phi_range)

    V_at_0 = V_scl(0.0)
    DeltaV = V_min - V_at_0
    nontrivial = abs(phi_min) > 1e-10

    if nontrivial and DeltaV != 0:
        abs_dv = abs(DeltaV)
        log_ratio = math.log10(abs_dv / RHO_LAMBDA)
        results_fine.append((r, a_s, a_p, y_s, y_p, DeltaV, log_ratio, phi_min))

        if abs(log_ratio) < abs(best_fine_lr):
            best_fine_lr = log_ratio
            best_fine_r = r

# Print every 5th point
for i, (r, a_s, a_p, y_s, y_p, DeltaV, log_ratio, phi_min) in enumerate(results_fine):
    if i % 5 == 0 or abs(r - 1.0) < 0.01 or abs(r - best_fine_r) < 0.01:
        mark = ""
        if abs(r - 1.0) < 0.02:
            mark = " <-- r=1 (natural)"
        if abs(r - best_fine_r) < 0.02:
            mark = " <-- MINIMUM"
        print(f"  {r:8.4f}  {a_s/alpha:8.4f}  {a_p/alpha:8.4f}  "
              f"{(a_s+a_p)/alpha:10.4f}  {y_s:10.6f}  {y_p:10.6f}  "
              f"{DeltaV:14.4e}  {log_ratio:14.2f}  {phi_min:12.4e}{mark}")

print()
print(f"  Fine scan minimum: r = {best_fine_r:.4f},  log₁₀|ΔV/ρΛ| = {best_fine_lr:.2f}")
print()

# ═══════════════════════════════════════════════════════════
#  Part 4: WHY does the minimum occur there?
# ═══════════════════════════════════════════════════════════

print(sep)
print("  PHYSICAL INTERPRETATION")
print(sep)
print()
print("  The effective potential V_eff = V_tree + δV_CW depends on r through:")
print("    1. V_tree does NOT depend on r (it's pure φ self-interaction)")
print("    2. δV_CW depends on M²_eff(φ) = (m_χ + y_s φ/2)² + (y_p φ/2)²")
print()
print("  The VEV ⟨φ⟩ is determined by V'_eff(⟨φ⟩) = 0.")
print("  Since V_tree dominates at large φ and δV_CW is a small correction,")
print("  the VEV is approximately set by V_tree, and ΔV ≈ δV_CW(⟨φ_tree⟩).")
print()
print("  At the VEV, δV_CW depends on:")
print("    M²_eff(⟨φ⟩) = m²_χ[1 + y_s⟨φ⟩/m_χ + (y²_s + y²_p)⟨φ⟩²/(4m²_χ)]")
print("                 ≈ m²_χ[1 + y_s⟨φ⟩/m_χ]  (since ⟨φ⟩ ~ 10⁻⁶ GeV ≪ m_χ)")
print()
print("  So δV_CW ∝ M⁴_eff ln(M²_eff/m²_χ), and the leading r-dependence")
print("  comes from the y_s⟨φ⟩/m_χ cross-term, which is ∝ √(r).")
print()
print("  The competition:")
print("    - Larger r → larger y_s → CW pushes φ further from 0 → deeper ΔV")
print("    - Larger r → larger α_s+α_p → CW mass correction larger → shallower ΔV")
print("    - The two effects balance at r_min")
print()

# Check: is r_min independent of (μ₃, λ₄)?
print(sep)
print("  CHECK: Does r_min depend on (μ₃/m_φ, λ₄)?")
print(sep)
print()

test_params = [
    (1.7, 1e-4),
    (1.7, 1e-2),
    (4.86, 8.46e-4),  # best from coarse scan
    (10.0, 1e-2),
    (40.0, 1.0),
    (80.0, 4*math.pi),
]

print(f"  {'μ₃/m_φ':>8s}  {'λ₄':>10s}  {'r_min':>8s}  {'log₁₀|ΔV/ρΛ|':>14s}  {'⟨φ⟩':>12s}")
print("  " + "-" * 60)

for mu3r_test, lam4_test in test_params:
    mu3_test = mu3r_test * m_phi
    best_lr_t = 999
    best_r_t = 0
    best_phi_t = 0

    for r in np.linspace(0.3, 3.0, 100):
        a_s = r * alpha
        a_p = alpha**2 / (8.0 * a_s)
        y_s = math.sqrt(4.0 * math.pi * a_s)
        y_p = math.sqrt(4.0 * math.pi * a_p)

        V_vec, V_scl = build_V_eff(m_chi, m_phi, y_s, y_p, mu3_test, lam4_test)
        phi_range = max(100.0, 5.0 * m_chi)
        phi_min, V_min, stable = find_minimum(V_vec, V_scl, phi_range)

        V_at_0 = V_scl(0.0)
        DeltaV = V_min - V_at_0
        nontrivial = abs(phi_min) > 1e-10

        if nontrivial and DeltaV != 0:
            abs_dv = abs(DeltaV)
            lr = math.log10(abs_dv / RHO_LAMBDA)
            if abs(lr) < abs(best_lr_t):
                best_lr_t = lr
                best_r_t = r
                best_phi_t = phi_min

    print(f"  {mu3r_test:8.2f}  {lam4_test:10.2e}  {best_r_t:8.4f}  "
          f"{best_lr_t:14.2f}  {best_phi_t:12.4e}")

print()
print(sep)
print("  CONCLUSION")
print(sep)
print()
print("  If r_min ≈ 1.17 for ALL (μ₃, λ₄), it's a CW-loop property:")
print("    → intrinsic to the M²_eff structure, not the tree-level potential.")
print("    → physical meaning: the ratio y_s/y_p that minimizes |ΔV|.")
print()
print("  If r_min varies with (μ₃, λ₄), it's a tree-CW interference effect:")
print("    → the minimum is an accident of the specific parameter point.")
print("    → no deep physical significance.")
