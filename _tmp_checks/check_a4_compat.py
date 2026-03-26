#!/usr/bin/env python3
"""Quick check: how many of the 122 relic-viable points are 
compatible with the A4 constraint alpha_p = alpha_s / 8?

Logic:
  - v27 uses <sigma v> = pi alpha^2 / (4 m_chi^2) with single alpha
  - For mixed coupling: <sigma v> = pi alpha_s alpha_p / (4 m_chi^2) 
  - So the relic product P = alpha_s * alpha_p = alpha_old^2 (approximately)
  - BUT the CP-separation config gives relic_product = 1.387474e-7 for BP1 
    while alpha_BP1^2 = 1.098e-6 — factor ~8 difference
  - This factor arises because for Majorana: only the y_s^2 y_p^2 cross-term
    contributes to s-wave. The single-coupling formula uses y^4, which is 
    the CP-symmetric limit where y_s = y_p = y, giving y_s^2 y_p^2 = y^4/4.
    Combined with the 1/2 Majorana factor: relic_product = alpha_old^2 / 8
  - We verify: (1.048e-3)^2 / 8 = 1.373e-7 ≈ 1.387e-7 ✓ (1% from m_phi effects)
  
  For each 122 point: 
    P(m_chi) = alpha_old^2 / 8  (where alpha_old is the "alpha" in v31 CSV)
    A4 constraint: alpha_p = alpha_s / 8, so alpha_s^2 / 8 = P
    => alpha_s_A4 = alpha_old  (!)
    
  Wait that gives alpha_s_A4 = alpha_old. Let me re-derive:
    P = alpha_old^2 / 8
    A4: alpha_s * alpha_s/8 = P => alpha_s^2 = 8P = alpha_old^2
    => alpha_s_A4 = alpha_old
    
  So A4 predicts alpha_s = alpha_old for every point. The SIDM cross sections 
  would be IDENTICAL to what's already in the CSV (since sigma_T depends on alpha_s).
  
  This means ALL 122 points are A4-compatible by construction — the alpha in the CSV
  IS the alpha_s that A4 would predict!
  
  ... unless the factor-of-8 relationship is not exact. Let me check more carefully.
"""
import sys, os, csv, math
import numpy as np

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

from output_manager import get_latest

# Load the 122-point CSV
v31_path = get_latest("v31_true_viable_points")
print(f"Loading: {v31_path}")

rows = []
with open(v31_path) as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append({
            'm_chi': float(r['m_chi_GeV']),
            'm_phi': float(r['m_phi_MeV']) * 1e-3,  # GeV
            'alpha': float(r['alpha']),
            'lambda': float(r['lambda']),
            'sigma_30': float(r['sigma_m_30']),
            'sigma_1000': float(r['sigma_m_1000']),
        })

print(f"Loaded {len(rows)} relic-viable points\n")

# Known relic product for BP1 masses
# CP-separation config: relic_product = 1.387474e-7
# BP1: m_chi = 20.6914, m_phi = 11.3379 MeV, alpha = 1.048412e-3
P_BP1 = 1.387474e-7
alpha_BP1 = 1.048412e-3
m_chi_BP1 = 20.6914

# Check the factor-of-8 hypothesis
print("=== Factor-of-8 Check ===")
print(f"alpha_BP1^2 = {alpha_BP1**2:.6e}")
print(f"alpha_BP1^2 / 8 = {alpha_BP1**2 / 8:.6e}")
print(f"relic_product_BP1 = {P_BP1:.6e}")
print(f"ratio alpha^2 / (8 * P) = {alpha_BP1**2 / (8*P_BP1):.4f}")
print()

# If P = alpha^2 / 8 exactly, then A4 gives alpha_s = alpha.
# But there's a ~1% correction from m_phi. Let's check:
# For A4: alpha_s^2 / 8 = P => alpha_s = sqrt(8P)
# For BP1: alpha_s_A4 = sqrt(8 * 1.387e-7) = sqrt(1.110e-6) = 1.0535e-3
alpha_s_A4_BP1 = math.sqrt(8 * P_BP1)
print(f"A4 alpha_s for BP1 = {alpha_s_A4_BP1:.6e}")
print(f"Actual alpha_BP1 = {alpha_BP1:.6e}")
print(f"Deviation: {abs(alpha_s_A4_BP1 - alpha_BP1)/alpha_BP1 * 100:.2f}%")
print()

# Now: for the general case, the relic product P(m_chi, m_phi) scales approximately
# as m_chi^2 (s-wave thermal relic). But there are corrections from g_*(T_fo), m_phi, etc.
# 
# Without re-running Boltzmann for each point, we can estimate:
# P(m_chi) ≈ P_BP1 * (m_chi / m_chi_BP1)^2
#
# Then A4 alpha_s = sqrt(8 * P(m_chi)) = sqrt(8 * P_BP1) * (m_chi / m_chi_BP1)
# = alpha_s_A4_BP1 * (m_chi / m_chi_BP1)
#
# The actual alpha in CSV is the one that gives omega_h2 = 0.12 with the OLD formula.
# If alpha_old^2/8 = P(m_chi), then alpha_old = sqrt(8*P(m_chi)).
#
# But P(m_chi) was computed by Boltzmann with the full formula, not the m_chi^2 approx.
# 
# Bottom line: alpha_s_A4 ≈ alpha_old * correction_factor, where correction comes from 
# the ratio between the exact relic product and the alpha^2/8 approximation.
# This correction is ~0.5% for BP1 and may vary across mass range.

# For each of 122 points, estimate whether A4 is compatible:
# We use the m_chi^2 scaling approximation
print("=== 122-Point A4 Compatibility ===")
print(f"{'m_chi':>8s} {'m_phi_MeV':>9s} {'alpha_CSV':>12s} {'alpha_A4_est':>12s} {'dev%':>8s} {'lambda':>8s} {'lam_A4':>8s} {'sigma30':>8s} {'Fornax12':>8s}")
print("-" * 100)

fornax_fail_a4_compat = 0
fornax_fail_a4_incompat = 0
fornax_pass_a4_compat = 0
fornax_pass_a4_incompat = 0

a4_compatible = []
a4_incompatible = []

for row in rows:
    m_chi = row['m_chi']
    m_phi = row['m_phi']
    alpha_csv = row['alpha']
    lam = row['lambda']
    
    # Estimate relic product using m_chi^2 scaling
    P_est = P_BP1 * (m_chi / m_chi_BP1)**2
    
    # A4 constraint: alpha_s = sqrt(8 * P)
    alpha_A4 = math.sqrt(8 * P_est)
    
    # Deviation
    dev = (alpha_csv - alpha_A4) / alpha_A4 * 100
    
    # Lambda at A4 point
    lam_A4 = alpha_A4 * m_chi / m_phi
    
    # Estimate sigma/m at 12 km/s (we don't have it in CSV, use lambda as proxy)
    # Points with lambda > ~15 tend to fail Fornax
    
    # Define "A4 compatible" as |dev| < 15%
    is_compat = abs(dev) < 15
    
    if is_compat:
        a4_compatible.append(row)
    else:
        a4_incompatible.append(row)

# Print summary stats
print(f"\nA4-compatible (|dev| < 15%): {len(a4_compatible)} / {len(rows)}")
print(f"A4-incompatible: {len(a4_incompatible)} / {len(rows)}")

# Check lambda distribution
lam_compat = [r['lambda'] for r in a4_compatible]
lam_incompat = [r['lambda'] for r in a4_incompatible]

print(f"\nA4-compatible lambda range: [{min(lam_compat):.1f}, {max(lam_compat):.1f}]" if lam_compat else "")
print(f"A4-incompatible lambda range: [{min(lam_incompat):.1f}, {max(lam_incompat):.1f}]" if lam_incompat else "")

# More useful: what percentage of A4-compatible points have sigma_30 in [0.5, 5]?
# And what's their lambda range vs the Fornax-failing ones?
print("\n=== Detailed Breakdown ===")

# Check the deviation distribution
devs = []
for row in rows:
    m_chi = row['m_chi']
    P_est = P_BP1 * (m_chi / m_chi_BP1)**2
    alpha_A4 = math.sqrt(8 * P_est)
    dev = (row['alpha'] - alpha_A4) / alpha_A4 * 100
    devs.append(dev)

devs = np.array(devs)
print(f"Deviation stats: mean={devs.mean():.1f}%, std={devs.std():.1f}%, min={devs.min():.1f}%, max={devs.max():.1f}%")

# Actually more useful: check which points sit near A4 and which don't
# Group by mass point (m_chi, m_phi) — multiple alpha values per mass pair
from collections import defaultdict
mass_groups = defaultdict(list)
for row in rows:
    key = (row['m_chi'], row['m_phi'])
    mass_groups[key].append(row)

print(f"\nUnique mass points: {len(mass_groups)}")
print(f"Mass points with >1 alpha value: {sum(1 for v in mass_groups.values() if len(v) > 1)}")

# For each mass pair, what's the range of alpha values?
print(f"\n{'m_chi':>8s} {'m_phi_MeV':>9s} {'n_pts':>5s} {'alpha_min':>12s} {'alpha_max':>12s} {'alpha_A4':>12s} {'A4_in_range':>11s}")
for (mc, mp), group in sorted(mass_groups.items()):
    alphas = [r['alpha'] for r in group]
    P_est = P_BP1 * (mc / m_chi_BP1)**2
    a4 = math.sqrt(8 * P_est)
    in_range = min(alphas) <= a4 <= max(alphas)
    # Only a4 is meaningful if it's actually achievable
    print(f"{mc:8.2f} {mp*1e3:9.2f} {len(group):5d} {min(alphas):12.6e} {max(alphas):12.6e} {a4:12.6e} {'YES' if in_range else 'no':>11s}")
