#!/usr/bin/env python3
"""
Quick analysis: r = α_s/α vs SIDM viability vs vacuum energy.
V2 — adapted for corrected hyperbolic parameterization.
"""
import csv, math
from collections import defaultdict

rows = []
with open('dark_energy_exploration/vev_scan_results.csv', 'r') as f:
    for r in csv.DictReader(f):
        rows.append(r)

print(f'Total rows: {len(rows)}')

sep = '=' * 80
print()
print(sep)
print('  ANALYSIS V2: r = α_s/α vs SIDM viability vs vacuum energy')
print('  Parameterization: α_s = r α,  α_p = α/(8r),  α_s α_p = α²/8')
print(sep)

# Group by (BP, r)
bp_r = defaultdict(list)
for r in rows:
    bp_r[(r['BP'], float(r['r_alpha_s_over_alpha']))].append(r)

bps_seen = list(dict.fromkeys(r['BP'] for r in rows))

print()
print('  r = α_s/α  |  α_s/α_p ratio  |  SIDM strength  |  CW strength (α_s+α_p)')
print('  ' + '-' * 72)
r_cp = 1.0 / math.sqrt(8)
for r_val in [0.01, 0.1, r_cp, 0.5, 1.0, 2.0, 5.0, 10.0]:
    a_s_over_a_p = 8 * r_val**2   # since α_s/α_p = (r α)/(α/(8r)) = 8r²
    sum_frac = r_val + 1.0 / (8 * r_val)  # (α_s + α_p)/α
    if r_val < 0.1:
        strength = 'DEAD (Born)'
    elif r_val < 0.3:
        strength = 'WEAK'
    elif r_val < 2.0:
        strength = 'VIABLE'
    else:
        strength = 'STRONG'
    mark = ' <-- CP-sym' if abs(r_val - r_cp) < 0.01 else (' <-- natural' if abs(r_val - 1.0) < 0.01 else '')
    print(f'  {r_val:8.3f}   |  {a_s_over_a_p:10.3f}      |  {strength:12s}  |  {sum_frac:.3f}α{mark}')

print()
print('  Key: SIDM ∝ α_s = r × α.  CW loop ∝ y_s² + y_p² = 4π(α_s + α_p).')
print(f'       r = 1 → natural point.  r = {r_cp:.3f} → CP-symmetric (α_s = α_p).')
print('       Sum α_s + α_p has minimum at r = 1/√8 (AM-GM).')

# Per-BP: best r for closest-to-Lambda among nontrivial+stable
print()
print(sep)
print('  Per-BP: best r for closest-to-Lambda (nontrivial + stable)')
print(sep)
print()

for bp_label in bps_seen:
    best = None
    best_ratio = 999
    for key, rlist in bp_r.items():
        if key[0] != bp_label:
            continue
        for row in rlist:
            if row['nontrivial'] != 'True' or row['stable'] != 'True':
                continue
            lr = float(row['log10_ratio_to_Lambda'])
            if abs(lr) < abs(best_ratio):
                best_ratio = lr
                best = row
    if best:
        r_val = float(best['r_alpha_s_over_alpha'])
        alpha_s = float(best['alpha_s'])
        alpha = float(best['alpha'])
        lam_vpm = float(best['lambda_VPM'])
        regime = best['sidm_regime']
        mu3r = best['mu3_over_mphi']
        lam4 = float(best['lambda4'])
        phi_min = best['phi_min_GeV']
        print(f'  {bp_label:12s}  r={r_val:.3f}  '
              f'λ_VPM={lam_vpm:.2f} ({regime})  '
              f'log₁₀(|ΔV|/ρ_Λ)={best_ratio:.1f}  '
              f'μ₃/m_φ={mu3r}  λ₄={lam4:.2e}  '
              f'⟨φ⟩={phi_min} GeV')

# Per-BP tradeoff across all r
for bp_label in ['BP1', 'MAP']:
    print()
    print(sep)
    print(f'  {bp_label} TRADEOFF: min |log₁₀(ΔV/ρ_Λ)| across all r')
    print(sep)
    print()
    print('  r=α_s/α   λ_VPM     best_log10  SIDM_regime   Note')
    print('  ' + '-' * 68)

    r_vals_seen = sorted(set(
        float(row['r_alpha_s_over_alpha'])
        for row in rows if row['BP'] == bp_label
    ))
    for r_val in r_vals_seen:
        key = (bp_label, r_val)
        if key not in bp_r:
            continue
        rlist = bp_r[key]
        best_lr = 999
        best_row = None
        for row in rlist:
            if row['nontrivial'] == 'True' and row['stable'] == 'True':
                lr = abs(float(row['log10_ratio_to_Lambda']))
                if lr < best_lr:
                    best_lr = lr
                    best_row = row

        lam_vpm = float(rlist[0]['lambda_VPM'])
        regime = rlist[0]['sidm_regime']
        note = ''
        if abs(r_val - 1.0) < 0.05:
            note = '<-- NATURAL (α_s=α)'
        elif abs(r_val - r_cp) < 0.05:
            note = '<-- CP-SYMMETRIC'
        trivial = '(trivial)' if best_lr > 900 else ''
        print(f'  {r_val:8.4f}  {lam_vpm:8.3f}  {best_lr:10.1f}      {regime:12s} {note} {trivial}')

# Check: σv is constant across r
print()
print(sep)
print('  CONSISTENCY CHECK: ⟨σv⟩₀ should be constant across r')
print(sep)
print()
for bp_label in bps_seen:
    svs = set()
    for row in rows:
        if row['BP'] == bp_label:
            svs.add(float(row['sigma_v_s_wave_GeV2']))
    sv_list = sorted(svs)
    if len(sv_list) > 1:
        spread = (sv_list[-1] - sv_list[0]) / sv_list[0] * 100
        print(f'  {bp_label:12s}  ⟨σv⟩ range: [{sv_list[0]:.4e}, {sv_list[-1]:.4e}]  spread: {spread:.2e}%')
    else:
        print(f'  {bp_label:12s}  ⟨σv⟩ = {sv_list[0]:.4e}  (constant ✓)')

print()
print(sep)
print('  VERDICT')
print(sep)
print()
print('  The tradeoff on the relic hyperbola:')
print('  - r ≪ 1 → weak SIDM (Born), large total coupling (CW grows)')
print('  - r ≈ 1/√8 → CP-symmetric, minimal total coupling (CW minimal)')
print('  - r = 1 → natural point (SIDM unchanged, α_p = α/8)')
print('  - r ≫ 1 → strong SIDM, large total coupling (CW grows)')
print()
print('  Unlike V1, ΔV should now VARY with r because α_s+α_p changes')
print('  (minimum at CP-symmetric point, grows as r departs from 1/√8).')
print()
print('  ALL points on hyperbola have identical relic density.')
print('  Only SIDM and vacuum structure change with r.')
