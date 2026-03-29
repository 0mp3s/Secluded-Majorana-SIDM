#!/usr/bin/env python3
# ⚠️ SUPERSEDED — This script applies (8/9) reduction to α_CSV,
# which is already α_s. See test_alpha_convention.py (Test 12).
# Also: m_φ(BP1) corrected from 11.34 → 10.83 MeV (CSV source).
# Buggy output archived in archived_buggy/resonance_bp1_BUGGY.txt
"""
Resonance structure near BP1: σ/m(30) vs λ around the first Yukawa resonance.
BP1: m_chi=20.69 GeV, m_phi=11.34 MeV, λ_new=1.70 (near λ_crit=1.68).
"""
import sys, math
import numpy as np

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

try:
    from numba import jit
except ImportError:
    def jit(**kwargs):
        def decorator(func): return func
        return decorator

GEV2_TO_CM2 = 3.8938e-28
GEV_IN_G    = 1.78266e-24
C_KM_S      = 299792.458

@jit(nopython=True, cache=True)
def sph_jn_numba(l, z):
    if z < 1e-30: return 1.0 if l == 0 else 0.0
    j0 = math.sin(z) / z
    if l == 0: return j0
    j1 = math.sin(z) / (z * z) - math.cos(z) / z
    if l == 1: return j1
    jp, jc = j0, j1
    for n in range(1, l):
        jn = (2*n+1)/z*jc - jp; jp, jc = jc, jn
        if abs(jc) < 1e-300: return 0.0
    return jc

@jit(nopython=True, cache=True)
def sph_yn_numba(l, z):
    if z < 1e-30: return -1e300
    y0 = -math.cos(z) / z
    if l == 0: return y0
    y1 = -math.cos(z) / (z * z) - math.sin(z) / z
    if l == 1: return y1
    yp, yc = y0, y1
    for n in range(1, l):
        yn = (2*n+1)/z*yc - yp; yp, yc = yc, yn
        if abs(yc) > 1e200: return yc
    return yc

@jit(nopython=True, cache=True)
def _vpm_rhs(l, kappa, lam, x, delta):
    if x < 1e-20: return 0.0
    z = kappa * x
    if z < 1e-20: return 0.0
    jl, nl = sph_jn_numba(l, z), sph_yn_numba(l, z)
    j_hat, n_hat = z * jl, -z * nl
    cd, sd = math.cos(delta), math.sin(delta)
    b = j_hat * cd - n_hat * sd
    if not math.isfinite(b): return 0.0
    v = lam * math.exp(-x) / (kappa * x) * b * b
    return v if math.isfinite(v) else 0.0

@jit(nopython=True, cache=True)
def vpm_phase_shift(l, kappa, lam, x_max=50.0, N=4000):
    if lam < 1e-30 or kappa < 1e-30: return 0.0
    xm = max(1e-5, 0.05/(kappa+0.01))
    if l > 0:
        xb = l/kappa
        if xb > xm: xm = xb
    h = (x_max - xm) / N
    d = 0.0
    for i in range(N):
        x = xm + i*h
        k1 = _vpm_rhs(l, kappa, lam, x, d)
        k2 = _vpm_rhs(l, kappa, lam, x+.5*h, d+.5*h*k1)
        k3 = _vpm_rhs(l, kappa, lam, x+.5*h, d+.5*h*k2)
        k4 = _vpm_rhs(l, kappa, lam, x+h, d+h*k3)
        d += h*(k1+2*k2+2*k3+k4)/6.0
    return d

@jit(nopython=True, cache=True)
def sigma_T_vpm(m_chi, m_phi, alpha, v_km_s):
    v = v_km_s / C_KM_S
    k = m_chi/2.0 * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15: return 0.0
    if kappa < 5:     xm, N = 50.0, 4000
    elif kappa < 50:  xm, N = 80.0, 8000
    else:             xm, N = 100.0, 12000
    lmax = min(max(3, min(int(kappa*xm), int(kappa)+int(lam)+20)), 500)
    s, pk, ns = 0.0, 0.0, 0
    for l in range(lmax+1):
        d = vpm_phase_shift(l, kappa, lam, xm, N)
        w = 1.0 if l%2==0 else 3.0
        c = w*(2*l+1)*math.sin(d)**2
        s += c
        if c > pk: pk = c
        if pk > 0 and c/pk < 1e-4:
            ns += 1
            if ns >= 5: break
        else: ns = 0
    sg = 2.0*math.pi*s/(k*k)
    return sg * GEV2_TO_CM2 / (m_chi * GEV_IN_G)

# ==============================================================
#  BP1 parameters
# ==============================================================
M_CHI = 20.69       # GeV
M_PHI = 10.83e-3    # GeV  (corrected from 11.34; CSV source)
ALPHA_ORIG = 1.048e-3
LAM_ORIG = ALPHA_ORIG * M_CHI / M_PHI  # 1.91

print("=" * 80)
print("  RESONANCE STRUCTURE NEAR BP1")
print("=" * 80)
print(f"  BP1: m_χ = {M_CHI} GeV, m_φ = {M_PHI*1e3} MeV")
print(f"  λ_orig = {LAM_ORIG:.3f}, λ_decomposed = {LAM_ORIG*8/9:.3f}")
print(f"  First Yukawa resonance: λ_crit ≈ 1.68")
print()

# JIT warmup
_ = sigma_T_vpm(20.0, 0.01, 1e-3, 30.0)

# Scan λ from 0.5 to 3.5 by varying α (keeping m_chi, m_phi fixed)
print("  σ/m(30 km/s) vs λ  [m_χ, m_φ fixed, α varied]")
print("  " + "-" * 75)

lam_vals = np.linspace(0.5, 3.5, 100)
sm_vals = []

for lam in lam_vals:
    alpha = lam * M_PHI / M_CHI
    sm = sigma_T_vpm(M_CHI, M_PHI, alpha, 30.0)
    sm_vals.append(sm)

sm_vals = np.array(sm_vals)

# Find key points
idx_max = np.argmax(sm_vals)
lam_peak = lam_vals[idx_max]
sm_peak = sm_vals[idx_max]

# Print ASCII plot
print()
print(f"  Peak: λ={lam_peak:.2f}, σ/m={sm_peak:.3f} cm²/g")
print()

# Mark special λ values
lam_decomp = LAM_ORIG * 8/9  # 1.70
alpha_decomp = lam_decomp * M_PHI / M_CHI

# Fine scan near resonance
print("  FINE SCAN near first resonance (λ = 1.4 to 2.2):")
print(f"  {'λ':>6}  {'α':>10}  {'σ/m(30)':>10}  {'notes':>30}")
print("  " + "-" * 65)

fine_lams = np.linspace(1.4, 2.2, 50)
fine_sms = []
for lam in fine_lams:
    alpha = lam * M_PHI / M_CHI
    sm = sigma_T_vpm(M_CHI, M_PHI, alpha, 30.0)
    fine_sms.append(sm)
    
    notes = ""
    if abs(lam - 1.68) < 0.02: notes = "← λ_crit (1st resonance)"
    elif abs(lam - lam_decomp) < 0.02: notes = "← BP1 decomposed (8/9)"
    elif abs(lam - LAM_ORIG) < 0.02: notes = "← BP1 original"
    
    if notes:
        tag = "***" if sm >= 0.5 else "   "
        print(f"  {lam:6.3f}  {alpha:.3e}  {sm:10.4f}  {notes} {tag}")

fine_sms = np.array(fine_sms)

# ASCII chart
print()
print("  σ/m(30) profile near resonance:")
print("  " + "-" * 65)
max_bar = 50
for i, (lam, sm) in enumerate(zip(fine_lams, fine_sms)):
    if i % 2 != 0 and abs(lam - 1.68) > 0.02 and abs(lam - lam_decomp) > 0.02 and abs(lam - LAM_ORIG) > 0.02:
        continue
    bar_len = int(sm / max(fine_sms) * max_bar)
    bar = "█" * bar_len
    marker = ""
    if abs(lam - 1.68) < 0.02: marker = " ◄ λ_crit"
    elif abs(lam - lam_decomp) < 0.02: marker = " ◄ BP1(8/9)"
    elif abs(lam - LAM_ORIG) < 0.02: marker = " ◄ BP1(orig)"
    thresh = "│" if bar_len >= int(0.5 / max(fine_sms) * max_bar) else " "
    print(f"  λ={lam:.2f} {bar}{marker}")

# The physical question: how far from resonance peak?
print()
print("=" * 80)
print("  QUANTUM RESONANCE ANALYSIS")
print("=" * 80)
print()

# Compute phase shift δ₀ at 30 km/s for decomposed BP1
v = 30.0 / C_KM_S
k = M_CHI / 2.0 * v
kappa = k / M_PHI

# At decomposed λ
lam_d = alpha_decomp * M_CHI / M_PHI
delta_0_decomp = vpm_phase_shift(0, kappa, lam_d)
delta_0_orig = vpm_phase_shift(0, kappa, LAM_ORIG * M_PHI / M_CHI * M_CHI / M_PHI)

# At resonance peak
alpha_peak = lam_peak * M_PHI / M_CHI
delta_0_peak = vpm_phase_shift(0, kappa, lam_peak)

print(f"  s-wave phase shift δ₀ at v=30 km/s:")
print(f"    BP1 original  (λ={LAM_ORIG:.2f}):  δ₀ = {delta_0_orig:.4f} rad = {math.degrees(delta_0_orig):.2f}°")
print(f"    BP1 decomposed(λ={lam_decomp:.2f}):  δ₀ = {delta_0_decomp:.4f} rad = {math.degrees(delta_0_decomp):.2f}°")
print(f"    Peak          (λ={lam_peak:.2f}):  δ₀ = {delta_0_peak:.4f} rad = {math.degrees(delta_0_peak):.2f}°")
print()
print(f"  sin²(δ₀): original = {math.sin(delta_0_orig)**2:.4f}")
print(f"  sin²(δ₀): decomposed = {math.sin(delta_0_decomp)**2:.4f}")
print(f"  sin²(δ₀): peak = {math.sin(delta_0_peak)**2:.4f}")
print(f"  Resonance peak: sin²δ₀ = 1 when δ₀ = π/2 = {math.pi/2:.4f} rad")
print()

# Distance to resonance in phase-shift space
dist_to_res = abs(delta_0_decomp - math.pi/2)
print(f"  Distance to resonance: |δ₀ - π/2| = {dist_to_res:.4f} rad = {math.degrees(dist_to_res):.2f}°")
print()

# What α would put us exactly on resonance?
print("  Searching for α that maximizes σ/m(30) [resonance peak]...")
best_sm = 0
best_alpha = 0
for lam_try in np.linspace(1.5, 2.0, 500):
    a = lam_try * M_PHI / M_CHI
    sm = sigma_T_vpm(M_CHI, M_PHI, a, 30.0)
    if sm > best_sm:
        best_sm = sm
        best_alpha = a
        best_lam = lam_try

print(f"  Resonance peak: λ = {best_lam:.4f}, α = {best_alpha:.4e}")
print(f"  σ/m(30) at peak = {best_sm:.3f} cm²/g")
print(f"  σ/m(30) at decomposed BP1 = {sigma_T_vpm(M_CHI, M_PHI, alpha_decomp, 30.0):.3f} cm²/g")
print(f"  Ratio: peak / decomposed = {best_sm / sigma_T_vpm(M_CHI, M_PHI, alpha_decomp, 30.0):.2f}×")
print()

# Could quantum tunneling through the resonance boost the effective cross section?
print("  PHYSICAL INTERPRETATION:")
print("  " + "-" * 65)
print(f"  BP1 decomposed sits at λ = {lam_decomp:.2f}")
print(f"  Nearest resonance peak at λ = {best_lam:.4f}")
print(f"  Distance: Δλ = {abs(lam_decomp - best_lam):.4f}")
print(f"  Resonance width (FWHM) from scan: ", end="")

# Estimate resonance width
half_max = best_sm / 2
above_half = fine_lams[fine_sms > half_max] if any(fine_sms > half_max) else []
if len(above_half) >= 2:
    width = above_half[-1] - above_half[0]
    print(f"Δλ_FWHM ≈ {width:.2f}")
    print(f"  BP1 is {abs(lam_decomp - best_lam)/width:.1f} × FWHM from peak")
else:
    print("(no clear resonance peak in this range)")

print()
