#!/usr/bin/env python3
# ⚠️ SUPERSEDED — This script applies (8/9) reduction to α_CSV,
# which is already α_s. See test_alpha_convention.py (Test 12).
# Qualitative trends valid; absolute σ/m values are ~11% low.
# Buggy output archived in archived_buggy/free_theta_scan_BUGGY.txt
"""
Free θ scan: σ/m(30) as function of θ for all 17 BPs.
Question: Does Yukawa resonance structure prefer a specific θ?
"""
import sys, math, os
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
M_PL        = 1.220890e19
S_0         = 2891.2
RHO_CRIT_H2 = 1.0539e-5

# VPM (same as before, compact)
@jit(nopython=True, cache=True)
def sph_jn(l, z):
    if z < 1e-30: return 1.0 if l == 0 else 0.0
    j0 = math.sin(z)/z
    if l == 0: return j0
    j1 = math.sin(z)/(z*z) - math.cos(z)/z
    if l == 1: return j1
    jp, jc = j0, j1
    for n in range(1,l):
        jn = (2*n+1)/z*jc - jp; jp, jc = jc, jn
        if abs(jc) < 1e-300: return 0.0
    return jc

@jit(nopython=True, cache=True)
def sph_yn(l, z):
    if z < 1e-30: return -1e300
    y0 = -math.cos(z)/z
    if l == 0: return y0
    y1 = -math.cos(z)/(z*z) - math.sin(z)/z
    if l == 1: return y1
    yp, yc = y0, y1
    for n in range(1,l):
        yn = (2*n+1)/z*yc - yp; yp, yc = yc, yn
        if abs(yc) > 1e200: return yc
    return yc

@jit(nopython=True, cache=True)
def _rhs(l, kappa, lam, x, delta):
    if x < 1e-20: return 0.0
    z = kappa*x
    if z < 1e-20: return 0.0
    jl, nl = sph_jn(l,z), sph_yn(l,z)
    jh, nh = z*jl, -z*nl
    cd, sd = math.cos(delta), math.sin(delta)
    b = jh*cd - nh*sd
    if not math.isfinite(b): return 0.0
    v = lam*math.exp(-x)/(kappa*x)*b*b
    return v if math.isfinite(v) else 0.0

@jit(nopython=True, cache=True)
def phase_shift(l, kappa, lam, xm=50.0, N=4000):
    if lam < 1e-30 or kappa < 1e-30: return 0.0
    x0 = max(1e-5, 0.05/(kappa+0.01))
    if l > 0:
        xb = l/kappa
        if xb > x0: x0 = xb
    h = (xm-x0)/N; d = 0.0
    for i in range(N):
        x = x0+i*h
        k1 = _rhs(l,kappa,lam,x,d)
        k2 = _rhs(l,kappa,lam,x+.5*h,d+.5*h*k1)
        k3 = _rhs(l,kappa,lam,x+.5*h,d+.5*h*k2)
        k4 = _rhs(l,kappa,lam,x+h,d+h*k3)
        d += h*(k1+2*k2+2*k3+k4)/6.0
    return d

@jit(nopython=True, cache=True)
def sigma_T(m_chi, m_phi, alpha, v_km):
    v = v_km/C_KM_S; k = m_chi/2*v; kappa = k/m_phi
    lam = alpha*m_chi/m_phi
    if kappa < 1e-15: return 0.0
    if kappa < 5:     xm, N = 50.0, 4000
    elif kappa < 50:  xm, N = 80.0, 8000
    else:             xm, N = 100.0, 12000
    lx = min(max(3, min(int(kappa*xm), int(kappa)+int(lam)+20)), 500)
    s, pk, ns = 0.0, 0.0, 0
    for l in range(lx+1):
        d = phase_shift(l,kappa,lam,xm,N)
        w = 1.0 if l%2==0 else 3.0
        c = w*(2*l+1)*math.sin(d)**2
        s += c
        if c > pk: pk = c
        if pk > 0 and c/pk < 1e-4:
            ns += 1
            if ns >= 5: break
        else: ns = 0
    sg = 2.0*math.pi*s/(k*k)
    return sg*GEV2_TO_CM2/(m_chi*GEV_IN_G)

# Relic
_GT = np.array([[1e4,106.75,106.75],[200,106.75,106.75],[80,86.25,86.25],
    [10,86.25,86.25],[1,75.75,75.75],[0.3,61.75,61.75],[0.2,17.25,17.25],
    [0.15,14.25,14.25],[0.1,10.75,10.75],[0.01,10.75,10.75],
    [0.001,10.75,10.75],[0.0005,10.75,10.75],[0.0001,3.36,3.91],
    [1e-5,3.36,3.91],[1e-8,3.36,3.91]])
_LT = np.log(_GT[:,0]); _GR = _GT[:,1]; _GS = _GT[:,2]
def gsr(T): return float(np.interp(math.log(T) if T>0 else -50, _LT[::-1], _GR[::-1]))
def gss(T): return float(np.interp(math.log(T) if T>0 else -50, _LT[::-1], _GS[::-1]))

def omega_h2(m_chi, alpha_s):
    sv = math.pi*alpha_s**2/(4*m_chi**2)
    xf = 20.0
    for _ in range(50):
        gr = gsr(m_chi/xf)
        a = 0.0764*2*M_PL*m_chi*sv/math.sqrt(gr*xf)
        if a <= 0: break
        xn = math.log(a) - 0.5*math.log(xf)
        if abs(xn-xf) < 0.001: xf = xn; break
        xf = 0.5*xf + 0.5*xn
    gs = gss(m_chi/xf); gr = gsr(m_chi/xf)
    lv = math.sqrt(math.pi/45)*gs/math.sqrt(gr)*M_PL*m_chi
    Yi = xf/(lv*sv)
    return m_chi*Yi*S_0/RHO_CRIT_H2

# Load 17 BPs
import csv
CSV = os.path.join(os.path.dirname(__file__), "..",
    "Secluded-Majorana-SIDM", "predictions", "output", "sweep_17bp_results.csv")
bps = []
with open(CSV, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        bps.append({
            "label": row["BP"], "m_chi": float(row["m_chi_GeV"]),
            "m_phi": float(row["m_phi_MeV"])*1e-3,
            "alpha": float(row["alpha"]), "lam": float(row["lambda"]),
            "sm30": float(row["sigma_m_30"]),
        })

LAM_CRITS = [1.68, 6.45, 14.7, 26.0]

# ==============================================================
print("=" * 90)
print("  FREE θ SCAN — σ/m(30) and Ωh² as function of mixing angle")
print("=" * 90)
print()
print("  For each BP, α_total = α_original (fixed)")
print("  α_s(θ) = α_total × cos²θ   →  used for SIDM + relic")
print("  Scan θ from 0° to 45°")
print()

# JIT warmup
_ = sigma_T(20.0, 0.01, 1e-3, 30.0)
print("  [JIT done]")
print()

# Scan θ for selected BPs (representatives)
theta_deg = np.linspace(0, 45, 100)
theta_rad = np.radians(theta_deg)
THETA_RELIC = math.degrees(math.asin(1/3))  # 19.47°

selected_bps = [0, 2, 5, 6, 8, 13, 15]  # BP1,3,6,7,9,14,16

print("  Per-BP scan: θ → σ/m(30) [looking for resonance peaks]")
print("  " + "=" * 85)

for idx in selected_bps:
    bp = bps[idx]
    mc, mp, a0 = bp["m_chi"], bp["m_phi"], bp["alpha"]
    lam0 = a0 * mc / mp
    
    print(f"\n  {bp['label']}: m_χ={mc:.2f} GeV, λ_orig={lam0:.2f}")
    
    sms = []
    oms = []
    for th in theta_rad:
        alpha_s = a0 * math.cos(th)**2
        sm = sigma_T(mc, mp, alpha_s, 30.0)
        om = omega_h2(mc, alpha_s)
        sms.append(sm)
        oms.append(om)
    
    sms = np.array(sms)
    oms = np.array(oms)
    
    # Find peaks / resonances
    idx_max = np.argmax(sms)
    
    # Find θ where σ/m = 0.5 (SIDM threshold)
    above_05 = theta_deg[sms >= 0.5]
    theta_max_sidm = above_05[-1] if len(above_05) > 0 else 0
    
    # Find θ where Ωh² = 0.120
    cross_relic = None
    for i in range(len(oms)-1):
        if (oms[i] - 0.12) * (oms[i+1] - 0.12) < 0:
            # Linear interpolation
            t = (0.12 - oms[i]) / (oms[i+1] - oms[i])
            cross_relic = theta_deg[i] + t * (theta_deg[i+1] - theta_deg[i])
            break
    
    # Value at θ_relic
    th_r_idx = np.argmin(np.abs(theta_deg - THETA_RELIC))
    sm_at_relic = sms[th_r_idx]
    om_at_relic = oms[th_r_idx]
    
    # Print key values
    print(f"    θ=0°:      σ/m={sms[0]:.4f}  Ωh²={oms[0]:.4f}  (original)")
    print(f"    θ=19.47°:  σ/m={sm_at_relic:.4f}  Ωh²={om_at_relic:.4f}  (enforced)")
    print(f"    θ_max_σ:   {theta_deg[idx_max]:.1f}° → σ/m={sms[idx_max]:.4f}")
    print(f"    θ_max_SIDM: ≤{theta_max_sidm:.1f}° (σ/m ≥ 0.5)")
    if cross_relic:
        sm_cr = sigma_T(mc, mp, a0*math.cos(math.radians(cross_relic))**2, 30.0)
        print(f"    θ_relic:   {cross_relic:.1f}° (Ωh²=0.120) → σ/m={sm_cr:.4f}")
    else:
        print(f"    θ_relic:   no crossing (Ωh² never = 0.120 in this range)")
    
    # ASCII chart
    mx = max(sms)
    for i in range(0, len(theta_deg), 5):
        th, sm = theta_deg[i], sms[i]
        bar = "█" * int(sm / mx * 40)
        mark = ""
        if abs(th - THETA_RELIC) < 2.5: mark = " ◄ θ_relic"
        if th == 0: mark = " ◄ original"
        thresh = "│" if sm >= 0.5 else " "
        print(f"    {th:5.1f}° {bar}{thresh}{mark}")

# ==============================================================
#  KEY ANALYSIS: Where do resonances fall as function of θ?
# ==============================================================
print()
print("=" * 90)
print("  RESONANCE MAP: λ(θ) = λ_orig × cos²θ vs Yukawa resonances")
print("=" * 90)
print()
print(f"  Yukawa quasi-bound states at λ_crit = {LAM_CRITS}")
print()

print(f"  {'BP':>4}  {'λ_orig':>7}  {'λ(0°)':>7}  {'λ(19°)':>7}  {'λ(30°)':>7}  {'λ(45°)':>7}  "
      f"{'nearest resonance at θ=19°':>30}")
print("  " + "-" * 95)

for bp in bps:
    l0 = bp["lam"]
    l19 = l0 * math.cos(math.radians(19.47))**2
    l30 = l0 * math.cos(math.radians(30))**2
    l45 = l0 * math.cos(math.radians(45))**2
    
    # Find nearest resonance at θ=19°
    dists = [(abs(l19 - lc), f"R{i+1}({lc})") for i, lc in enumerate(LAM_CRITS)]
    dists.sort()
    near = dists[0]
    
    # Does any θ in [0,45] cross a resonance?
    crosses = []
    for lc in LAM_CRITS:
        if l45 <= lc <= l0:
            # cos²θ = lc/l0 → θ = arccos(sqrt(lc/l0))
            ratio = lc / l0
            if 0 < ratio <= 1:
                th_cross = math.degrees(math.acos(math.sqrt(ratio)))
                crosses.append((th_cross, lc))
    
    cross_str = ""
    if crosses:
        cross_str = "  CROSSES: " + ", ".join(f"R({lc:.2f})@{th:.1f}°" for th, lc in crosses)
    
    print(f"  {bp['label']:>4}  {l0:7.2f}  {l0:7.2f}  {l19:7.2f}  {l30:7.2f}  {l45:7.2f}  "
          f"dist={near[0]:.2f} to {near[1]}{cross_str}")

# ==============================================================
#  The real question: σ/m at resonance-crossing θ
# ==============================================================
print()
print("=" * 90)
print("  RESONANCE CROSSING: BPs that cross a Yukawa resonance as θ varies")
print("=" * 90)
print()

print(f"  {'BP':>4}  {'λ_orig':>7}  {'λ_crit':>7}  {'θ_cross':>8}  "
      f"{'σ/m(30)':>10}  {'Ωh²':>8}  {'SIDM?':>5}  {'relic?':>6}")
print("  " + "-" * 75)

for bp in bps:
    mc, mp, a0 = bp["m_chi"], bp["m_phi"], bp["alpha"]
    l0 = bp["lam"]
    
    for i_res, lc in enumerate(LAM_CRITS):
        ratio = lc / l0
        if 0.1 < ratio <= 1.0:
            th_cross = math.degrees(math.acos(math.sqrt(ratio)))
            if th_cross <= 45:
                alpha_s = a0 * ratio  # = a0 * cos²θ
                sm = sigma_T(mc, mp, alpha_s, 30.0)
                om = omega_h2(mc, alpha_s)
                sidm_ok = "✓" if 0.5 <= sm <= 10 else "✗"
                rel_ok = "✓" if 0.115 <= om <= 0.125 else "✗"
                print(f"  {bp['label']:>4}  {l0:7.2f}  R{i_res+1}={lc:.2f}  {th_cross:7.1f}°  "
                      f"{sm:10.4f}  {om:8.4f}  {sidm_ok:>5}  {rel_ok:>6}")

print()
print("=" * 90)
print("  BOTTOM LINE: Is there a θ that makes both SIDM AND relic work via resonance?")
print("=" * 90)
print()

# For each BP, find θ range where BOTH σ/m≥0.5 AND 0.115≤Ωh²≤0.125
for bp in bps:
    mc, mp, a0 = bp["m_chi"], bp["m_phi"], bp["alpha"]
    
    good_thetas = []
    for th in np.linspace(0, 40, 200):
        cs2 = math.cos(math.radians(th))**2
        a_s = a0 * cs2
        sm = sigma_T(mc, mp, a_s, 30.0)
        sm1k = sigma_T(mc, mp, a_s, 1000.0)
        om = omega_h2(mc, a_s)
        
        if sm >= 0.5 and sm1k < 0.47 and 0.115 <= om <= 0.125:
            good_thetas.append(th)
    
    if good_thetas:
        print(f"  {bp['label']}: θ ∈ [{min(good_thetas):.1f}°, {max(good_thetas):.1f}°]  "
              f"(width = {max(good_thetas)-min(good_thetas):.1f}°)")
    else:
        print(f"  {bp['label']}: NO θ simultaneously satisfies SIDM + relic (KT)")

print()
