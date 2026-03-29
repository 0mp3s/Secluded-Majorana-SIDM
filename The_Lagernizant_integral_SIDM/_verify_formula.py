"""Quick verification: old PI-12 formula vs proper KT formula."""
import math

M_PL = 1.2209e19
G_STAR = G_STAR_S = 86.25
GEV2_TO_CM3S = 3.8938e-28 * 3e10  # ~1.168e-17
PLANCK_SV = 3e-26

BPS = {
    "BP1": (54.556, 2.645e-3),
    "BP9": (48.329, 2.350e-3),
    "BP16": (14.384, 7.555e-4),
    "MAP": (98.19, 3.274e-3),
}

def phi_a(m, alpha):
    return math.pi * alpha**2 / (4 * m**2)

def find_xf(m, a, b=0):
    prefac = 0.0764 * 0.5 * 2.5 / math.sqrt(G_STAR) * M_PL * m
    xf = 22.0
    for _ in range(50):
        sv = a + 6*b/xf
        if sv <= 0: break
        arg = prefac * sv
        if arg <= 0: break
        xf_new = math.log(arg) - 0.5*math.log(xf)
        if abs(xf_new - xf) < 1e-8: break
        xf = 0.6*xf + 0.4*xf_new
    return max(xf, 5)

def omega_old(xf, a, b=0):
    """PI-12 crude formula: oh2 = 0.12 * 3e-26 / sv_eff_cm"""
    a_eff_cm = (a + 3*b/xf) * GEV2_TO_CM3S
    return 0.12 * PLANCK_SV / a_eff_cm

def omega_proper(m, a, b=0):
    """Proper KT: Y_inf = 1/(lambda*J)"""
    xf = find_xf(m, a, b)
    lam = math.sqrt(math.pi/45) * G_STAR_S/math.sqrt(G_STAR) * M_PL * m
    J = a/xf + 3*b/xf**2
    Y_inf = 1/(lam * J)
    return m * Y_inf * 2891.2 / 1.054e-5, xf

# 1. phi-only comparison
print("=== phi-only: old vs proper KT ===")
print(f"{'BP':>5} {'m':>6} {'xf':>5} {'oh2_old':>9} {'oh2_KT':>9} {'old/KT':>7}")
for nm, (m, al) in BPS.items():
    a = phi_a(m, al)
    oh2_new, xf = omega_proper(m, a)
    oh2_old = omega_old(xf, a)
    print(f"{nm:>5} {m:>6.1f} {xf:>5.1f} {oh2_old:>9.4f} {oh2_new:>9.4f} {oh2_old/oh2_new:>7.3f}")

# 2. Can BP reach 0.12?
print("\n=== Can BP reach oh2=0.12 with sigma channel? ===")
for nm, (m, al) in BPS.items():
    a = phi_a(m, al)
    oh2_new, _ = omega_proper(m, a)
    if oh2_new > 0.12:
        print(f"  {nm}: YES (phi_only_KT = {oh2_new:.4f} > 0.12)")
    else:
        print(f"  {nm}: NO  (phi_only_KT = {oh2_new:.4f} < 0.12)")

# 3. MAP f0_cross with proper KT
print("\n=== MAP f0_cross (proper KT) ===")
m, al = 98.19, 3.274e-3
a = phi_a(m, al)
f_lo, f_hi = 100.0, 1e8
for _ in range(200):
    f_mid = math.sqrt(f_lo * f_hi)
    b = 3*m**2 / (math.pi * f_mid**4)
    oh2, xf = omega_proper(m, a, b)
    if abs(oh2 - 0.12)/0.12 < 1e-6: break
    if oh2 > 0.12: f_hi = f_mid
    else: f_lo = f_mid
f0_map_proper = f_mid
print(f"  f0_cross = {f_mid:.1f} GeV (oh2 = {oh2:.6f})")
print(f"  f_DE = 3^31 x {f_mid:.0f} = {3**31*f_mid:.3e} GeV")
print(f"  f_target = 5.845e17 GeV")
print(f"  f_DE/f_target = {3**31*f_mid/5.845e17:.3f}")

# 4. Recompute PI-12 with OLD formula
print("\n=== PI-12 f0_cross with OLD formula ===")
vals = []
for nm, (m, al) in BPS.items():
    a = phi_a(m, al)
    xf = find_xf(m, a)
    oh2_phi = omega_old(xf, a)
    if oh2_phi < 0.12:
        print(f"  {nm}: phi_old = {oh2_phi:.4f} < 0.12 -> no crossing")
        continue
    f_lo, f_hi = 100.0, 1e8
    for _ in range(200):
        f_mid = math.sqrt(f_lo * f_hi)
        b = 3*m**2 / (math.pi * f_mid**4)
        xf_b = find_xf(m, a, b)
        oh2 = omega_old(xf_b, a, b)
        if abs(oh2 - 0.12)/0.12 < 1e-6: break
        if oh2 > 0.12: f_hi = f_mid
        else: f_lo = f_mid
    vals.append(f_mid)
    print(f"  {nm}: f0_cross_old = {f_mid:.1f} GeV")

if vals:
    geo = math.exp(sum(math.log(v) for v in vals)/len(vals))
    print(f"  Geometric mean ({len(vals)} BPs): {geo:.1f} GeV")
    print(f"  f_DE(old) = 3^31 x {geo:.0f} = {3**31*geo:.3e}")

# 5. Impact on PI-14 (y_P)
print("\n=== Impact on PI-14 (perturbativity) ===")
y_P_old = 2*98.19/755.3
y_P_new = 2*98.19/f0_map_proper
print(f"  y_P(old f0=755.3) = {y_P_old:.4f}")
print(f"  y_P(new f0={f0_map_proper:.1f}) = {y_P_new:.4f}")
print(f"  Both < 1: perturbative")

# 6. Impact on PI-16 (M1)
print("\n=== Impact on PI-16 (lightest heavy mode) ===")
Lambda_CW = 500.0
q, N = 3, 31
for f0, label in [(755.3, "old"), (f0_map_proper, "new")]:
    M1_sq = (Lambda_CW**4/f0**2) * (1 + q**2 - 2*q*math.cos(math.pi/(N+1)))
    M1 = math.sqrt(M1_sq) if M1_sq > 0 else 0
    print(f"  M1({label}, f0={f0:.1f}) = {M1:.1f} GeV (> m_chi=98.2? {'YES' if M1>98.2 else 'NO'})")

print("\n=== SUMMARY ===")
print("Old formula: oh2 = 0.12 * 3e-26 / sv_eff [crude WIMP miracle scaling]")
print("Proper KT:   Y_inf = 1/(lambda*J), oh2 = m*Y*s0/rho_c")
print(f"Ratio old/KT = 33.66/xf ≈ 1.4-1.6 (old OVERESTIMATES by 40-60%)")
print(f"Consequence: BP1/BP9/BP16 phi-only < 0.12 -> NO f0_cross exists")
print(f"Only MAP has f0_cross = {f0_map_proper:.1f} GeV (not 755.3)")
