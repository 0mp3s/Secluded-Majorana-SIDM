#!/usr/bin/env python3
"""
Condition 1: Amplitude |M|² for χχ → φφ  (Mixed Majorana SIDM)
================================================================

Model:  L ⊃ (1/2) χ̄(yₛ + i yₚ γ₅) χ φ     [Majorana χ, real scalar φ]

Process: χ(p₁) χ(p₂) → φ(k₁) φ(k₂)  via t-channel + u-channel

CRITICAL:  The relative sign is  M = M_t + M_u  (PLUS)
because the two diagrams differ by k₁↔k₂ (identical bosons = Bose symmetry).

Method:
  - Numerical 4×4 gamma matrices (Dirac representation, +−−− metric)
  - Spin-summed |M|² via trace formula
  - Gauss-Legendre integration over cos(θ) in CM frame
  - Threshold expansion:  σv = a₀ + a₁ v² + O(v⁴)

Key results:
  ★ a₀(yₛ, 0)  = 0     CP conserved → Majorana ¹S₀ (CP=−1) ↛ φφ (CP=+1)
  ★ a₀(0, yₚ)  = 0     CP conserved → same selection rule
  ★ a₀(yₛ, yₚ) = yₛ²yₚ²/(8π m_χ²)   CP VIOLATED → s-wave allowed!

      equivalently:  a₀ = 2π αₛ αₚ / m_χ²    where αᵢ = yᵢ²/(4π)

Date: 23 March 2026
"""

import sys, os
import numpy as np

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)

# ═══════════════════════════════════════════════════════════════
#  4×4 Gamma Matrices  (Dirac representation, metric +−−−)
# ═══════════════════════════════════════════════════════════════
I4 = np.eye(4, dtype=complex)

g0 = np.array([[ 1, 0, 0, 0],
                [ 0, 1, 0, 0],
                [ 0, 0,-1, 0],
                [ 0, 0, 0,-1]], dtype=complex)

g1 = np.array([[ 0, 0, 0, 1],
                [ 0, 0, 1, 0],
                [ 0,-1, 0, 0],
                [-1, 0, 0, 0]], dtype=complex)

g2 = np.array([[ 0,   0,  0, -1j],
                [ 0,   0, 1j,   0],
                [ 0,  1j,  0,   0],
                [-1j,  0,  0,   0]], dtype=complex)

g3 = np.array([[ 0, 0, 1, 0],
                [ 0, 0, 0,-1],
                [-1, 0, 0, 0],
                [ 0, 1, 0, 0]], dtype=complex)

g5 = 1j * g0 @ g1 @ g2 @ g3     # diag(−1,−1,+1,+1) in Dirac rep

# ── Self-test: Clifford algebra {γ^μ, γ^ν} = 2 η^{μν} ──
_eta = np.diag([1., -1., -1., -1.])
_gam = [g0, g1, g2, g3]
for _mu in range(4):
    for _nu in range(4):
        _anti = _gam[_mu] @ _gam[_nu] + _gam[_nu] @ _gam[_mu]
        assert np.allclose(_anti, 2 * _eta[_mu, _nu] * I4), \
            f"Clifford check failed ({_mu},{_nu})"
assert np.allclose(g5 @ g5, I4), "γ₅² ≠ I"
for _mu in range(4):
    assert np.allclose(g5 @ _gam[_mu] + _gam[_mu] @ g5, 0), \
        f"{{γ₅, γ^{_mu}}} ≠ 0"


# ═══════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════

def slash(p):
    """/p = p⁰γ⁰ − p¹γ¹ − p²γ² − p³γ³   (p stored as [E, px, py, pz])"""
    return p[0]*g0 - p[1]*g1 - p[2]*g2 - p[3]*g3

def dot4(a, b):
    """Minkowski inner product a·b = a⁰b⁰ − a⃗·b⃗"""
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]


# ═══════════════════════════════════════════════════════════════
#  Spin-summed |M|²  via trace technology
# ═══════════════════════════════════════════════════════════════

def M2_spin_summed(p1, p2, k1, k2, m_chi, ys, yp):
    """
    Compute  Σ_{spins} |M|²  for  χ(p₁) χ(p₂) → φ(k₁) φ(k₂).

    M = M_t + M_u   (PLUS sign: Bose symmetry for identical φ's)

    M_t = v̄(p₂) Γ [(q̸_t + m)/(t − m²)] Γ u(p₁)
    M_u = v̄(p₂) Γ [(q̸_u + m)/(u − m²)] Γ u(p₁)

    Γ = yₛ I₄ + i yₚ γ₅     (vertex factor)
    q_t = p₁ − k₁,   q_u = p₁ − k₂

    The amplitude is symmetric under k₁↔k₂ as required by Bose stats.

    Trace formula (using Γ̄ = Γ):
      Σ|M|² = T_tt/D_t² + T_uu/D_u² + 2 Re(T_tu)/(D_t D_u)

      T_xy = Tr[ (p̸₂ − m) Γ (q̸_x + m) Γ (p̸₁ + m) Γ (q̸_y + m) Γ ]
      D_x  = q_x² − m²
    """
    G = ys * I4 + 1j * yp * g5

    P1 = slash(p1) + m_chi * I4          # Σ_s u ū = /p₁ + m
    P2 = slash(p2) - m_chi * I4          # Σ_s v v̄ = /p₂ − m

    qt = p1 - k1
    qu = p1 - k2

    Dt = dot4(qt, qt) - m_chi**2
    Du = dot4(qu, qu) - m_chi**2

    Nt = slash(qt) + m_chi * I4
    Nu = slash(qu) + m_chi * I4

    def chain(Nx, Ny):
        """Tr[ P₂ Γ Nₓ Γ P₁ Γ Nᵧ Γ ]"""
        return np.trace(P2 @ G @ Nx @ G @ P1 @ G @ Ny @ G)

    Ttt = chain(Nt, Nt).real
    Tuu = chain(Nu, Nu).real
    Ttu = chain(Nt, Nu)

    # +2 Re(T_tu) because M = M_t + M_u (Bose symmetry)
    return Ttt / Dt**2 + Tuu / Du**2 + 2.0 * Ttu.real / (Dt * Du)


# ═══════════════════════════════════════════════════════════════
#  σ × v_rel  in the CM frame
# ═══════════════════════════════════════════════════════════════

def sigma_v_rel(v_cm, m_chi, m_phi, ys, yp, n_gl=80):
    """
    σ × v_Möller in CM frame.

    Formula:  σv = k_f / (512 π E³) × ∫₋₁¹ d(cosθ) Σ|M|²

    Includes 1/4 spin average and S=1/2 for identical final φ's.
    """
    if v_cm < 1e-15:
        v_cm = 1e-15

    gamma = 1.0 / np.sqrt(1.0 - v_cm**2)
    E     = m_chi * gamma
    p_mag = E * v_cm
    omega = E

    if omega < m_phi:
        return 0.0
    kf = np.sqrt(omega**2 - m_phi**2)

    p1 = np.array([E, 0., 0.,  p_mag])
    p2 = np.array([E, 0., 0., -p_mag])

    nodes, weights = np.polynomial.legendre.leggauss(n_gl)

    integral = 0.0
    for ct, w in zip(nodes, weights):
        st = np.sqrt(max(1.0 - ct*ct, 0.0))
        k1 = np.array([omega, 0., kf*st,  kf*ct])
        k2 = np.array([omega, 0.,-kf*st, -kf*ct])
        integral += w * M2_spin_summed(p1, p2, k1, k2, m_chi, ys, yp)

    return kf * integral / (512.0 * np.pi * E**3)


# ═══════════════════════════════════════════════════════════════
#  Threshold expansion:  σv = a₀ + a₁ v² + O(v⁴)
# ═══════════════════════════════════════════════════════════════

def extract_threshold(m_chi, m_phi, ys, yp, v_max=0.05, n_v=25):
    """Fit σv(v) = a₀ + a₁ v² and return (a₀, a₁, v_arr, sv_arr)."""
    v_arr  = np.linspace(1e-6, v_max, n_v)
    sv_arr = np.array([sigma_v_rel(v, m_chi, m_phi, ys, yp) for v in v_arr])

    A = np.column_stack([np.ones(n_v), v_arr**2])
    coeffs, *_ = np.linalg.lstsq(A, sv_arr, rcond=None)
    return coeffs[0], coeffs[1], v_arr, sv_arr


# ═══════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("  Condition 1: |M|² for χχ → φφ   (Mixed Majorana SIDM)")
    print("  M = M_t + M_u  (Bose symmetry for identical φ's)")
    print("=" * 72)

    m_chi = 20.69
    m_phi = 0.01134
    r     = m_phi / m_chi
    print(f"\n  m_χ = {m_chi} GeV,  m_φ = {m_phi*1e3:.2f} MeV,  r = {r:.4e}")

    # ──────────────────────────────────────────────────────────
    #  Part 1: a₀ for various couplings
    # ──────────────────────────────────────────────────────────
    print("\n" + "─" * 72)
    print("  Part 1:  a₀ for different coupling combinations")
    print("─" * 72 + "\n")

    combos = [
        (1.0, 0.0, "Pure scalar      "),
        (0.0, 1.0, "Pure pseudoscalar"),
        (1.0, 1.0, "Mixed (1, 1)     "),
        (2.0, 1.0, "Mixed (2, 1)     "),
        (1.0, 2.0, "Mixed (1, 2)     "),
        (2.0, 2.0, "Mixed (2, 2)     "),
        (3.0, 1.0, "Mixed (3, 1)     "),
        (1.0, 3.0, "Mixed (1, 3)     "),
        (0.5, 0.5, "Mixed (0.5, 0.5) "),
    ]

    print(f"  {'Label':<22s} {'yₛ':>5s} {'yₚ':>5s}  {'a₀ [GeV⁻²]':>14s}  {'a₁ [GeV⁻²]':>14s}")
    print("  " + "─" * 66)

    results = {}
    for ys, yp, lab in combos:
        a0, a1, _, _ = extract_threshold(m_chi, m_phi, ys, yp)
        results[(ys, yp)] = (a0, a1)
        print(f"  {lab}  {ys:5.1f} {yp:5.1f}  {a0:14.6e}  {a1:14.6e}")

    # ──────────────────────────────────────────────────────────
    #  Part 2: Verify a₀ ∝ yₛ² yₚ²
    # ──────────────────────────────────────────────────────────
    print("\n" + "─" * 72)
    print("  Part 2:  Verify  a₀ = c × yₛ² × yₚ²")
    print("─" * 72 + "\n")

    a0_ref = results[(1.0, 1.0)][0]
    print(f"  Reference: a₀(1,1) = c = {a0_ref:.6e} GeV⁻²\n")

    print(f"  {'(yₛ,yₚ)':<12s}  {'yₛ²yₚ²':>8s}  {'a₀ num':>14s}  "
          f"{'c·yₛ²yₚ²':>14s}  {'ratio':>8s}  {'ok':>3s}")
    print("  " + "─" * 66)

    all_pass = True
    for ys, yp, _ in combos:
        a0_num  = results[(ys, yp)][0]
        ysyp    = ys**2 * yp**2
        a0_pred = a0_ref * ysyp

        if abs(a0_pred) > 1e-30:
            ratio = a0_num / a0_pred
            ok    = abs(ratio - 1.0) < 0.02
        elif abs(a0_num) < abs(a0_ref) * 1e-6:
            ratio = 0.0
            ok    = True    # Both zero = OK
        else:
            ratio = float('inf')
            ok    = False

        if not ok:
            all_pass = False
        mark = "✓" if ok else "✗"
        print(f"  ({ys:.1f},{yp:.1f}){'':<6s}  {ysyp:8.2f}  "
              f"{a0_num:14.6e}  {a0_pred:14.6e}  {ratio:8.5f}  {mark:>3s}")

    scaling_status = "PASS" if all_pass else "FAIL"
    print(f"\n  ★ Scaling test  a₀ ∝ yₛ² yₚ² :  {scaling_status}")

    # ──────────────────────────────────────────────────────────
    #  Part 3: Analytical comparison
    # ──────────────────────────────────────────────────────────
    print("\n" + "─" * 72)
    print("  Part 3:  Analytical comparison  (m_φ ≪ m_χ limit)")
    print("─" * 72 + "\n")

    # Analytical: a₀ = yₛ²yₚ² / (8π m_χ²) = 2π αₛ αₚ / m_χ²
    beta_f = np.sqrt(1.0 - r**2)
    c_exact = beta_f * m_chi / (2.0 * np.pi * (2*m_chi**2 - m_phi**2)**2)
    c_lo    = 1.0 / (8.0 * np.pi * m_chi**2)

    dev_exact = abs(a0_ref / c_exact - 1.0) if c_exact != 0 else float('inf')
    dev_lo    = abs(a0_ref / c_lo - 1.0)

    print(f"  c_exact = β_f m/(2π D²)    = {c_exact:.6e} GeV⁻²   |dev| = {dev_exact:.2e}")
    print(f"  c_LO    = 1/(8π m_χ²)      = {c_lo:.6e} GeV⁻²   |dev| = {dev_lo:.2e}")
    print(f"  c_num                       = {a0_ref:.6e} GeV⁻²")

    best_dev = min(dev_exact, dev_lo)
    if best_dev < 0.01:
        print(f"\n  ✓ Analytical match at {best_dev*100:.5f}% level")
    else:
        print(f"\n  ⚠ Deviation = {best_dev*100:.3f}% (investigating...)")

    # ──────────────────────────────────────────────────────────
    #  Part 4: Physical benchmark  (relic density target)
    # ──────────────────────────────────────────────────────────
    print("\n" + "─" * 72)
    print("  Part 4:  Physical benchmark  (relic density)")
    print("─" * 72 + "\n")

    alpha_tot = 1.048e-3
    alpha_s   = alpha_tot / 2.0
    alpha_p   = alpha_tot / 2.0
    ys_phys   = np.sqrt(4 * np.pi * alpha_s)
    yp_phys   = np.sqrt(4 * np.pi * alpha_p)

    a0_phys, a1_phys, _, _ = extract_threshold(m_chi, m_phi, ys_phys, yp_phys)
    a0_formula = 2.0 * np.pi * alpha_s * alpha_p / m_chi**2
    sv_canon   = 2.6e-9     # GeV⁻²  (≈ 3 × 10⁻²⁶ cm³/s)

    print(f"  α_total = {alpha_tot:.3e}  →  αₛ = αₚ = {alpha_s:.3e}")
    print(f"  yₛ = {ys_phys:.6f},  yₚ = {yp_phys:.6f}")
    print(f"\n  a₀ (numerical)       = {a0_phys:.4e} GeV⁻²")
    print(f"  a₀ (2π αₛ αₚ / m_χ²) = {a0_formula:.4e} GeV⁻²")
    print(f"  Canonical relic       = {sv_canon:.4e} GeV⁻²  (≈ 3×10⁻²⁶ cm³/s)")
    print(f"\n  num / 2π αₛαₚ/m²  = {a0_phys / a0_formula:.6f}")
    print(f"  a₀ / canonical     = {a0_phys / sv_canon:.4f}")

    # ──────────────────────────────────────────────────────────
    #  Part 5: σv(v) velocity profile
    # ──────────────────────────────────────────────────────────
    print("\n" + "─" * 72)
    print("  Part 5:  σv(v) velocity dependence  [mixed, physical couplings]")
    print("─" * 72 + "\n")

    v_list = [0.001, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3]
    print(f"  {'v':>8s}  {'σv [GeV⁻²]':>14s}  {'σv / a₀':>10s}")
    print("  " + "─" * 36)
    for v in v_list:
        sv = sigma_v_rel(v, m_chi, m_phi, ys_phys, yp_phys)
        r_sv = sv / a0_phys if abs(a0_phys) > 0 else 0
        print(f"  {v:8.3f}  {sv:14.6e}  {r_sv:10.6f}")

    # ──────────────────────────────────────────────────────────
    #  Part 6: Direct threshold (v=0) cross check
    # ──────────────────────────────────────────────────────────
    print("\n" + "─" * 72)
    print("  Part 6:  Direct threshold evaluation (v = 0 exactly)")
    print("─" * 72 + "\n")

    k0 = np.sqrt(m_chi**2 - m_phi**2)
    p1_thr = np.array([m_chi, 0., 0., 0.])
    p2_thr = np.array([m_chi, 0., 0., 0.])

    nodes, weights = np.polynomial.legendre.leggauss(80)

    # Physical couplings
    integral_phys = 0.0
    for ct, w in zip(nodes, weights):
        st = np.sqrt(max(1 - ct*ct, 0.0))
        k1 = np.array([m_chi, 0., k0*st,  k0*ct])
        k2 = np.array([m_chi, 0.,-k0*st, -k0*ct])
        integral_phys += w * M2_spin_summed(p1_thr, p2_thr, k1, k2,
                                            m_chi, ys_phys, yp_phys)
    a0_direct = k0 * integral_phys / (512.0 * np.pi * m_chi**3)
    print(f"  a₀ (direct, v=0)   = {a0_direct:.6e} GeV⁻²")
    print(f"  a₀ (from fit)      = {a0_phys:.6e} GeV⁻²")
    dev_fit = abs(a0_direct / a0_phys - 1)*100 if abs(a0_phys) > 0 else 0
    print(f"  agreement          = {dev_fit:.4f}%")

    # Pure scalar at threshold (must be zero)
    integral_s = 0.0
    for ct, w in zip(nodes, weights):
        st = np.sqrt(max(1 - ct*ct, 0.0))
        k1 = np.array([m_chi, 0., k0*st,  k0*ct])
        k2 = np.array([m_chi, 0.,-k0*st, -k0*ct])
        integral_s += w * M2_spin_summed(p1_thr, p2_thr, k1, k2,
                                         m_chi, 1.0, 0.0)
    a0_s_direct = k0 * integral_s / (512.0 * np.pi * m_chi**3)
    print(f"\n  a₀(1,0) direct     = {a0_s_direct:.6e}  (should be ~0)")

    # (1,1) at threshold
    integral_11 = 0.0
    for ct, w in zip(nodes, weights):
        st = np.sqrt(max(1 - ct*ct, 0.0))
        k1 = np.array([m_chi, 0., k0*st,  k0*ct])
        k2 = np.array([m_chi, 0.,-k0*st, -k0*ct])
        integral_11 += w * M2_spin_summed(p1_thr, p2_thr, k1, k2,
                                          m_chi, 1.0, 1.0)
    a0_11_direct = k0 * integral_11 / (512.0 * np.pi * m_chi**3)
    print(f"  a₀(1,1) direct     = {a0_11_direct:.6e} GeV⁻²")
    print(f"  a₀(1,1) from fit   = {a0_ref:.6e} GeV⁻²")

    # ──────────────────────────────────────────────────────────
    #  Summary
    # ──────────────────────────────────────────────────────────
    a0_s = results[(1., 0.)][0]
    a0_p = results[(0., 1.)][0]

    zero_s = abs(a0_s) < abs(a0_ref) * 1e-6
    zero_p = abs(a0_p) < abs(a0_ref) * 1e-6
    nonzero_mix = abs(a0_ref) > 1e-30

    condition1 = zero_s and zero_p and nonzero_mix and all_pass

    print("\n" + "=" * 72)
    print("  SUMMARY — CONDITION 1")
    print("=" * 72)
    print(f"""
  Pure scalar (yₛ, 0):        a₀ = {a0_s: .2e}  {"→ ZERO ✓" if zero_s else "→ NONZERO ✗"}
  Pure pseudoscalar (0, yₚ):   a₀ = {a0_p: .2e}  {"→ ZERO ✓" if zero_p else "→ NONZERO ✗"}
  Mixed (yₛ, yₚ):             a₀ = {a0_ref: .6e}  {"→ NONZERO ✓" if nonzero_mix else "→ ZERO ✗"}

  Scaling:  a₀ = c × yₛ² yₚ²   test: {scaling_status}

  Analytical formula:
    a₀ = yₛ² yₚ² / (8π m_χ²) = 2π αₛ αₚ / m_χ²

  CP physics:
    Majorana ¹S₀ → CP = −1.   Two identical scalars s-wave → CP = +1.
    → s-wave REQUIRES CP violation.
    → Only mixed coupling (yₛ ≠ 0 AND yₚ ≠ 0) breaks CP.
    → a₀ ∝ yₛ² yₚ² is the CP-violating cross-term.

  ★ CONDITION 1:  {"PASSED ✓" if condition1 else "NEEDS REVIEW ✗"} ★
""")
    return condition1


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    try:
        import sys as _sys, os as _os
        _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'core'))
        from tg_notify import notify
        notify("\u2705 condition1_amplitude done!")
    except Exception:
        pass
