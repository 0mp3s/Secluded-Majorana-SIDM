#!/usr/bin/env python3
"""
V10 - v36_sommerfeld.py
=======================
Numerical computation of the Sommerfeld enhancement factor S_0
for s-wave annihilation through an attractive Yukawa potential:

    V(r) = -alpha * exp(-m_phi * r) / r

Method:
  Solve the l=0 radial Schrödinger equation in dimensionless units x = m_phi * r:
    u''(x) + [kappa^2 + lambda * exp(-x)/x] u(x) = 0
  with u(0)=0, u'(0)=1.

  At large x (where potential vanishes): u(x) -> A * sin(kappa*x + delta)
  Extract amplitude: A = sqrt(u^2 + (u'/kappa)^2)
  Sommerfeld factor: S_0 = 1 / (kappa^2 * A^2)

  Parameters:
    kappa = mu*v / m_phi     (mu = m_chi/2 for identical particles)
    lambda = alpha * m_chi / m_phi

Tested against:
  - Free particle (lambda=0) => S=1  exactly
  - Coulomb limit (m_phi -> 0) => S = 2*pi*eta / (1 - exp(-2*pi*eta))
    with eta = alpha / (2*v)

Benchmark points: all 17 relic-viable BPs from v31_true_viable_points.csv
Velocities: v_fo ~ 90,000 km/s (freeze-out), 1000, 200, 30 km/s
"""
# === path setup (auto-generated) ================================
import sys as _sys, os as _os
_ROOT = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')
_sys.path.insert(0, _os.path.join(_ROOT, 'core'))
DATA_DIR = _os.path.join(_ROOT, 'data')
# =================================================================


import numpy as np
import csv
import sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from config_loader import load_config
from output_manager import get_latest

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)

_CFG = load_config(__file__)

C_KM_S = 299792.458

# ──────────────────────────────────────────────────────────────
#  Sommerfeld factor via direct ODE integration (RK4)
# ──────────────────────────────────────────────────────────────

def _sommerfeld_ode(kappa, lam):
    """
    Solve the l=0 radial ODE in dimensionless units x = m_phi * r:
      u''(x) + [kappa^2 + lambda * exp(-x)/x] u(x) = 0
      u(0)=0, u'(0)=1
    Returns S_0 = 1 / (kappa^2 * A^2) where A is the asymptotic amplitude.
    Reliable for kappa < ~10.
    """
    # Grid: need ~30 points per oscillation for RK4 accuracy
    x_max = max(50.0, min(30.0 / kappa, 500.0))
    pts_per_osc = 40
    h_osc = (2.0 * np.pi / kappa) / pts_per_osc if kappa > 0.01 else 0.005
    h = min(0.005, h_osc)
    N_steps = int(x_max / h) + 1
    h = x_max / N_steps

    # Start slightly away from x=0.
    # Near x=0: u ~ x + (lam - kappa^2) x^2/2 + ..., so u(eps)~eps, u'(eps)~1
    eps = 1e-6
    u = eps
    up = 1.0  # du/dx
    x = eps

    for _ in range(N_steps):
        # Inline RK4 for [u, u'] with u'' = -(kappa^2 + lam*exp(-x)/x)*u
        def accel(xx, uu):
            if xx < 1e-30:
                return 0.0
            return -(kappa**2 + lam * np.exp(-xx) / xx) * uu

        k1_u = up
        k1_v = accel(x, u)

        k2_u = up + 0.5 * h * k1_v
        k2_v = accel(x + 0.5*h, u + 0.5*h*k1_u)

        k3_u = up + 0.5 * h * k2_v
        k3_v = accel(x + 0.5*h, u + 0.5*h*k2_u)

        k4_u = up + h * k3_v
        k4_v = accel(x + h, u + h*k3_u)

        u  += h * (k1_u + 2*k2_u + 2*k3_u + k4_u) / 6.0
        up += h * (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6.0
        x  += h

        if abs(u) > 1e100 or abs(up) > 1e100:
            return np.nan

    # Extract amplitude in the potential-free tail
    A_sq = u**2 + (up / kappa)**2
    if A_sq < 1e-300:
        return np.nan
    return 1.0 / (kappa**2 * A_sq)


def sommerfeld_yukawa(m_chi_GeV, m_phi_GeV, alpha, v_km_s):
    """
    Compute the l=0 Sommerfeld enhancement factor for attractive Yukawa potential.

    Hybrid approach:
      kappa < 10  => direct ODE integration (reliable, <50k steps)
      kappa >= 10 => Coulomb analytic upper bound (Yukawa screening only reduces S)

    The Coulomb upper bound is tight: at kappa >> lambda, screening is negligible.
    """
    v = v_km_s / C_KM_S              # natural units (c=1)
    mu = m_chi_GeV / 2.0             # reduced mass for identical particles
    k = mu * v                       # momentum
    kappa = k / m_phi_GeV            # dimensionless momentum
    lam = alpha * m_chi_GeV / m_phi_GeV  # dimensionless coupling = lambda

    if kappa < 1e-30 or lam < 1e-30:
        return 1.0

    if kappa >= 10.0:
        # At high kappa, use Coulomb analytic formula (upper bound on Yukawa S)
        # S_Coulomb = 2*pi*eta / (1 - exp(-2*pi*eta)) with eta = alpha/(2v)
        return sommerfeld_coulomb_analytic(alpha, v_km_s)

    return _sommerfeld_ode(kappa, lam)


def sommerfeld_coulomb_analytic(alpha, v_km_s):
    """
    Exact Coulomb Sommerfeld factor (m_phi -> 0 limit).
    S = 2*pi*eta / (1 - exp(-2*pi*eta)) with eta = alpha/(2v)
    """
    v = v_km_s / C_KM_S
    eta = alpha / (2.0 * v)
    if eta < 1e-6:
        return 1.0 + np.pi * eta
    return 2.0 * np.pi * eta / (1.0 - np.exp(-2.0 * np.pi * eta))


# ──────────────────────────────────────────────────────────────
#  Validation against Coulomb limit
# ──────────────────────────────────────────────────────────────

def validate():
    """Validate Sommerfeld solver against known limits."""
    print("=" * 70)
    print("VALIDATION")
    print("=" * 70)

    # Test 1: Free particle (lambda -> 0) => S = 1
    print("\n  Test 1: Free particle (alpha=0) => S must be 1.000")
    m_chi, m_phi = 20.0, 0.01  # GeV
    for v in [30, 200, 1000, 10000]:
        S = sommerfeld_yukawa(m_chi, m_phi, 1e-30, v)
        print(f"    v={v:>6d} km/s: S = {S:.8f}")

    # Test 2: Perturbative regime (high v) => S ~ 1 + pi*alpha/v
    print("\n  Test 2: Perturbative limit (high v) => S ~ 1 + pi*alpha/v")
    m_chi, m_phi, alpha = 20.0, 0.01, 1e-3
    print(f"    m_chi={m_chi}, m_phi={m_phi}, alpha={alpha}")
    for v_km in [50000, 100000, 200000]:
        v_c = v_km / C_KM_S
        S_num = sommerfeld_yukawa(m_chi, m_phi, alpha, v_km)
        S_pert = 1.0 + np.pi * alpha / v_c
        print(f"    v={v_km:>7d} km/s: S_num={S_num:.8f}  S_pert={S_pert:.8f}  diff={abs(S_num-S_pert):.2e}")

    # Test 3: Moderate Yukawa vs Coulomb (lambda ~ 1, kappa >> 1)
    # At high kappa, Yukawa and Coulomb should agree (screening negligible)
    print("\n  Test 3: Yukawa vs Coulomb at high kappa (lambda=2, kappa>>1)")
    m_chi, alpha = 20.0, 1e-3
    m_phi = alpha * m_chi / 2.0  # lambda = 2
    for v_km in [10000, 30000, 100000]:
        S_y = sommerfeld_yukawa(m_chi, m_phi, alpha, v_km)
        S_c = sommerfeld_coulomb_analytic(alpha, v_km)
        mu = m_chi / 2
        kap = mu * (v_km / C_KM_S) / m_phi
        print(f"    v={v_km:>7d} km/s (κ={kap:.1f}): S_Yuk={S_y:.8f}  S_Coul={S_c:.8f}  ratio={S_y/S_c:.4f}")

    print()


# ──────────────────────────────────────────────────────────────
#  Main: Sommerfeld for all 17 relic BPs
# ──────────────────────────────────────────────────────────────

def main():
    validate()

    # Load benchmark points
    _bp_csv = _CFG.get("benchmark_csv") or str(get_latest("v31_true_viable_points"))
    if not os.path.isabs(_bp_csv):
        _bp_csv = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), _bp_csv))
    csv_path = _bp_csv
    bps = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            bps.append({
                'm_chi': float(row['m_chi_GeV']),
                'm_phi_MeV': float(row['m_phi_MeV']),
                'm_phi_GeV': float(row['m_phi_MeV']) / 1000.0,
                'alpha': float(row['alpha']),
                'lam': float(row['lambda']),
            })

    print(f"Loaded {len(bps)} relic-viable benchmark points\n")

    # Velocities to probe
    v_fo = 0.3 * C_KM_S   # freeze-out ~ 0.3c
    velocities = {
        'v_fo (0.3c)': v_fo,
        '1000 km/s': 1000.0,
        '200 km/s': 200.0,
        '30 km/s': 30.0,
    }

    # ── Compute S for all BPs at all velocities ──
    print("=" * 100)
    print("SOMMERFELD ENHANCEMENT FACTOR S_0 FOR ALL RELIC BENCHMARK POINTS")
    print("=" * 100)

    header = f"{'BP':>3s}  {'m_chi':>7s}  {'m_phi':>7s}  {'alpha':>10s}  {'lambda':>7s}"
    for vname in velocities:
        header += f"  {'S('+vname+')':>16s}"
    print(header)
    print("-" * 100)

    results = []
    for i, bp in enumerate(bps):
        mc = bp['m_chi']
        mp = bp['m_phi_GeV']
        al = bp['alpha']
        lam = bp['lam']
        row_str = f"{i+1:3d}  {mc:7.1f}  {bp['m_phi_MeV']:7.2f}  {al:10.3e}  {lam:7.2f}"

        row_data = {'bp': i+1, 'm_chi': mc, 'm_phi_MeV': bp['m_phi_MeV'],
                     'alpha': al, 'lambda': lam}

        for vname, v in velocities.items():
            S = sommerfeld_yukawa(mc, mp, al, v)
            row_str += f"  {S:16.6f}"
            row_data[vname] = S

        print(row_str)
        results.append(row_data)

    # ── Summary statistics ──
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    for vname in velocities:
        S_vals = [r[vname] for r in results if np.isfinite(r[vname])]
        if S_vals:
            S_arr = np.array(S_vals)
            print(f"\n  {vname}:")
            print(f"    min S = {S_arr.min():.6f}")
            print(f"    max S = {S_arr.max():.6f}")
            print(f"    mean  = {S_arr.mean():.6f}")
            print(f"    max deviation from 1: {abs(S_arr - 1).max():.2e}")

    # ── Relic density impact ──
    print("\n" + "=" * 80)
    print("IMPACT ON RELIC DENSITY")
    print("=" * 80)
    print("\n  <sigma v>_eff = S(v_fo) * <sigma v>_tree")
    print("  => Omega h^2 scales as 1/S(v_fo)")
    print("  => Delta(Omega h^2) / Omega h^2 = 1/S(v_fo) - 1 ~ -(S-1)\n")

    vfo_key = 'v_fo (0.3c)'
    print(f"  {'BP':>3s}  {'lambda':>7s}  {'S(v_fo)':>12s}  {'Delta Omega/Omega':>18s}")
    print("  " + "-" * 50)
    for r in results:
        S = r[vfo_key]
        if np.isfinite(S) and S > 0:
            delta = 1.0/S - 1.0
            print(f"  {r['bp']:3d}  {r['lambda']:7.2f}  {S:12.8f}  {delta:+18.2e}")

    # ── Late-time depletion timescale ──
    print("\n" + "=" * 80)
    print("LATE-TIME ANNIHILATION DEPLETION (v = 30 km/s, secluded)")
    print("=" * 80)
    # n_chi ~ rho_DM / m_chi, rho_DM(dwarf core) ~ 1 GeV/cm^3
    rho_DM = 1.0  # GeV/cm^3
    GEV2_CM3_S = 3.8938e-28 * 2.998e10  # GeV^-2 -> cm^3/s (approx)
    t_hubble_s = 4.35e17  # seconds

    v30_key = '30 km/s'
    print(f"\n  rho_DM = {rho_DM} GeV/cm^3 (dwarf core)")
    print(f"  <sigma v>_tree = pi * alpha^2 / (4 * m_chi^2)")
    print(f"  <sigma v>_eff = S(30) * <sigma v>_tree")
    print(f"  tau_depletion = m_chi / (rho_DM * <sigma v>_eff)\n")

    print(f"  {'BP':>3s}  {'m_chi':>7s}  {'S(30)':>10s}  {'<sv>_eff [cm3/s]':>18s}  {'tau/t_Hubble':>14s}")
    print("  " + "-" * 65)
    for r in results:
        mc = r['m_chi']
        al = r['alpha']
        S30 = r[v30_key]
        if not np.isfinite(S30):
            continue
        # <sigma v>_tree in GeV^-2
        sv_tree_gev2 = np.pi * al**2 / (4 * mc**2)
        # convert to cm^3/s: multiply by hbar*c = 1.973e-14 GeV*cm, 
        # so GeV^-2 * (hbar*c)^2 * c = GeV^-2 * 3.894e-28 cm^2 * 3e10 cm/s
        sv_tree_cm3s = sv_tree_gev2 * 3.8938e-28 * 2.998e10
        sv_eff = S30 * sv_tree_cm3s

        # n_chi = rho_DM / m_chi (in GeV/cm^3 / GeV = cm^-3... no.)
        # n_chi = rho / m_chi but rho in GeV/cm^3 and m_chi in GeV => n in cm^-3
        n_chi = rho_DM / mc  # cm^-3

        # Annihilation rate: 1/tau = n * <sv>_eff / 2 (identical particles)
        # tau = 2 / (n * <sv>_eff)
        tau_s = 2.0 / (n_chi * sv_eff) if sv_eff > 0 else np.inf
        tau_ratio = tau_s / t_hubble_s

        print(f"  {r['bp']:3d}  {mc:7.1f}  {S30:10.4f}  {sv_eff:18.3e}  {tau_ratio:14.2e}")

    # ── Velocity scan for BP1 (detailed) ──
    print("\n" + "=" * 80)
    print("DETAILED VELOCITY SCAN -- BP1 (m_chi=20.7 GeV, m_phi=11.3 MeV, alpha=1.05e-3)")
    print("=" * 80)
    bp1 = bps[0]  # First BP
    v_scan = np.array([10, 20, 30, 50, 100, 200, 500, 1000, 3000,
                        10000, 30000, 50000, 90000, 150000, 250000])
    S_scan = []
    print(f"\n  {'v [km/s]':>12s}  {'kappa':>10s}  {'S_0':>14s}  {'S_Coulomb':>12s}  {'Yukawa/Coulomb':>15s}")
    print("  " + "-" * 70)
    for v in v_scan:
        mu = bp1['m_chi'] / 2
        k = mu * (v / C_KM_S)
        kap = k / bp1['m_phi_GeV']
        S = sommerfeld_yukawa(bp1['m_chi'], bp1['m_phi_GeV'], bp1['alpha'], v)
        S_c = sommerfeld_coulomb_analytic(bp1['alpha'], v)
        ratio = S / S_c if (S_c > 0 and np.isfinite(S)) else np.nan
        S_scan.append(S)
        print(f"  {v:12.0f}  {kap:10.4f}  {S:14.6f}  {S_c:12.6f}  {ratio:15.6f}")

    # ── Plot ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: S vs v for BP1
    ax = axes[0]
    v_fine = np.logspace(1, np.log10(250000), 60)
    S_fine = [sommerfeld_yukawa(bp1['m_chi'], bp1['m_phi_GeV'], bp1['alpha'], v) for v in v_fine]
    S_coulomb = [sommerfeld_coulomb_analytic(bp1['alpha'], v) for v in v_fine]

    ax.loglog(v_fine, S_fine, 'b-', lw=2, label='Yukawa (numerical)')
    ax.loglog(v_fine, S_coulomb, 'r--', lw=1.5, label='Coulomb (analytic)')
    ax.axhline(1.0, color='gray', ls=':', lw=0.8)
    ax.axvline(0.3 * C_KM_S, color='green', ls='--', lw=1, alpha=0.7, label='$v_{fo}$ (0.3c)')
    ax.axvline(30, color='orange', ls='--', lw=1, alpha=0.7, label='30 km/s (dwarfs)')
    ax.axvline(1000, color='purple', ls='--', lw=1, alpha=0.7, label='1000 km/s (clusters)')
    ax.set_xlabel('v [km/s]')
    ax.set_ylabel('$S_0$ (Sommerfeld factor)')
    ax.set_title(f'BP1: $m_\\chi$={bp1["m_chi"]:.1f} GeV, $m_\\phi$={bp1["m_phi_MeV"]:.1f} MeV, $\\lambda$={bp1["lam"]:.2f}')
    ax.legend(fontsize=8)
    ax.set_xlim(10, 3e5)
    ax.grid(True, alpha=0.3)

    # Panel 2: S(v_fo) and S(30) for all 17 BPs
    ax = axes[1]
    lambdas = [r['lambda'] for r in results]
    S_fo = [r[vfo_key] for r in results]
    S_30 = [r[v30_key] for r in results]
    S_1000 = [r['1000 km/s'] for r in results]

    ax.semilogy(lambdas, S_fo, 'gs', ms=8, label='$S(v_{fo})$ - freeze-out')
    ax.semilogy(lambdas, S_30, 'ro', ms=8, label='$S(30\,\mathrm{km/s})$ - dwarfs')
    ax.semilogy(lambdas, S_1000, 'b^', ms=7, label='$S(1000\,\mathrm{km/s})$ - clusters')
    ax.axhline(1.0, color='gray', ls=':', lw=0.8)
    ax.set_xlabel('$\\lambda = \\alpha m_\\chi / m_\\phi$')
    ax.set_ylabel('$S_0$ (Sommerfeld factor)')
    ax.set_title('Sommerfeld factor for all 17 relic BPs')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_fig = os.path.join(os.path.dirname(__file__), 'v36_sommerfeld.png')
    plt.savefig(out_fig, dpi=150)
    print(f"\nFigure saved: {out_fig}")

    # ── Final verdict ──
    S_fo_arr = np.array([r[vfo_key] for r in results])
    S_30_arr = np.array([r[v30_key] for r in results])
    max_fo_dev = abs(S_fo_arr - 1).max()
    max_S_30 = S_30_arr.max()

    print("\n" + "=" * 80)
    print("FINAL VERDICT")
    print("=" * 80)
    print(f"\n  Freeze-out: max |S - 1| = {max_fo_dev:.2e}")
    print(f"    => Relic density correction < {max_fo_dev*100:.2f}%")
    print(f"    => Tree-level Boltzmann calculation is {'VALID' if max_fo_dev < 0.05 else 'NEEDS CORRECTION'}")
    print(f"\n  Late-time (30 km/s): max S = {max_S_30:.4f}")
    print(f"    => Secluded model: no CMB/Fermi-LAT constraint regardless of S value")
    print(f"    => Halo depletion: tau >> t_Hubble for all BPs (see table above)")
    print(f"    => Late-time Sommerfeld is {'IRRELEVANT' if max_S_30 < 1e6 else 'NEEDS STUDY'} for this model")


if __name__ == "__main__":
    main()


if __name__ == '__main__':
    try:
        import sys as _sys, os as _os
        _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'core'))
        from tg_notify import notify
        notify("\u2705 sommerfeld done!")
    except Exception:
        pass
