#!/usr/bin/env python3
"""
model_validations/smbh_seeds/predict_smbh_seeds.py
===================================================
§7.4: SMBH Seed Formation at High Redshift via Gravothermal Collapse

Physics:
  At high-z, SIDM halos can undergo gravothermal collapse, forming
  a dense central mass concentration — a seed for SMBH growth.

  Key ingredients:
  1. Halo concentration c(M, z) from Correa+2015 (valid z ≤ 20)
  2. NFW ρ_s from c, M200
  3. σ/m(v_vir) from VPM
  4. Relaxation time: t_relax = 1/(ρ_s × σ/m × σ_v)
  5. Collapse time: t_gc ≈ 150 × t_relax (Balberg+2002)
  6. Seed mass: M_seed ≈ f_c × M200 with f_c ≈ 0.01 (Pollack+2015)
  7. Comparison: t_gc < t_universe(z) → collapse possible?

References:
  - Correa+2015 (MNRAS 452, 1217): c(M, z) relation
  - Balberg+2002: gravothermal collapse timescale
  - Pollack+2015: SIDM seed formation
  - Harikane+2023, Maiolino+2024: JWST high-z SMBH observations
"""
import sys, os, math
import numpy as np

# ---------- path bootstrap ----------
_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
# ------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config_loader import load_config
from v22_raw_scan import sigma_T_vpm

# Warm up JIT
sigma_T_vpm(20.0, 10e-3, 1e-3, 100.0)

_DIR = os.path.dirname(os.path.abspath(__file__))

# Constants
MSUN_G = 1.989e33
KPC_CM = 3.086e21
GYR_S = 3.156e16
KM_S_CM_S = 1e5
G_NEWTON = 6.674e-8  # cm³ g⁻¹ s⁻²


def correa_concentration(M200_Msun, z, h=0.674):
    """
    Correa+2015 (MNRAS 452, 1217) concentration-mass-redshift relation.
    Valid for 10^8 < M200/M_sun < 10^15, 0 < z < 20.
    
    c(M, z) = A(z) × (M / 10^12 h^{-1} M_sun)^{B(z)}
    A(z) = 5.71 × (1+z)^{-0.47}
    B(z) = -0.084 - 0.025 × ln(1+z)
    """
    A = 5.71 * (1 + z)**(-0.47)
    B = -0.084 - 0.025 * math.log(1 + z)
    M_pivot = 1e12 / h  # h^{-1} M_sun → M_sun
    c = A * (M200_Msun / M_pivot)**B
    return max(c, 1.5)  # floor at 1.5


def nfw_params(M200_Msun, c200, h=0.674):
    """Compute NFW rho_s [M_sun/kpc³], r_s [kpc], R200 [kpc] from M200 and c200."""
    rho_crit = 277.5 * h**2  # M_sun/kpc³
    R200 = (3 * M200_Msun / (4 * math.pi * 200 * rho_crit))**(1.0/3.0)
    r_s = R200 / c200
    gc = math.log(1 + c200) - c200 / (1 + c200)
    rho_s = M200_Msun / (4 * math.pi * r_s**3 * gc)
    return rho_s, r_s, R200


def virial_velocity(M200_Msun, R200_kpc):
    """Virial velocity [km/s] from M200 and R200."""
    # v_vir² = G M200 / R200
    G_kpc = 4.302e-3  # (km/s)² kpc / M_sun
    return math.sqrt(G_kpc * M200_Msun / R200_kpc)


def age_of_universe(z):
    """Approximate age of universe at redshift z [Gyr] (matter-dominated + Λ)."""
    # Using Planck 2018: H0 = 67.4, Ωm = 0.315, ΩΛ = 0.685
    H0_Gyr = 67.4 * 1.022e-3  # km/s/Mpc → Gyr⁻¹ (1 km/s/Mpc ≈ 1.022e-3 Gyr⁻¹)
    # Numerical integration of dt/dz = -1/((1+z) H(z))
    from scipy.integrate import quad
    Om, OL = 0.315, 0.685
    def integrand(zp):
        return 1.0 / ((1 + zp) * H0_Gyr * math.sqrt(Om * (1 + zp)**3 + OL))
    t, _ = quad(integrand, z, 1000)
    return t


def relaxation_time(rho_s_Msun_kpc3, sigma_over_m_cm2_g, sigma_v_km_s):
    """
    SIDM relaxation time [Gyr].
    t_relax = 1 / (ρ_s × σ/m × σ_v)
    with proper unit conversions.
    """
    rho_cgs = rho_s_Msun_kpc3 * MSUN_G / KPC_CM**3  # g/cm³
    v_cgs = sigma_v_km_s * KM_S_CM_S  # cm/s
    # Γ = σ/m × ρ × v [s⁻¹]
    gamma = sigma_over_m_cm2_g * rho_cgs * v_cgs
    if gamma <= 0:
        return float('inf')
    t_relax_s = 1.0 / gamma
    return t_relax_s / GYR_S


def main():
    cfg = load_config(__file__)
    out_dir = os.path.join(_DIR, cfg.get('output_dir', 'output'))
    os.makedirs(out_dir, exist_ok=True)

    bps = cfg['benchmark_points']
    redshifts = cfg['redshifts']
    halo_masses_log = cfg['halo_masses_log10']
    h = cfg.get('h', 0.674)
    collapse_factor = cfg.get('collapse_factor', 150)
    f_seed = cfg.get('seed_fraction', 0.01)

    print("=" * 100)
    print("  §7.4: SMBH Seed Formation at High Redshift via Gravothermal Collapse")
    print("=" * 100)
    print()
    print(f"  Cosmology: h = {h}, Ωm = 0.315, ΩΛ = 0.685")
    print(f"  Collapse: t_gc = {collapse_factor} × t_relax (Balberg+2002)")
    print(f"  Seed fraction: f_c = {f_seed} (Pollack+2015)")
    print()

    # Concentration table
    print("  ─── Correa+2015 c(M, z) ───")
    header = f"  {'log10(M/M☉)':>14}"
    for z in redshifts:
        header += f"  {'z='+str(z):>8}"
    print(header)
    print("  " + "─" * (14 + 10 * len(redshifts)))
    for lm in halo_masses_log:
        M = 10**lm
        row = f"  {lm:>14.1f}"
        for z in redshifts:
            c = correa_concentration(M, z, h)
            row += f"  {c:>8.2f}"
        print(row)
    print()

    # Universe age
    print("  ─── Age of Universe ───")
    for z in redshifts:
        t = age_of_universe(z)
        print(f"    z = {z:>3}: t_univ = {t:.3f} Gyr = {t*1000:.0f} Myr")
    print()

    # Main computation: for each BP × z × M
    for bp in bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        lam = alpha * m_chi / m_phi

        print(f"  ═══ {label}: m_χ={m_chi} GeV, m_φ={bp['m_phi_MeV']} MeV, "
              f"α={alpha:.3e}, λ={lam:.2f} ═══")
        print()

        print(f"  {'z':>4}  {'log(M)':>7}  {'c':>6}  {'v_vir':>7}  {'σ/m':>8}  "
              f"{'t_relax':>9}  {'t_gc':>9}  {'t_univ':>8}  {'t_gc<t_u?':>10}  "
              f"{'M_seed':>10}  {'Status':>12}")
        print("  " + "─" * 105)

        results = []
        for z in redshifts:
            t_univ = age_of_universe(z)
            for lm in halo_masses_log:
                M200 = 10**lm
                c = correa_concentration(M200, z, h)
                rho_s, r_s, R200 = nfw_params(M200, c, h)
                v_vir = virial_velocity(M200, R200)

                sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_vir)
                t_relax = relaxation_time(rho_s, sigma_m, v_vir / math.sqrt(3))
                t_gc = collapse_factor * t_relax

                can_collapse = t_gc < t_univ
                M_seed = f_seed * M200 if can_collapse else 0

                status = "COLLAPSE" if can_collapse else "NO COLLAPSE"
                if can_collapse and M_seed >= 1e5:
                    status = "SEED ✓"
                elif can_collapse:
                    status = "SMALL SEED"

                results.append({
                    'z': z, 'log_M': lm, 'c': c, 'v_vir': v_vir,
                    'sigma_m': sigma_m, 't_relax': t_relax, 't_gc': t_gc,
                    't_univ': t_univ, 'collapse': can_collapse,
                    'M_seed': M_seed, 'status': status
                })

                M_seed_str = f"{M_seed:.1e}" if M_seed > 0 else "—"
                print(f"  {z:>4}  {lm:>7.1f}  {c:>6.2f}  {v_vir:>6.1f}  "
                      f"{sigma_m:>8.4f}  {t_relax:>8.3f}  {t_gc:>8.2f}  "
                      f"{t_univ:>7.3f}  {'YES' if can_collapse else 'no':>10}  "
                      f"{M_seed_str:>10}  {status:>12}")

            print()  # gap between redshifts

        # Summary for this BP
        n_collapse = sum(1 for r in results if r['collapse'])
        n_seed = sum(1 for r in results if r['M_seed'] >= 1e5)
        print(f"  {label} summary: {n_collapse} collapse cases, "
              f"{n_seed} with M_seed ≥ 10⁵ M☉ (out of {len(results)} total)")
        print()

    # ──── JWST comparison ────
    print("=" * 100)
    print("  JWST Comparison: Can our model produce observed high-z SMBHs?")
    print("=" * 100)
    print()
    print("  JWST observations (Harikane+2023, Maiolino+2024):")
    print("    GN-z11 (z ≈ 10.6): M_BH ≈ 1.6 × 10⁶ M☉")
    print("    UHZ1 (z ≈ 10.1): M_BH ≈ 4 × 10⁷ M☉ (lensed)")
    print("    CEERS-1019 (z ≈ 8.7): M_BH ≈ 9 × 10⁶ M☉")
    print()
    print("  Required seed mass (assuming Eddington growth from z_seed):")
    print("    M_seed(z=15) → M_BH(z=10) with t_grow ≈ 0.17 Gyr:")
    print("    Salpeter time t_S ≈ 45 Myr × (ε/0.1)")
    print("    Growth factor: e^{t_grow/t_S} ≈ e^{3.8} ≈ 45×")
    print("    → M_seed ≈ 10⁴–10⁵ M☉ sufficient for 10⁶–10⁷ BH at z=10")
    print()

    for bp in bps:
        label = bp['label']
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']

        # Check z=15, M=10^{10.5} case
        M_test = 10**10.5
        c = correa_concentration(M_test, 15, h)
        rho_s, r_s, R200 = nfw_params(M_test, c, h)
        v_vir = virial_velocity(M_test, R200)
        sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_vir)
        t_relax = relaxation_time(rho_s, sigma_m, v_vir / math.sqrt(3))
        t_gc = collapse_factor * t_relax
        t_univ = age_of_universe(15)
        M_seed = f_seed * M_test if t_gc < t_univ else 0

        print(f"  {label}: z=15, M=10^10.5 M☉ → c={c:.1f}, v_vir={v_vir:.0f} km/s, "
              f"σ/m={sigma_m:.3f}")
        print(f"    t_gc = {t_gc:.2f} Gyr {'<' if t_gc < t_univ else '>'} "
              f"t_univ = {t_univ:.3f} Gyr → "
              f"{'SEED = ' + f'{M_seed:.1e} M☉' if M_seed > 0 else 'NO SEED'}")
    print()

    # ──── Plot ────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = {'BP1': 'steelblue', 'BP9': 'seagreen', 'MAP': 'firebrick'}
    markers = ['o', 's', '^']

    # Plot 1: t_gc vs z for M = 10^10 M_sun
    ax = axes[0]
    z_arr = np.array(redshifts)
    t_univ_arr = np.array([age_of_universe(z) for z in z_arr])
    ax.plot(z_arr, t_univ_arr, 'k--', lw=2, label='t_universe')

    for k, bp in enumerate(bps):
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        t_gc_arr = []
        for z in z_arr:
            M_test = 1e10
            c = correa_concentration(M_test, z, h)
            rho_s, r_s, R200 = nfw_params(M_test, c, h)
            v_vir = virial_velocity(M_test, R200)
            sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_vir)
            t_relax = relaxation_time(rho_s, sigma_m, v_vir / math.sqrt(3))
            t_gc_arr.append(collapse_factor * t_relax)
        ax.semilogy(z_arr, t_gc_arr, marker=markers[k],
                     color=colors.get(bp['label'], 'gray'), lw=2,
                     label=f"{bp['label']} t_gc")

    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Time [Gyr]')
    ax.set_title('Gravothermal Collapse Time vs Universe Age\n(M = 10¹⁰ M☉)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.invert_xaxis()

    # Plot 2: Seed mass vs halo mass at z=15
    ax = axes[1]
    M_arr = np.logspace(9, 11.5, 50)
    t_univ_15 = age_of_universe(15)

    for k, bp in enumerate(bps):
        m_chi = bp['m_chi_GeV']
        m_phi = bp['m_phi_MeV'] / 1000.0
        alpha = bp['alpha']
        seeds = []
        for M in M_arr:
            c = correa_concentration(M, 15, h)
            rho_s, r_s, R200 = nfw_params(M, c, h)
            v_vir = virial_velocity(M, R200)
            sigma_m = sigma_T_vpm(m_chi, m_phi, alpha, v_vir)
            t_relax = relaxation_time(rho_s, sigma_m, v_vir / math.sqrt(3))
            t_gc = collapse_factor * t_relax
            if t_gc < t_univ_15:
                seeds.append(f_seed * M)
            else:
                seeds.append(0)
        seeds = np.array(seeds)
        valid = seeds > 0
        if np.any(valid):
            ax.loglog(M_arr[valid], seeds[valid], marker=markers[k],
                      color=colors.get(bp['label'], 'gray'), lw=2,
                      label=bp['label'], markersize=4)

    # JWST requirement
    ax.axhline(1e5, color='orange', ls='--', lw=2, label='M_seed = 10⁵ M☉ (JWST need)')
    ax.axhline(1e7, color='red', ls=':', lw=1.5, label='M_seed = 10⁷ M☉')
    ax.set_xlabel('Halo Mass M₂₀₀ [M☉]')
    ax.set_ylabel('Seed Mass [M☉]')
    ax.set_title('SIDM Seed Mass at z = 15')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig_path = os.path.join(out_dir, 'smbh_seeds_prediction.png')
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"  Plot saved: {fig_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
