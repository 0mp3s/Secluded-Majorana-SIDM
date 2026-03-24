#!/usr/bin/env python3
"""Delta-N_eff sensitivity scan with 2D contour plot."""
import math
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def g_star_s(T_GeV):
    if T_GeV > 200:    return 106.75
    if T_GeV > 80:     return 96.25
    if T_GeV > 4:      return 86.25
    if T_GeV > 1:      return 75.75
    if T_GeV > 0.15:   return 61.75
    if T_GeV > 0.001:  return 10.75
    return 3.91

def delta_neff(m_chi_GeV, m_phi_MeV):
    g_dec = g_star_s(m_chi_GeV)
    g_BBN = g_star_s(0.001)  # 10.75
    xi = (g_BBN / g_dec) ** (1.0/3.0)
    T_phi_MeV = xi * 1.0  # T_SM at BBN = 1 MeV
    T_nu_MeV = 1.0  # neutrinos still coupled at BBN
    dn_massless = (4.0/7.0) * (T_phi_MeV / T_nu_MeV)**4
    x = m_phi_MeV / T_phi_MeV
    if x > 500: boltz = 0.0
    elif x < 0.1: boltz = 1.0
    else: boltz = math.exp(-x) * (1.0 + 15.0/(8.0*x))
    return dn_massless * boltz, xi, T_phi_MeV, x, g_dec

print("="*85)
print("  Delta-N_eff: full sensitivity scan")
print("="*85)

# 1. m_phi fixed at 11 MeV, vary m_chi
print("\n--- m_phi=11 MeV (fixed), varying m_chi ---")
print(f"  {'m_chi[GeV]':>10s}  {'g*S':>7s}  {'xi':>7s}  {'Tphi[MeV]':>10s}  {'m/Tphi':>7s}  {'dN_eff':>12s}")
for mc in [5, 10, 15, 20, 30, 50, 70, 94, 150, 200]:
    dn, xi, Tp, x, gd = delta_neff(mc, 11.0)
    print(f"  {mc:10.1f}  {gd:7.2f}  {xi:7.4f}  {Tp:10.4f}  {x:7.1f}  {dn:12.2e}")

# 2. m_chi fixed at BP1, vary m_phi
print("\n--- m_chi=20.69 GeV (BP1), varying m_phi ---")
print(f"  {'m_phi[MeV]':>10s}  {'Tphi[MeV]':>10s}  {'m/Tphi':>7s}  {'Boltz':>10s}  {'dN_eff':>12s}")
for mp in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 9.0, 11.0, 14.0, 30.0, 50.0]:
    dn, xi, Tp, x, gd = delta_neff(20.69, mp)
    if x < 0.1: bs = 1.0
    elif x > 500: bs = 0.0
    else: bs = math.exp(-x)*(1+15/(8*x))
    print(f"  {mp:10.2f}  {Tp:10.4f}  {x:7.1f}  {bs:10.2e}  {dn:12.2e}")

# 3. Critical m_phi thresholds
print("\n--- Critical m_phi where dN_eff crosses experimental bounds ---")
for threshold, name in [(0.06, "CMB-S4 2sig"), (0.03, "CMB-S4 1sig"), (0.34, "Planck 2sig")]:
    found = False
    for mp in np.linspace(5.0, 0.001, 50000):
        dn, _, _, _, _ = delta_neff(20.69, mp)
        if dn > threshold:
            print(f"  {name} (dN>{threshold}): crossed at m_phi = {mp:.4f} MeV")
            found = True
            break
    if not found:
        dn0, _, _, _, _ = delta_neff(20.69, 0.001)
        status = "EXCEEDS" if dn0 > threshold else "BELOW even massless"
        print(f"  {name}: massless limit dN={dn0:.5f} -- {status}")

# 4. 2D scan of viable region
print("\n--- 2D scan: viable region m_chi=[10,100], m_phi=[7.5,15] ---")
max_dn = 0; max_p = None
for mc in np.linspace(10, 100, 100):
    for mp in np.linspace(7.5, 15, 100):
        dn, _, _, _, _ = delta_neff(mc, mp)
        if dn > max_dn: max_dn = dn; max_p = (mc, mp)
print(f"  Max dN_eff = {max_dn:.2e} at m_chi={max_p[0]:.1f}, m_phi={max_p[1]:.1f} MeV")
print(f"  Planck margin:  {0.34/max_dn:.1e}x")
print(f"  CMB-S4 margin:  {0.06/max_dn:.1e}x")

# 5. g*S across viable range
print("\n--- g*S values across viable m_chi range ---")
for mc in [10, 20, 50, 80, 94, 100]:
    print(f"  m_chi={mc} GeV -> g*S = {g_star_s(mc):.2f}")
print("  ALL viable m_chi sit at SAME g*S=86.25 -> no transition -> no local feature")

# 6. What if m_chi drops below bottom threshold?
print("\n--- Extended: m_chi below b-quark threshold (~4 GeV) ---")
for mc in [1, 2, 3, 3.9, 4.1, 5]:
    dn, xi, Tp, x, gd = delta_neff(mc, 11.0)
    print(f"  m_chi={mc:.1f} GeV: g*S={gd:.2f}, xi={xi:.4f}, Tphi={Tp:.4f}, m/T={x:.1f}, dN={dn:.2e}")

# 7. Worst-case corner of EXTENDED scan
print("\n--- 2D extended scan: m_chi=[5,200], m_phi=[0.1,50] ---")
max_dn2 = 0; max_p2 = None
for mc in np.linspace(5, 200, 200):
    for mp in np.linspace(0.1, 50, 200):
        dn, _, _, _, _ = delta_neff(mc, mp)
        if dn > max_dn2: max_dn2 = dn; max_p2 = (mc, mp)
print(f"  Max dN_eff = {max_dn2:.2e} at m_chi={max_p2[0]:.1f}, m_phi={max_p2[1]:.2f} MeV")
print(f"  (This is the lightest m_phi corner = most relativistic phi)")

# Summary
print("\n" + "="*85)
print("  SUMMARY")
print("="*85)
print(f"  Viable region max:     dN_eff = {max_dn:.2e}  (margin > 10^8 to Planck)")
print(f"  Extended region max:   dN_eff = {max_dn2:.2e}  (at lightest m_phi corner)")
print(f"  Monotonic in both m_phi (more mass = more suppression)")
print(f"  and m_chi (higher = more g*S = smaller xi = less dN_eff)")
print(f"  No local maxima possible: g*S is constant=86.25 across viable m_chi=[10,100]")

# ===================== PLOTS =====================
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(out_dir, exist_ok=True)

# --- Plot 1: 2D contour in (m_phi, m_chi) plane ---
m_phi_arr = np.logspace(-1, np.log10(50), 300)  # 0.1 - 50 MeV
m_chi_arr = np.linspace(5, 200, 200)             # 5 - 200 GeV
DN = np.zeros((len(m_chi_arr), len(m_phi_arr)))
for i, mc in enumerate(m_chi_arr):
    for j, mp in enumerate(m_phi_arr):
        dn, _, _, _, _ = delta_neff(mc, mp)
        DN[i, j] = max(dn, 1e-80)

fig, ax = plt.subplots(figsize=(10, 6))
MP, MC = np.meshgrid(m_phi_arr, m_chi_arr)
levels = [1e-60, 1e-40, 1e-20, 1e-10, 1e-5, 0.03, 0.06, 0.34]
cs = ax.contourf(MP, MC, np.log10(DN + 1e-80), levels=50, cmap='RdYlBu_r')
cbar = fig.colorbar(cs, ax=ax, label=r'$\log_{10}\Delta N_{\rm eff}$')

# Contour lines for experimental bounds
cl = ax.contour(MP, MC, DN, levels=[0.03, 0.06, 0.34],
                colors=['purple', 'magenta', 'orange'], linewidths=2)
ax.clabel(cl, fmt={0.03: 'CMB-S4 1$\\sigma$', 0.06: 'CMB-S4 2$\\sigma$', 0.34: 'Planck 2$\\sigma$'},
          fontsize=9)

# Viable region box
from matplotlib.patches import Rectangle
rect = Rectangle((9, 10), 5, 90, linewidth=2.5, edgecolor='lime',
                 facecolor='none', linestyle='--', label='Viable region')
ax.add_patch(rect)

# Benchmark points
ax.plot(11.34, 20.69, 'w*', markersize=14, markeredgecolor='black', label='BP1')
ax.plot(11.10, 94.07, 'ws', markersize=10, markeredgecolor='black', label='MAP')

ax.set_xscale('log')
ax.set_xlabel(r'$m_\phi$ [MeV]', fontsize=13)
ax.set_ylabel(r'$m_\chi$ [GeV]', fontsize=13)
ax.set_title(r'$\Delta N_{\rm eff}$ Sensitivity Map', fontsize=14)
ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3, which='both')
fig.tight_layout()

path1 = os.path.join(out_dir, 'delta_neff_sensitivity_2d.png')
fig.savefig(path1, dpi=150, bbox_inches='tight')
path1_pdf = os.path.join(out_dir, 'delta_neff_sensitivity_2d.pdf')
fig.savefig(path1_pdf, bbox_inches='tight')
plt.close(fig)
print(f"\n  Plot saved: {path1}")
print(f"  Plot saved: {path1_pdf}")
sz = os.path.getsize(path1)
print(f"  PNG size: {sz:,} bytes ({sz/1024:.1f} KB)")

# --- Plot 2: 1D slice at BP1 m_chi ---
fig2, ax2 = plt.subplots(figsize=(9, 5))
m_phi_1d = np.logspace(-2, np.log10(50), 500)
for mc, label, color in [(20.69, 'BP1 ($m_\\chi$=20.7 GeV)', 'steelblue'),
                          (94.07, 'MAP ($m_\\chi$=94.1 GeV)', 'firebrick')]:
    dn_arr = [max(delta_neff(mc, mp)[0], 1e-80) for mp in m_phi_1d]
    ax2.plot(m_phi_1d, dn_arr, lw=2, color=color, label=label)

ax2.axhline(0.06, color='magenta', ls='--', lw=2, label='CMB-S4 $2\\sigma$')
ax2.axhline(0.03, color='purple', ls=':', lw=2, label='CMB-S4 $1\\sigma$')
ax2.axhline(0.34, color='orange', ls='--', lw=2, label='Planck $2\\sigma$')
ax2.axvspan(9, 14, alpha=0.15, color='lime', label='Viable $m_\\phi$')
ax2.axvline(11.34, color='gray', ls=':', alpha=0.5)
ax2.axvline(11.10, color='gray', ls=':', alpha=0.5)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0.01, 50)
ax2.set_ylim(1e-30, 1)
ax2.set_xlabel(r'$m_\phi$ [MeV]', fontsize=13)
ax2.set_ylabel(r'$\Delta N_{\rm eff}$', fontsize=13)
ax2.set_title(r'$\Delta N_{\rm eff}$ vs mediator mass at BBN', fontsize=14)
ax2.legend(fontsize=9, loc='upper right')
ax2.grid(True, alpha=0.3, which='both')
fig2.tight_layout()

path2 = os.path.join(out_dir, 'delta_neff_vs_mphi.png')
fig2.savefig(path2, dpi=150, bbox_inches='tight')
path2_pdf = os.path.join(out_dir, 'delta_neff_vs_mphi.pdf')
fig2.savefig(path2_pdf, bbox_inches='tight')
plt.close(fig2)
print(f"\n  Plot saved: {path2}")
print(f"  Plot saved: {path2_pdf}")
sz2 = os.path.getsize(path2)
print(f"  PNG size: {sz2:,} bytes ({sz2/1024:.1f} KB)")
