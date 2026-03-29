"""
Test 42 — Degeneracy Ridge Plot: (ε, θ_i) Parameter Space
===========================================================

Uses Test 39 CSV data to create a 2D contour plot showing:
  - Σχ² surface in (ε, θ_i) space
  - The degeneracy ridge θ_i* ≈ 2.95 + 0.7|ε|
  - 1σ, 2σ, 3σ contours (Δχ² = 2.30, 6.18, 11.83 for 2 dof)
  - Best-fit point marked

Produces both PNG and PDF (paper-ready).
"""

import numpy as np
import os, sys

sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata

# ── Load Test 39 results ────────────────────────────────────────────────
csv_path = os.path.join(os.path.dirname(__file__), 'test39_fine_results.csv')

data = np.genfromtxt(csv_path, delimiter=',', names=True)

eps_data   = data['eps']
theta_data = data['theta_i']
chi2_data  = data['chi2_tot']
w0_data    = data['w0']
wa_data    = data['wa']

# Best-fit point
best_idx = np.argmin(chi2_data)
chi2_min = chi2_data[best_idx]
eps_best = eps_data[best_idx]
theta_best = theta_data[best_idx]
w0_best = w0_data[best_idx]
wa_best = wa_data[best_idx]

print(f"Best fit: ε={eps_best:.3f}, θ_i={theta_best:.4f}, Σχ²={chi2_min:.2f}")
print(f"          w₀={w0_best:.4f}, w_a={wa_best:.4f}")

# ── Create interpolated 2D grid ────────────────────────────────────────
eps_grid   = np.linspace(eps_data.min(), eps_data.max(), 100)
theta_grid = np.linspace(theta_data.min(), theta_data.max(), 100)
EPS, THETA = np.meshgrid(eps_grid, theta_grid)

# Δχ² relative to minimum
delta_chi2 = chi2_data - chi2_min

CHI2_GRID = griddata(
    (eps_data, theta_data), delta_chi2,
    (EPS, THETA), method='cubic', fill_value=np.nan
)

# Also interpolate w0, wa for contour labels
W0_GRID = griddata(
    (eps_data, theta_data), w0_data,
    (EPS, THETA), method='cubic', fill_value=np.nan
)

# ── Confidence levels for 2 dof ────────────────────────────────────────
# For 2 free parameters (ε, θ_i), Δχ² thresholds:
CL_1sigma = 2.30     # 68.3% CL
CL_2sigma = 6.18     # 95.4% CL  
CL_3sigma = 11.83    # 99.7% CL

# ── Figure 1: Σχ² contour in (ε, θ_i) space ────────────────────────────
fig, ax = plt.subplots(figsize=(9, 7))

# Filled contours of Δχ²
levels = np.array([0, 0.5, CL_1sigma, 4.0, CL_2sigma, 8.0, CL_3sigma, 15, 20, 30])
cf = ax.contourf(-EPS, THETA, CHI2_GRID, levels=levels,
                 cmap='RdYlGn_r', alpha=0.8, extend='max')
cbar = plt.colorbar(cf, ax=ax, label=r'$\Delta\chi^2 = \Sigma\chi^2 - \Sigma\chi^2_{\min}$')

# Contour lines for CL
cs = ax.contour(-EPS, THETA, CHI2_GRID,
                levels=[CL_1sigma, CL_2sigma, CL_3sigma],
                colors=['green', 'orange', 'red'],
                linewidths=[2.0, 1.5, 1.0],
                linestyles=['-', '--', ':'])
ax.clabel(cs, fmt={CL_1sigma: '1σ', CL_2sigma: '2σ', CL_3sigma: '3σ'},
          fontsize=11, inline=True)

# Best-fit point
ax.plot(-eps_best, theta_best, '*', color='black', markersize=15,
        markeredgecolor='white', markeredgewidth=1.0, zorder=10,
        label=f'Best fit: ε=−{abs(eps_best):.2f}, θ_i={theta_best:.3f}')

# Degeneracy ridge: θ_i* ≈ 2.95 + 0.7|ε|
eps_ridge = np.linspace(0.03, 0.13, 50)
theta_ridge = 2.95 + 0.7 * eps_ridge
ax.plot(eps_ridge, theta_ridge, 'k--', linewidth=1.5, alpha=0.7,
        label=r'Ridge: $\theta_i^* \approx 2.95 + 0.7|\varepsilon|$')

# Labels
ax.set_xlabel(r'$|\varepsilon|$', fontsize=14)
ax.set_ylabel(r'$\theta_i$ [rad]', fontsize=14)
ax.set_title(r'Parameter degeneracy: $\Sigma\chi^2(\varepsilon, \theta_i)$ vs DESI DR2'
             f'\nBest fit: Σχ²={chi2_min:.2f}, w₀={w0_best:.4f}, w_a={wa_best:.4f}',
             fontsize=12)
ax.legend(loc='upper left', fontsize=10, framealpha=0.9)
ax.tick_params(labelsize=11)

# Inset: w₀ contours
w0_levels = [-0.76, -0.74, -0.72, -0.70]
cs_w0 = ax.contour(-EPS, THETA, W0_GRID, levels=w0_levels,
                   colors='blue', linewidths=0.8, linestyles='-', alpha=0.5)
ax.clabel(cs_w0, fmt={v: f'w₀={v:.2f}' for v in w0_levels},
          fontsize=8, inline=True, colors='blue')

plt.tight_layout()

# Save
fig_path = os.path.join(os.path.dirname(__file__), 'test42_degeneracy_ridge.png')
fig.savefig(fig_path, dpi=200, bbox_inches='tight')
print(f"PNG saved: {fig_path}")

pdf_path = fig_path.replace('.png', '.pdf')
fig.savefig(pdf_path, bbox_inches='tight')
print(f"PDF saved: {pdf_path}")

plt.close()

# ── Figure 2: 1D slice along ridge ─────────────────────────────────────
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Left: Σχ² vs |ε| (at closest θ_i to ridge)
eps_unique = np.unique(eps_data)
chi2_on_ridge = []
theta_on_ridge = []
w0_on_ridge = []
wa_on_ridge = []

for e in eps_unique:
    mask = eps_data == e
    idx_best_e = np.argmin(chi2_data[mask])
    chi2_on_ridge.append(chi2_data[mask][idx_best_e])
    theta_on_ridge.append(theta_data[mask][idx_best_e])
    w0_on_ridge.append(w0_data[mask][idx_best_e])
    wa_on_ridge.append(wa_data[mask][idx_best_e])

ax1.plot(-eps_unique, chi2_on_ridge, 'ro-', markersize=7, linewidth=2)
ax1.axhline(chi2_min + CL_1sigma, color='green', linestyle='--', alpha=0.7, label='1σ')
ax1.axhline(chi2_min + CL_2sigma, color='orange', linestyle='--', alpha=0.7, label='2σ')
ax1.set_xlabel(r'$|\varepsilon|$', fontsize=13)
ax1.set_ylabel(r'$\Sigma\chi^2_{\min}$ (optimized over $\theta_i$)', fontsize=12)
ax1.set_title(r'$\Sigma\chi^2$ along the ridge', fontsize=12)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.tick_params(labelsize=11)

# Right: w₀ and w_a along ridge
ax2_twin = ax2.twinx()
l1, = ax2.plot(-eps_unique, w0_on_ridge, 'b^-', markersize=7, linewidth=2, label='w₀')
l2, = ax2_twin.plot(-eps_unique, wa_on_ridge, 'rs-', markersize=7, linewidth=2, label='w_a')

# DESI DR2 bands
ax2.axhspan(-0.838-0.055, -0.838+0.055, alpha=0.15, color='blue', label='DESI PP w₀ ±1σ')
ax2_twin.axhspan(-0.62-0.205, -0.62+0.205, alpha=0.15, color='red', label='DESI PP w_a ±1σ')

ax2.set_xlabel(r'$|\varepsilon|$', fontsize=13)
ax2.set_ylabel(r'$w_0$ (CPL)', fontsize=13, color='blue')
ax2_twin.set_ylabel(r'$w_a$ (CPL)', fontsize=13, color='red')
ax2.set_title(r'CPL parameters along the ridge', fontsize=12)
ax2.tick_params(axis='y', colors='blue', labelsize=11)
ax2_twin.tick_params(axis='y', colors='red', labelsize=11)
ax2.tick_params(axis='x', labelsize=11)

lines = [l1, l2]
ax2.legend(lines, [l.get_label() for l in lines], loc='upper left', fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()

fig2_path = os.path.join(os.path.dirname(__file__), 'test42_ridge_slices.png')
fig2.savefig(fig2_path, dpi=200, bbox_inches='tight')
print(f"Ridge slices PNG: {fig2_path}")

fig2_pdf = fig2_path.replace('.png', '.pdf')
fig2.savefig(fig2_pdf, bbox_inches='tight')
print(f"Ridge slices PDF: {fig2_pdf}")

plt.close()

# ── Print ridge data ────────────────────────────────────────────────────
print(f"\n{'='*70}")
print(f"  Ridge Summary (best θ_i for each ε)")
print(f"{'='*70}")
print(f"  {'|ε|':>6}  {'θ_i*':>6}  {'w₀':>8}  {'w_a':>8}  {'Σχ²':>7}")
print(f"  " + "-" * 45)
for i, e in enumerate(eps_unique):
    print(f"  {-e:6.3f}  {theta_on_ridge[i]:6.3f}  {w0_on_ridge[i]:8.4f}  {wa_on_ridge[i]:8.4f}  {chi2_on_ridge[i]:7.2f}")

# Linear fit to ridge
from numpy.polynomial import polynomial as P
coeffs = P.polyfit(-eps_unique, theta_on_ridge, 1)
print(f"\n  Linear fit: θ_i* = {coeffs[0]:.3f} + {coeffs[1]:.3f}·|ε|")
print(f"  (Claimed: θ_i* ≈ 2.95 + 0.7|ε|)")

print(f"\n  Δχ² along ridge: {max(chi2_on_ridge) - min(chi2_on_ridge):.2f}")
print(f"  → The ridge is FLAT: Σχ² varies by only ~{max(chi2_on_ridge) - min(chi2_on_ridge):.1f}")
print(f"     over |ε| = {-eps_unique[-1]:.2f}–{-eps_unique[0]:.2f}")
print(f"  → ε and θ_i are DEGENERATE along this line")
