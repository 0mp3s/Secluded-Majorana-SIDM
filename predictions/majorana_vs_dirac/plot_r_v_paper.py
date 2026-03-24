#!/usr/bin/env python3
"""
predictions/majorana_vs_dirac/plot_r_v_paper.py
================================================
Publication-quality 2-panel figure for §7.7:
  Panel (a): R(v) = sigma_T^Maj / sigma_T^Dir for BP1 and MAP
  Panel (b): f_odd(v) = odd-ell fraction

Uses VPM phase-shift solver from core/v22_raw_scan.py.
"""
import sys, os
import numpy as np

_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from v22_raw_scan import sigma_T_vpm
from analyze_ratio_function import sigma_T_dirac, partial_wave_decomposition

# ---- Benchmarks — loaded from global_config.json ----
from global_config import GC
_bp1 = GC.benchmark("BP1")
_map = GC.benchmark("MAP")
BP1 = {"label": r"BP1 ($\lambda=%.1f$)" % (_bp1["alpha"] * _bp1["m_chi_GeV"] / (_bp1["m_phi_MeV"] * 1e-3)),
       "m_chi": _bp1["m_chi_GeV"], "m_phi_MeV": _bp1["m_phi_MeV"], "alpha": _bp1["alpha"], "color": "#1565C0"}
MAP = {"label": r"MAP ($\lambda=%.1f$)" % (_map["alpha"] * _map["m_chi_GeV"] / (_map["m_phi_MeV"] * 1e-3)),
       "m_chi": _map["m_chi_GeV"], "m_phi_MeV": _map["m_phi_MeV"], "alpha": _map["alpha"], "color": "#C62828"}

V_GRID = np.logspace(np.log10(2.0), np.log10(3000.0), 10000)
V_SCALES = {'dSph\n30 km/s': 30.0, 'MW\n220 km/s': 220.0, 'Cluster\n1200 km/s': 1200.0}


def compute_all(bp):
    m_chi = bp['m_chi']
    m_phi = bp['m_phi_MeV'] / 1000.0
    alpha = bp['alpha']
    R     = np.zeros(len(V_GRID))
    f_odd = np.zeros(len(V_GRID))
    for i, v in enumerate(V_GRID):
        sm = sigma_T_vpm(m_chi, m_phi, alpha, v)
        sd = sigma_T_dirac(m_chi, m_phi, alpha, v)
        R[i] = sm / sd if sd > 1e-20 else np.nan
        _, fo, _ = partial_wave_decomposition(m_chi, m_phi, alpha, v)
        f_odd[i] = fo
    return R, f_odd


def main():
    print("Computing BP1...")
    R_bp1, f_bp1 = compute_all(BP1)
    print("Computing MAP...")
    R_map, f_map = compute_all(MAP)
    print("Plotting...")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                                    gridspec_kw={'hspace': 0.08})

    # ---- Panel (a): R(v) ----
    ax1.plot(V_GRID, R_bp1, color=BP1['color'], lw=2, label=BP1['label'])
    ax1.plot(V_GRID, R_map, color=MAP['color'], lw=2, label=MAP['label'])
    ax1.axhline(1.0, color='k', ls='--', lw=0.8, zorder=0)

    # Shade R > 1 region for MAP
    mask = R_map > 1.0
    if np.any(mask):
        ax1.fill_between(V_GRID, 1.0, R_map, where=mask,
                         color=MAP['color'], alpha=0.15, label=r'$R > 1$ (MAP)')

    for label, vv in V_SCALES.items():
        ax1.axvline(vv, color='gray', ls=':', lw=0.6, alpha=0.7)
        ax1.text(vv, 1.19, label, ha='center', va='top', fontsize=7,
                 color='gray', style='italic')

    ax1.set_ylabel(r'$R(v) = \sigma_T^{\rm Maj} / \sigma_T^{\rm Dir}$', fontsize=12)
    ax1.set_ylim(0.4, 1.22)
    ax1.set_xlim(2, 3000)
    ax1.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax1.text(0.02, 0.95, '(a)', transform=ax1.transAxes, fontsize=13, fontweight='bold', va='top')

    # Inset zoom on R > 1 peak
    axins = inset_axes(ax1, width="40%", height="45%", loc='center left',
                       bbox_to_anchor=(0.08, 0.05, 1, 1), bbox_transform=ax1.transAxes)
    axins.plot(V_GRID, R_bp1, color=BP1['color'], lw=1.5)
    axins.plot(V_GRID, R_map, color=MAP['color'], lw=1.5)
    axins.axhline(1.0, color='k', ls='--', lw=0.6)
    if np.any(mask):
        axins.fill_between(V_GRID, 1.0, R_map, where=mask,
                           color=MAP['color'], alpha=0.15)
    axins.set_xlim(25, 200)
    axins.set_ylim(0.82, 1.15)
    axins.set_xscale('log')
    axins.tick_params(labelsize=7)
    axins.set_xlabel(r'$v$ [km/s]', fontsize=7)
    # Mark peak
    i_peak = np.nanargmax(R_map[(V_GRID > 25) & (V_GRID < 200)])
    v_sub = V_GRID[(V_GRID > 25) & (V_GRID < 200)]
    R_sub = R_map[(V_GRID > 25) & (V_GRID < 200)]
    axins.annotate(f'$R = {R_sub[i_peak]:.2f}$\n$v = {v_sub[i_peak]:.0f}$ km/s',
                   xy=(v_sub[i_peak], R_sub[i_peak]),
                   xytext=(120, 1.10), fontsize=7,
                   arrowprops=dict(arrowstyle='->', color=MAP['color'], lw=0.8),
                   color=MAP['color'])

    # ---- Panel (b): f_odd(v) ----
    ax2.plot(V_GRID, f_bp1, color=BP1['color'], lw=2, label=BP1['label'])
    ax2.plot(V_GRID, f_map, color=MAP['color'], lw=2, label=MAP['label'])
    ax2.axhline(0.5, color='k', ls='--', lw=0.8, zorder=0)

    mask_f = f_map > 0.5
    if np.any(mask_f):
        ax2.fill_between(V_GRID, 0.5, f_map, where=mask_f,
                         color=MAP['color'], alpha=0.15, label=r'$f_{\rm odd} > 1/2$ (MAP)')

    for label, vv in V_SCALES.items():
        ax2.axvline(vv, color='gray', ls=':', lw=0.6, alpha=0.7)

    ax2.set_xlabel(r'$v_{\rm rel}$ [km/s]', fontsize=12)
    ax2.set_ylabel(r'$f_{\rm odd}(v)$', fontsize=12)
    ax2.set_xscale('log')
    ax2.set_ylim(0, 0.7)
    ax2.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax2.text(0.02, 0.95, '(b)', transform=ax2.transAxes, fontsize=13, fontweight='bold', va='top')

    # Annotation: R > 1 iff f_odd > 1/2
    ax2.annotate(r'$R > 1 \Leftrightarrow f_{\rm odd} > 1/2$',
                 xy=(57, 0.606), xytext=(300, 0.65), fontsize=9,
                 arrowprops=dict(arrowstyle='->', color='gray', lw=0.8),
                 color='gray', style='italic')

    ax1.set_xscale('log')

    fig.subplots_adjust(hspace=0.08, left=0.12, right=0.95, top=0.97, bottom=0.07)
    os.makedirs(os.path.join(_DIR, 'output'), exist_ok=True)
    for ext in ['png', 'pdf']:
        path = os.path.join(_DIR, 'output', f'r_v_ratio_paper_10k.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        size_bytes = os.path.getsize(path)
        print(f"Saved: {path}  [{size_bytes:,} bytes = {size_bytes/1024:.1f} KB]")
    plt.close()


if __name__ == '__main__':
    main()
