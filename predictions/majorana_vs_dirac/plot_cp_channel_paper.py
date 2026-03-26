#!/usr/bin/env python3
"""
Generate paper-quality CP-channel figure for the preprint.

2-panel figure (MAP_relic only):
  Top: R_triplet(v) with observation markers and tension band
  Bottom: σ_Maj/σ_Dir ratio with observation markers

Saves to: arxiv/figures/cp_channel_paper.pdf
"""
import sys, os
import numpy as np

_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_DIR, '..', '..')
sys.path.insert(0, os.path.join(_ROOT, 'core'))
sys.path.insert(0, os.path.join(_ROOT, 'cross_checks'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

from global_config import GC
from cp_channel_tension import sigma_channels

# ── Parameters ──
bp = GC.benchmark("MAP_relic")
mc = bp['m_chi_GeV']
mp = bp['m_phi_MeV'] * 1e-3  # GeV
al = bp['alpha']

OBS = GC.observations_as_tuples()
obs_names = [o[0] for o in OBS]
obs_v = np.array([o[1] for o in OBS])

# ── Compute dense grid ──
N = 500
V = np.logspace(np.log10(3.0), np.log10(5000.0), N)

sig_even = np.empty(N)
sig_odd = np.empty(N)
for i, v in enumerate(V):
    se, so = sigma_channels(mc, mp, al, float(v))
    sig_even[i] = se
    sig_odd[i] = so

sig_maj = sig_even + 3 * sig_odd
sig_dir = sig_even + sig_odd

R_triplet = np.where(sig_maj > 0, 3 * sig_odd / sig_maj, 0) * 100  # percent
R_ratio = np.where(sig_dir > 0, sig_maj / sig_dir, 1.0)

# ── Observation markers ──
obs_R_trip = []
obs_ratio = []
for v in obs_v:
    se, so = sigma_channels(mc, mp, al, float(v))
    sm = se + 3 * so
    sd = se + so
    obs_R_trip.append(3 * so / sm * 100 if sm > 0 else 0)
    obs_ratio.append(sm / sd if sd > 0 else 1)
obs_R_trip = np.array(obs_R_trip)
obs_ratio = np.array(obs_ratio)

# ── Short labels for observations ──
short_names = {
    'Draco dSph': 'Draco',
    'Fornax dSph': 'Fornax',
    'NGC 2976': 'NGC 2976',
    'NGC 1560': 'NGC 1560',
    'IC 2574': 'IC 2574',
    'NGC 720 (group)': 'NGC 720',
    'NGC 1332 (group)': 'NGC 1332',
    'Abell 611': 'A611',
    'Abell 2537': 'A2537',
    'Diverse RC band': 'Div. RC',
    'Bullet Cluster': 'Bullet',
    '72 cluster mergers': '72 mergers',
    'TBTF dwarfs': 'TBTF',
}

# ── Figure ──
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.375, 4.5), sharex=True,
                                gridspec_kw={'hspace': 0.08})

# Colors
c_main = '#c0392b'
c_band = '#f39c12'
c_obs = '#2c3e50'
c_stat = '#7f8c8d'

# ── Top: R_triplet ──
ax1.semilogx(V, R_triplet, color=c_main, lw=1.3)
ax1.axhspan(0, 100, xmin=0, xmax=1, alpha=0)  # dummy for xlim
ax1.axvspan(40, 80, alpha=0.12, color=c_band, zorder=0)
ax1.axhline(75, color=c_stat, ls=':', lw=0.7, alpha=0.6)
ax1.text(3500, 76.5, 'statistical\nequipartition', fontsize=5.5,
         color=c_stat, ha='right', va='bottom')

# Observation markers
for i, (v, rt, nm) in enumerate(zip(obs_v, obs_R_trip, obs_names)):
    label = short_names.get(nm, nm)
    ax1.plot(v, rt, 'o', ms=3, color=c_obs, zorder=5)
    # Selective labeling to avoid clutter
    if nm in ('NGC 1560', 'Fornax dSph', 'Abell 611', 'Bullet Cluster', 'TBTF dwarfs'):
        offset_y = 4 if rt < 70 else -6
        offset_x = 1.15 if v > 100 else 0.85
        ax1.annotate(label, (v, rt), fontsize=4.5,
                     xytext=(v * offset_x, rt + offset_y),
                     ha='center', color=c_obs, alpha=0.8)

ax1.set_ylabel(r'$R_{\rm triplet}$ [%]', fontsize=9)
ax1.set_ylim(-3, 100)
ax1.set_xlim(3, 5000)
ax1.text(0.03, 0.92, r'MAP$_{\rm relic}$: $\lambda=38.3$',
         transform=ax1.transAxes, fontsize=7, va='top',
         bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.8))
ax1.text(52, 85, 'tension\nregion', fontsize=5.5, ha='center',
         color='#e67e22', alpha=0.8)

# ── Bottom: Maj/Dir ratio ──
ax2.semilogx(V, R_ratio, color='#2c3e50', lw=1.3)
ax2.axhline(1.0, color='gray', ls='--', lw=0.6)
ax2.axhline(2.0, color=c_stat, ls=':', lw=0.6, alpha=0.5)
ax2.axvspan(40, 80, alpha=0.12, color=c_band, zorder=0)

for i, (v, rr, nm) in enumerate(zip(obs_v, obs_ratio, obs_names)):
    label = short_names.get(nm, nm)
    ax2.plot(v, rr, 'o', ms=3, color=c_obs, zorder=5)
    if nm in ('NGC 1560', 'Fornax dSph', 'Abell 611', 'Bullet Cluster'):
        offset_y = 0.06 if rr < 2 else -0.08
        ax2.annotate(label, (v, rr), fontsize=4.5,
                     xytext=(v * 1.15, rr + offset_y),
                     ha='center', color=c_obs, alpha=0.8)

ax2.set_ylabel(r'$\sigma_T^{\rm Maj}/\sigma_T^{\rm Dir}$', fontsize=9)
ax2.set_xlabel(r'$v$ [km/s]', fontsize=9)
ax2.set_ylim(0.9, 2.35)
ax2.set_xlim(3, 5000)

# Annotation: peak ratio
i_peak = np.argmax(R_ratio)
ax2.annotate(f'$R = {R_ratio[i_peak]:.1f}$',
             xy=(V[i_peak], R_ratio[i_peak]),
             xytext=(V[i_peak] * 3, R_ratio[i_peak] - 0.15),
             fontsize=6, color=c_main,
             arrowprops=dict(arrowstyle='->', color=c_main, lw=0.7))

for ax in (ax1, ax2):
    ax.tick_params(labelsize=7, which='both', direction='in', top=True, right=True)

plt.tight_layout()

out_dir = os.path.join(_ROOT, 'arxiv', 'figures')
os.makedirs(out_dir, exist_ok=True)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(out_dir, f'cp_channel_paper.{ext}'),
                dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Saved to arxiv/figures/cp_channel_paper.{{pdf,png}}")
