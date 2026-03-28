#!/usr/bin/env python3
"""Quick run of condition3 with NEW BP1 parameters (m_χ=54.556, m_φ=12.975 MeV)."""
import math, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'core'))

import condition3_cannibal_sensitivity as c3

# Override with new BP1
c3.M_CHI = 54.556
c3.M_PHI = 12.975e-3  # 12.975 MeV
c3.ALPHA_S = 2.645e-3
c3.ALPHA_P = 2.645e-3
c3.SV0 = 2.0 * math.pi * c3.ALPHA_S * c3.ALPHA_P / c3.M_CHI**2

print("=" * 70)
print("CONDITION 3 — NEW BP1")
print(f"  m_chi = {c3.M_CHI} GeV,  m_phi = {c3.M_PHI*1e3} MeV")
print(f"  alpha = {c3.ALPHA_S:.3e},  SV0 = {c3.SV0:.3e} GeV^-2")
print("=" * 70)

ratios = [0.5, 0.65, 1.0, 1.3, 1.5, 1.7, 2.0, 3.0, 5.0]
print(f"{'mu3/mphi':>10}  {'mu3 [GeV]':>12}  {'Omega_phi':>12}  {'xi_bbn':>8}  {'T_cann_fo':>12}  Status")
print("-" * 75)

for r in ratios:
    mu3 = r * c3.M_PHI
    ana = c3.cannibal_analysis(mu3)
    Y = ana['Y_phi_bbn']
    xi = ana['xi_bbn']
    Tcf = ana['T_cann_fo']
    omega = c3.Y_to_omega(Y, c3.M_PHI)

    if omega > 0.12:
        st = "OVERCLOSURE"
    elif omega > 0.001:
        st = "subdominant"
    else:
        st = "depleted"
    Ts = f"{Tcf*1e3:.2f} MeV" if Tcf > 0 else "---"
    print(f"{r:10.2f}  {mu3:12.4e}  {omega:12.4e}  {xi:8.4f}  {Ts:>12}  {st}")
