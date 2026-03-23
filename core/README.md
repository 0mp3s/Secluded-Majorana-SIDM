# core/ — Shared Numerical Solvers

## v22_raw_scan.py — VPM Phase-Shift Solver

The Variable Phase Method (VPM) solver for self-interacting dark matter
transfer cross sections. Solves the radial Schrödinger equation for Yukawa
potential with partial-wave summation up to convergence.

**Key function:**
```python
sigma_T_vpm(m_chi_GeV, m_phi_GeV, alpha, v_km_s) → σ/m [cm²/g]
```

**Parameters:**
| Parameter   | Unit | Description |
|-------------|------|-------------|
| `m_chi_GeV` | GeV  | Dark matter mass |
| `m_phi_GeV` | GeV  | Mediator mass (⚠ GeV, not MeV) |
| `alpha`     | —    | Yukawa coupling α_χ |
| `v_km_s`    | km/s | Relative velocity |

**Also exports:** `vpm_phase_shift`, `C_KM_S`, `GEV2_TO_CM2`, `GEV_IN_G`, `M_CHI_VALS`, `M_PHI_VALS`, `LAM_CRITS`

When run as `__main__`, performs a full parameter-space scan and saves:
- `all_viable_raw_v8.csv` (80,142 SIDM-viable points)
- `all_viable_representative_v8.csv` (1 per cell)

---

## v27_boltzmann_relic.py — Boltzmann Relic Density Solver

Numerical Boltzmann equation solver for s-wave Majorana DM annihilation 
through a scalar mediator. Uses RK4 with g_*(T) tables.

**Key functions:**
```python
solve_boltzmann(m_chi, sigma_v0, ...)  → Y_final
Y_to_omega_h2(Y_inf, m_chi)           → Ω h²
kolb_turner_swave(m_chi, sigma_v0)     → Ω h² (analytic approximation)
```
