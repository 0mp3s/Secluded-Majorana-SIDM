# core/ — Shared Numerical Solvers

## config_loader.py — Configuration Utility

Loads `config.json` files for all project scripts.

```python
from config_loader import load_config
cfg = load_config(__file__)                # loads config.json next to the script
cfg = load_config(__file__, "my_cfg.json") # explicit path
# Or via CLI: python script.py --config my_cfg.json
```

Priority: `--config` CLI argument > explicit `config_path` > `config.json` in script's directory.

Returns empty dict `{}` if no config file found (scripts fall back to hardcoded defaults).

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

**Also exports:** `vpm_phase_shift`, `C_KM_S`, `GEV2_TO_CM2`, `GEV_IN_G`

When run as `__main__`, performs a full parameter-space scan (grid configurable via
`--config` or `config.json`):

```jsonc
{
    "grid": {
        "n_chi": 50,       // m_χ grid points
        "n_phi": 70,       // m_φ grid points
        "n_res": 4,        // resonance bands
        "n_alpha": 200,    // α points per band
        "m_chi_range_GeV": [0.1, 100.0],
        "m_phi_range_GeV": [0.1e-3, 200e-3],
        "lambda_crits": [1.68, 6.45, 14.7, 26.0]
    }
}
```

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
