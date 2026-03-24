# core/ — Shared Numerical Solvers

## Mathematical Background

### Yukawa Potential

The DM self-interaction is mediated by a light scalar φ with attractive Yukawa potential:

$$V(r) = -\alpha \frac{e^{-m_\phi r}}{r}$$

where α = y²/(4π) is the Yukawa coupling and m_φ is the mediator mass.

### Variable Phase Method (VPM)

The transfer cross section is computed by solving the radial Schrödinger equation
via the Variable Phase Method. The VPM ODE for partial-wave phase shift δ_l(x) is:

$$\frac{d\delta_l}{dx} = \frac{\lambda\, e^{-x}}{\kappa\, x} \left[\hat{j}_l(\kappa x)\cos\delta_l - \hat{n}_l(\kappa x)\sin\delta_l\right]^2$$

where the dimensionless variables are:

| Symbol | Definition | Description |
|--------|-----------|-------------|
| x | m_φ r | Dimensionless radial coordinate |
| κ | k / m_φ = μv / m_φ | Dimensionless momentum (μ = m_χ/2) |
| λ | α m_χ / m_φ | Dimensionless coupling (**no factor of 2**) |

ĵ_l = z j_l(z) and n̂_l = −z y_l(z) are Riccati–Bessel functions,
with j_l, y_l the standard spherical Bessel functions.

The ODE is integrated from x_min to x_max = 50–100 using 4th-order Runge–Kutta
(Numba JIT-compiled), with initial condition δ_l(x_min) = 0.

### Transfer Cross Section

For **identical Majorana fermions**, only even (l = 0, 2, …) partial waves
contribute to the direct amplitude, and odd (l = 1, 3, …) to the exchange
amplitude, with statistical weight 1 and 3 respectively:

$$\frac{\sigma_T}{m_\chi} = \frac{2\pi}{k^2\, m_\chi} \sum_{l=0}^{l_{\max}} w_l\,(2l+1)\sin^2\delta_l$$

where w_l = 1 for even l, w_l = 3 for odd l. The factor 2π (not 4π) accounts
for identical-particle symmetry. The sum is truncated when the fractional
contribution drops below 10⁻³.

Unit conversion: σ_T/m [cm²/g] = (σ [GeV⁻²] × 3.8938×10⁻²⁸ cm²/GeV⁻²) / (m_χ [GeV] × 1.783×10⁻²⁴ g/GeV).

### Boltzmann Equation (Relic Density)

The freeze-out yield Y = n/s evolves as:

$$\frac{dY}{dx} = -\sqrt{\frac{\pi}{45}}\, g_{\mathrm{eff}}\, M_{\mathrm{Pl}}\, m_\chi\, \frac{\langle\sigma v\rangle}{x^2}\,(Y^2 - Y_{\mathrm{eq}}^2)$$

where x = m_χ/T and g_eff = g_*^S √g_*^ρ (1 + T/(3g_*^S) dg_*^S/dT).

For **s-wave** Majorana annihilation χχ → φφ via t/u-channel:

$$\langle\sigma v\rangle = \frac{\pi\alpha^2}{4\, m_\chi^2}$$

The equilibrium yield (non-relativistic limit):

$$Y_{\mathrm{eq}} = \frac{45}{4\pi^4}\,\frac{g_\chi}{g_{*S}}\, x^{3/2}\, e^{-x}$$

with g_χ = 2 (Majorana). The relic density is:

$$\Omega h^2 = \frac{m_\chi\, s_0\, Y_\infty}{\rho_{\mathrm{crit}}/h^2}$$

where s_0 = 2891.2 cm⁻³ and ρ_crit/h² = 1.054×10⁻⁵ GeV/cm³.

---

## API Reference

### config_loader.py — Configuration Utility

Loads `config.json` files for all project scripts.

```python
from config_loader import load_config
cfg = load_config(__file__)                # loads config.json next to the script
cfg = load_config(__file__, "my_cfg.json") # explicit path
# Or via CLI: python script.py --config my_cfg.json
```

Priority: `--config` CLI argument > explicit `config_path` > `config.json` in script's directory.

Returns empty dict `{}` if no config file found (scripts fall back to hardcoded defaults).

### v22_raw_scan.py — VPM Phase-Shift Solver

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
