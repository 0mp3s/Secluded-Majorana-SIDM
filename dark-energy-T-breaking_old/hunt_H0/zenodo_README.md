# Layer 8: Deriving H₀ from a Dark Matter Lagrangian via Coupled σ–Friedmann ODE

## Description

Numerical solver for the coupled Klein-Gordon + Friedmann system that derives the Hubble constant H₀ as an **output** of a self-interacting dark matter (SIDM) Lagrangian extended with a dark QCD axion sector.

The dark axion field σ — a pseudo-Nambu-Goldstone boson of dark SU(2) chiral symmetry breaking — evolves in the potential V(σ) = Λ_d⁴(1 − cos(σ/f)) from reheating to today. The Hubble parameter H₀ = H(t₀) emerges from the numerical solution without being assumed.

## Key Result

With 6 Lagrangian parameters (m_χ, m_φ, α_D, f, Λ_d, θ_i) and observed Ω_χh² = 0.120:

| θ_i (rad) | H₀ (km/s/Mpc) | Ω_DE | w_σ |
|---|---|---|---|
| 3.0 | 71.1 | 0.72 | −0.86 |
| 3.1 | 73.1 | 0.73 | −0.99 |

The field naturally produces **dynamical dark energy** (w ≠ −1), distinguishable from ΛCDM by DESI/Euclid/CMB-S4.

## Physics

- The same Lagrangian that solves SIDM phenomenology (velocity-dependent cross sections, relic density Ωh² = 0.120, A₄ CP angle θ = 19.47°) also contains a dark axion that drives late-time acceleration.
- Λ_d ~ 2 meV = √(H₀ M_Pl) — the neutrino mass scale coincidence.
- Hilltop quintessence: θ_i ≈ π required for observed H₀.

## Files

- `layer8_cosmic_ode.py` — Main solver: coupled ODE (σ, dσ/dN) via scipy.integrate.solve_ivp
- `README.md` — 7-path strategy document for deriving H₀
- `layer8_evolution.png` — Four-panel diagnostic plot (θ(N), H(N), ρ(N), w(N))

## Requirements

- Python ≥ 3.9
- numpy, scipy, matplotlib
- Depends on `core/v27_boltzmann_relic.py` from the parent repository for Layer 7 (relic density)

## Usage

```python
from layer8_cosmic_ode import solve_layer8

result = solve_layer8(
    m_chi=98.19,        # GeV
    m_phi_GeV=9.66e-3,  # GeV
    alpha_d=3.274e-3,
    f=0.27 * 2.435e18,  # 0.27 M_Pl in GeV
    Lambda_d=2.0e-12,    # GeV (2 meV)
    theta_i=3.0,         # rad (near π)
    omega_chi_h2=0.120,
)

print(f"H₀ = {result.H0_kms:.2f} km/s/Mpc")
print(f"w_σ = {result.w_sigma:.4f}")
```

## Citation

If you use this code, please cite:

> Omer P., "Layer 8: Deriving H₀ from coupled dark axion–Friedmann dynamics in secluded Majorana SIDM", Zenodo (2026). DOI: [to be assigned]

## License

CC-BY-4.0
