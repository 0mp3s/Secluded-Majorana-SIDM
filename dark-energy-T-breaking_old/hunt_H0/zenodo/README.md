# Layer 8: Deriving H₀ from a Dark Matter Lagrangian via Coupled σ–Friedmann ODE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19274637.svg)](https://doi.org/10.5281/zenodo.19274637)

## Description

Numerical solver and analysis suite for deriving the Hubble constant H₀ as an **output** of a self-interacting dark matter (SIDM) Lagrangian extended with a dark QCD axion sector.

The dark axion field σ — a pseudo-Nambu-Goldstone boson of dark SU(2) chiral symmetry breaking — evolves in the potential V(σ) = Λ_d⁴(1 − cos(σ/f)) from reheating to today. The Hubble parameter H₀ = H(t₀) emerges from the numerical solution without being assumed.

**v2** adds RG transmutation analysis, DESI DR1 comparison, gauge coupling unification check, and a derivation of f from first principles — reducing the model to **2 free parameters**.

## Key Results

### H₀ from the ODE solver

With 6 Lagrangian parameters (m_χ, m_φ, α_D, f, Λ_d, θ_i) and observed Ω_χh² = 0.120:

| θ_i (rad) | H₀ (km/s/Mpc) | Ω_DE | w_σ |
|---|---|---|---|
| 3.0 | 71.1 | 0.72 | −0.86 |
| 3.1 | 73.1 | 0.73 | −0.99 |
| π | 73.3 | 0.73 | −1.00 |

### Parameter reduction (v2)

| Analysis | Result | Parameters removed |
|---|---|---|
| **RG transmutation** | Λ_d = μ·exp(−2π/(b₀·α_d(μ))), b₀ = 19/3 | Λ_d → derived from α_d |
| **f from first principles** | f = √(3Ω_DE)·M_Pl / [c·√(1−cosθ_i)], c ~ O(few) | f → derived from {Λ_d, θ_i} |
| **Total** | {f, Λ_d, θ_i} → {α_d, θ_i} | **3 → 2 free DE parameters** |

### DESI DR1 comparison (v2)

| θ_i | w₀ (CPL fit) | wₐ | H₀ | Match |
|---|---|---|---|---|
| 2.887 | **−0.727** | −0.49 | 66.5 | ← exact DESI w₀ |
| 2.92 | −0.793 | −0.37 | 68.0 | Planck+DESI sweet spot |
| π | −1.000 | 0.000 | 73.3 | ΛCDM limit |

DESI DR1 (2024): w₀ = −0.727 ± 0.067 → model predicts θ_i ≈ 2.89.

### Benchmark points (v2)

| Benchmark | θ_i | Λ_d [meV] | H₀ | w₀ |
|---|---|---|---|---|
| Planck + DESI | 2.92 | 2.0 | 68.0 | −0.79 |
| SH0ES + DESI | 2.96 | 2.1 | 73.3 | −0.80 |
| Planck + ΛCDM | π | 1.9 | 68.2 | −1.00 |
| SH0ES + ΛCDM | π | 2.0 | 73.3 | −1.00 |

## Physics

- The same Lagrangian that solves SIDM phenomenology (velocity-dependent cross sections, relic density Ωh² = 0.120, A₄ CP angle θ = 19.47°) also contains a dark axion that drives late-time acceleration.
- Λ_d ~ 2 meV = √(H₀ M_Pl) — the neutrino mass scale coincidence.
- Hilltop quintessence: θ_i ≈ π required for observed H₀.
- σ is a **dark axion** (f >> Λ_d), not a dark pion: f/Λ_d ~ 10²⁹, same hierarchy as f_a/Λ_QCD in QCD.
- Weak Gravity Conjecture satisfied: f = 0.27 M_Pl < M_Pl (with S_inst ~ O(1) at strong coupling).
- Gauge coupling unification with SM does not work (Δ(1/α) = 27.5 at M_GUT) — α_d remains free.

## Files

| File | Description |
|---|---|
| `layer8_cosmic_ode.py` | **Core solver** — coupled Klein-Gordon + Friedmann ODE. Self-contained. |
| `rg_transmutation.py` | **v2**: RG running of α_d → Λ_d transmutation. QCD analogy. Self-contained. |
| `unification_check.py` | **v2**: SM + dark coupling running, unification check at M_GUT. Self-contained. |
| `desi_comparison.py` | **v2**: CPL w₀, wₐ fit extraction, DESI DR1 comparison. Requires `layer8_cosmic_ode.py`. |
| `f_from_first_principles.py` | **v2**: Derivation of f ~ M_Pl from cosmological self-consistency. Requires `layer8_cosmic_ode.py`. |
| `layer8_evolution.png` | Four-panel diagnostic plot: θ(N), H(N), ρ(N), w(N) |
| `README.md` | This file |

## Requirements

```
Python >= 3.9
numpy
scipy
matplotlib  (optional, for plotting)
```

No other dependencies. All scripts are self-contained — run without the parent repository.

## Quick Start

```bash
# Core ODE solver
python layer8_cosmic_ode.py

# RG transmutation (v2)
python rg_transmutation.py

# DESI comparison (v2) — requires layer8_cosmic_ode.py in same directory
python desi_comparison.py

# f from first principles (v2) — requires layer8_cosmic_ode.py in same directory
python f_from_first_principles.py

# Unification check (v2)
python unification_check.py
```

## Usage as Library

```python
from layer8_cosmic_ode import solve_layer8

result = solve_layer8(
    m_chi=98.19,        # GeV
    m_phi_GeV=9.66e-3,  # GeV
    alpha_d=3.274e-3,
    f=0.27 * 2.435e18,  # 0.27 M_Pl in GeV
    Lambda_d=2.0e-12,   # GeV (2 meV)
    theta_i=3.0,        # rad (near π)
    omega_chi_h2=0.120,
)

print(f"H₀ = {result.H0_kms:.2f} km/s/Mpc")
print(f"w_σ = {result.w_sigma:.4f}")
```

## Related

- **Paper 1 (SIDM):** [Secluded-Majorana-SIDM](https://github.com/0mp3s/Secluded-Majorana-SIDM) (DOI: [10.5281/zenodo.19225823](https://doi.org/10.5281/zenodo.19225823)) — Full numerical pipeline for velocity-dependent SIDM cross sections.
- **Paper 2 (in preparation):** Dark energy from the same Lagrangian — Layer 8 connects dark matter and dark energy.

## Citation

```bibtex
@software{layer8_cosmic_ode:2026,
    author  = {Pesach, Omer},
    title   = {{Deriving H₀ from a Secluded Majorana SIDM Lagrangian:
               Coupled Dark Axion–Friedmann Solver}},
    version = {v2},
    year    = {2026},
    publisher = {Zenodo},
    doi     = {10.5281/zenodo.19274637},
    url     = {https://doi.org/10.5281/zenodo.19274637}
}
```

## Changelog

### v2 (2026-03-28)
- **Added** `rg_transmutation.py`: α_d → Λ_d via dimensional transmutation (SU(2)_d, b₀ = 19/3)
- **Added** `unification_check.py`: SM gauge coupling running, crossing analysis at M_GUT
- **Added** `desi_comparison.py`: CPL (w₀, wₐ) fit extraction, DESI DR1 w₀ = −0.727 exact match at θ_i = 2.887
- **Added** `f_from_first_principles.py`: f ~ M_Pl from GMOR + Friedmann + quintessence self-consistency
- **Result**: Parameter reduction from {f, Λ_d, θ_i} to {α_d, θ_i} — 2 free parameters

### v1 (2026-03-28)
- Initial release: `layer8_cosmic_ode.py` solver, diagnostic plot, README

## License

CC-BY-4.0
