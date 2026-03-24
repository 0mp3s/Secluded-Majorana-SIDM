# stats_mcmc/ — MCMC Posterior Sampling

Bayesian posterior sampling of the (m_χ, m_φ, α) parameter space
using emcee affine-invariant ensemble sampler.

## run_mcmc.py (v38)

MCMC configuration:
- **Sampler:** emcee (32 walkers, 300 burn-in, 2000 production)
- **Parallelism:** multiprocessing pool (12 workers)
- **Priors (uniform in log₁₀):**
  - log₁₀(m_χ/GeV) ∈ [log₁₀(5), log₁₀(200)]
  - log₁₀(m_φ/MeV) ∈ [log₁₀(3), log₁₀(30)]
  - log₁₀(α) ∈ [log₁₀(10⁻⁵), log₁₀(0.05)]
- **Likelihood:** χ² from 13 astrophysical systems (same as chi2_fit.py)
- **Seeds:** 17 relic benchmark points + v34 best-fit point

**Input:** `data/v31_true_viable_points.csv`, core VPM solver  
**Output:** `output/v38_corner.png`, `output/v38_mcmc_chains.png`, 
`output/v38_lambda_posterior.png`, `output/v38_mcmc_samples.csv`

## Results (run: 24 March 2026, 5000 production steps)

| Metric | Value |
|--------|-------|
| Total samples | 160,000 |
| Acceptance fraction | ~0.52 |
| Max autocorrelation | ~75 steps |
| N/τ | 66.6 (>50, converged) |
| Effective samples | ~2,132 |

### Best-fit (MAP)

| Parameter | Value |
|-----------|-------|
| m_χ | 94.07 GeV |
| m_φ | 11.10 MeV |
| α | 5.734×10⁻³ |
| λ = αm_χ/m_φ | 48.6 |
| χ²/dof | **0.1575** (10 dof) |

### Parameter Constraints (median ± 68% CI)

| Parameter | Median | 16th %ile | 84th %ile |
|-----------|--------|-----------|-----------|
| m_χ [GeV] | 72.98 | 16.97 | 115.91 |
| m_φ [MeV] | 10.83 | 6.51 | 15.59 |
| α | 0.004 | 0.001 | 0.022 |
| λ (derived) | 59.0 | 3.6 | 272.2 |

### Relic BPs vs Posterior

All **17/17** relic benchmark points lie within the 95% credible region.
Best relic BP: BP16 (m_χ=20.7, m_φ=9.91, α=1.048×10⁻³) with χ²=5.396.

## config.json Schema

```json
{
    "mcmc": {
        "n_walkers": 32,        // Number of emcee walkers
        "n_burn": 300,          // Burn-in steps
        "n_production": 2000,   // Production steps
        "n_workers": 12         // Multiprocessing pool size
    },
    "priors": {
        "m_chi_GeV_min": 5.0,   // Lower bound on m_chi [GeV]
        "m_chi_GeV_max": 200.0, // Upper bound on m_chi [GeV]
        "m_phi_MeV_min": 3.0,   // Lower bound on m_phi [MeV]
        "m_phi_MeV_max": 30.0,  // Upper bound on m_phi [MeV]
        "alpha_min": 1e-5,      // Lower bound on alpha
        "alpha_max": 0.05       // Upper bound on alpha
    },
    "relic_bp_csv": "../data/v31_true_viable_points.csv",
    "observations": [...]       // Same format as observations/config.json
}
```

## Output Files

| File | Description |
|------|-------------|
| `output/v38_corner.png` | 3×3 corner plot with 68%/95% contours |
| `output/v38_mcmc_chains.png` | Walker trace plots for convergence check |
| `output/v38_lambda_posterior.png` | Posterior density of derived λ parameter |
| `output/v38_mcmc_samples.csv` | Full MCMC chain (64,000 rows × 8 columns) |
