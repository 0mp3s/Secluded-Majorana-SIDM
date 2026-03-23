# stats_mcmc/ — MCMC Posterior Sampling

Bayesian posterior sampling of the (m_χ, m_φ, α) parameter space
using emcee affine-invariant ensemble sampler.

## run_mcmc.py (v38)

*Will be added after the current MCMC run completes.*

MCMC configuration:
- **Sampler:** emcee (32 walkers, 300 burn-in, 2000 production)
- **Parallelism:** multiprocessing pool (12 workers)
- **Priors (uniform):**
  - log₁₀(m_χ/GeV) ∈ [log₁₀(5), log₁₀(200)]
  - log₁₀(m_φ/MeV) ∈ [log₁₀(3), log₁₀(30)]
  - log₁₀(α) ∈ [log₁₀(10⁻⁵), log₁₀(0.05)]
- **Likelihood:** χ² from 13 astrophysical systems (same as chi2_fit.py)
- **Seeds:** 17 relic benchmark points + v34 best-fit point

**Input:** `data/v31_true_viable_points.csv`, core VPM solver  
**Output:** `output/v38_corner.png`, `output/v38_mcmc_chains.png`, 
`output/v38_lambda_posterior.png`, `output/v38_mcmc_samples.csv`
