# Secluded Majorana SIDM

**Self-Interacting Dark Matter via a Secluded Scalar Mediator:
Majorana Fermion with Velocity-Dependent Cross Sections**

A complete numerical pipeline for exploring self-interacting dark matter (SIDM)
in a secluded Majorana model with a light scalar mediator. The code computes
velocity-dependent transfer cross sections, solves the Boltzmann equation for
relic density, and fits the model to astrophysical observations.

---

## Repository Structure

```
Secluded-Majorana-SIDM/
├── README.md                 # This file
├── LICENSE                   # MIT License
├── requirements.txt          # Python dependencies
│
├── core/                     # Shared numerical solvers
│   ├── v22_raw_scan.py       # VPM phase-shift solver (sigma_T_vpm)
│   ├── v27_boltzmann_relic.py # Boltzmann relic density solver
│   └── __init__.py
│
├── data/                     # Shared input/output data
│   ├── all_viable_raw_v8.csv           # 80k SIDM-viable parameter points
│   ├── all_viable_representative_v8.csv # Representative subset
│   ├── v31_true_viable_points.csv      # 17 relic+SIDM viable benchmarks
│   ├── v31_all_relic_points.csv        # All relic-correct points
│   ├── v30_perfect_benchmarks.csv      # Top-5 KT benchmarks
│   └── v34_results.csv                 # chi-squared fit results
│
├── vpm_scan/                 # VPM solver validation
│   ├── born_validation.py    # Born-limit accuracy test
│   └── error_budget.py       # Numerical error analysis
│
├── relic_density/            # Boltzmann solver & benchmark extraction
│   ├── benchmark_extractor.py # Extract top benchmarks from scan
│   ├── boltzmann_correction.py # Verify via exact numerical Boltzmann
│   ├── smart_scan.py          # Full grid scan -> 17 viable BPs
│   └── plot_island.py         # Viability island visualization
│
├── cosmology/                # Mediator cosmology & bound states
│   ├── mediator_cosmology.py  # BBN, Delta-N_eff, Higgs portal
│   └── bsf_estimate.py       # Bound state formation estimate
│
├── observations/             # Astrophysical comparison & chi-squared fit
│   ├── observational_comparison.py  # BP1 vs 13 systems (Figure 3)
│   └── chi2_fit.py           # Full chi-squared fit (multiprocessing)
│
├── cross_checks/             # Independent validation & theory checks
│   ├── blind_sanity.py       # Off-grid blind validation
│   ├── blind_large.py        # Large-scale random validation
│   ├── literature_crosscheck.py # vs Tulin-Yu-Zurek 2013 (6/6 PASSED)
│   ├── tyz_comparison.py     # VPM vs Born transfer cross section
│   ├── sommerfeld.py         # Sommerfeld enhancement computation
│   └── velocity_averaged.py  # Maxwell-Boltzmann averaging
│
├── stats_mcmc/               # MCMC posterior sampling
│   └── run_mcmc.py           # emcee sampler + corner plots (pending)
│
└── docs/                     # Preprint & peer reviews
    ├── preprint_draft_v10.md
    ├── research_journal_v10.md
    └── peer_reviews/
```

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run VPM validation
cd vpm_scan && python born_validation.py

# Run chi-squared fit (uses 12 CPU cores)
cd observations && python chi2_fit.py

# Run MCMC posterior sampling
cd stats_mcmc && python run_mcmc.py
```

## Physics Summary

| Component | Method | Key Result |
|-----------|--------|------------|
| Cross section | VPM partial-wave solver | 80,142 SIDM-viable points |
| Relic density | Numerical Boltzmann (RK4) | 17 benchmark points with Omega h^2 = 0.120 |
| Chi-squared fit | 13 astrophysical systems | chi^2/dof = 0.54 (BP1) |
| Sommerfeld | RK4 Schrodinger solver | S(freeze-out) = 1.003-1.025 (tree-level OK) |
| Velocity averaging | MB integration | Delta(chi^2) < 9% vs fixed-v |
| Literature | vs TYZ (2013) | 6/6 tests PASSED |
| MCMC posterior | emcee (32 walkers) | Pending... |

## Model

Secluded Majorana fermion dark matter (chi) with a light scalar mediator (phi):

- Yukawa potential: V(r) = -alpha * exp(-m_phi * r) / r
- Annihilation: chi chi -> phi phi (s-wave, secluded)
- SIDM: chi chi -> chi chi via phi exchange (velocity-dependent)
- Parameter space: m_chi in [5, 200] GeV, m_phi in [3, 30] MeV, alpha in [1e-5, 0.05]

## Dependencies

- Python >= 3.10
- NumPy, SciPy, Matplotlib
- Numba >= 0.60 (JIT-compiled VPM solver)
- emcee >= 3.1, corner >= 2.2 (MCMC)
- pandas (data handling)

## License

MIT
