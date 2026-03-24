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
│   ├── run_mcmc.py           # emcee sampler + corner plots
│   └── config.json           # MCMC configuration
│
├── mixed_coupling/           # Mixed scalar/pseudoscalar coupling analysis
│   └── opusB/                # Condition validation (C1–C4)
│
├── predictions/              # Testable predictions vs published data
│   ├── gravothermal/         # dSph gravothermal evolution
│   ├── rotation_curves/      # SPARC diversity (V(2kpc) vs V_max)
│   ├── cluster_offsets/      # Cluster merger σ/m bounds
│   ├── delta_neff/           # ΔN_eff from light mediator
│   ├── fornax_jeans/         # Fornax Jeans analysis (SIDM + feedback + β)
│   ├── multi_dsph_jeans/     # Multi-dSph cross-validation (5 dSphs)
│   └── majorana_vs_dirac/    # Majorana vs Dirac σ_T(v) fingerprint
│
├── discussion/               # A–B discussion notes
│   └── MixedMajorana.md
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

# Run testable predictions
cd predictions/gravothermal && python predict_gravothermal.py
cd ../rotation_curves && python predict_core_sizes.py
cd ../cluster_offsets && python predict_offsets.py
cd ../delta_neff && python predict_neff.py

# Fornax Jeans analysis (SIDM + feedback + Osipkov-Merritt anisotropy)
cd ../fornax_jeans && python predict_fornax_jeans_aniso.py

# Multi-dSph cross-validation (Fornax, Sculptor, Draco, Carina, Sextans)
cd ../multi_dsph_jeans && python predict_multi_dsph_jeans.py

# Majorana vs Dirac σ_T(v) fingerprint
cd ../majorana_vs_dirac && python predict_maj_vs_dir.py
```

## Configuration

All scripts support external configuration via JSON files, so you can run them
with **new data** without editing source code.

Each folder has a `config.json` with the default parameter values.
Override them by editing the file or using the `--config` CLI flag:

```bash
# Use custom config
python observations/chi2_fit.py --config my_observations.json

# Or edit the default config in-place
nano observations/config.json
```

Key configurable parameters per folder:
- **observations/**: observational data (systems, velocities, σ/m bounds), input CSVs, worker count
- **relic_density/**: scan grid (m_χ, m_φ ranges, resolution), Boltzmann solver settings, SIDM cuts
- **cross_checks/**: benchmark CSV path, test velocities
- **cosmology/**: benchmark point (m_χ, m_φ, α), Higgs-portal coupling
- **core/**: VPM scan grid (when run as standalone script)

See each folder's `README.md` for the full config.json schema and expected data formats.

## Physics Summary

| Component | Method | Key Result |
|-----------|--------|------------|
| Cross section | VPM partial-wave solver | 80,142 SIDM-viable points |
| Relic density | Numerical Boltzmann (RK4) | 17 benchmark points with Omega h^2 = 0.120 |
| Chi-squared fit | 13 astrophysical systems | chi^2/dof = 0.54 (BP1) |
| Sommerfeld | RK4 Schrodinger solver | S(freeze-out) = 1.003-1.025 (tree-level OK) |
| Velocity averaging | MB integration | Delta(chi^2) < 9% vs fixed-v |
| Literature | vs TYZ (2013) | 6/6 tests PASSED |
| MCMC posterior | emcee (32 walkers, 64k samples) | chi^2/dof = 0.16 (MAP), 17/17 BPs within 95% CI |
| **Gravothermal** | dSph core formation | MAP: 6/8 CORED (Fornax, Sculptor, Carina, Sextans, Leo I, Leo II) |
| **Cluster mergers** | Harvey+2015 bounds | ALL PASS (8 systems × 3 BPs) |
| **SPARC rotation** | Baryons + SIDM fit | DDO_154: Υ_* = 0.32 (physical range) |
| **ΔN_eff** | Light mediator BBN | ≈ 0 (Boltzmann-suppressed) |
| **Fermi-LAT** | dSph indirect detection | Tree-level = 0 (secluded); loop 16 orders below UL |
| **Fornax Jeans** | σ_los(R) profile fit | MAP+fb+β: **χ²/dof = 1.7/7** (all residuals < 1σ) |
| **Multi-dSph** | 5 classical dSphs | **4/5 with χ²/dof < 2** — not a statistical fluke |
| **Majorana vs Dirac** | σ_T(v) ratio comparison | Ratio oscillates 0.5–1.13: smoking-gun Majorana fingerprint |

### Peer Review Status

Four peer reviews received (Claude Opus 4.6, GPT-5.4 ×2, Gemini 3.1 Pro).
Technical bugs fixed (NFW ρ_s normalization, CSV columns, λ convention).
Two physics issues identified (s-wave for Majorana; φ overclosure)
— resolved via **mixed scalar/pseudoscalar Yukawa coupling** (§7.1):
$\mathcal{L} \supset \frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\phi$.
Relic density depends on $\alpha_s \times \alpha_p$ (s-wave annihilation),
SIDM depends on $\alpha_s$ only (elastic scattering). See `mixed_coupling/` for validation.

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
