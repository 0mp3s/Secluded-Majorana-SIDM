# predictions/ — Testable Predictions

Quantitative predictions derived from the MCMC posterior (§4.7) and the
17 relic benchmark points, tested against published observational data.

Each sub-folder contains:
- **data/** or a CSV with published observational measurements
- **predict_*.py** — script that computes our model prediction and compares to data
- **config.json** — configurable parameters
- **README.md** — description, data sources, expected outcome

## Predictions

| # | Prediction | Data Source | Key Test |
|---|-----------|-------------|----------|
| 1 | [Gravothermal collapse timescales in dSphs](gravothermal/) | Read+2019, Walker+2009 (GAIA DR3) | λ~333 (MAP) predicts collapse; λ~4 (BP1) predicts stable cores |
| 2 | [SPARC rotation curve diversity](rotation_curves/) | Lelli+2016 (SPARC), Oman+2015 | Core size predictions vs measured inner slopes |
| 3 | [DM–galaxy offsets in cluster mergers](cluster_offsets/) | Harvey+2015 (72 clusters) | σ/m(v_cluster) → predicted offset vs upper bound |
| 4 | [Extra radiation ΔN_eff](delta_neff/) | Planck 2018, CMB-S4 forecast | ΔN_eff ≈ 0 (m_φ ≫ T at BBN/CMB → Boltzmann-suppressed) |
| 5 | [Fornax stellar dispersion (Jeans)](fornax_jeans/) | Walker+2009 (2633 stars) | SIDM core flattens σ_los(R); NFW cusp ruled out |

## Usage

```bash
# Run all predictions
cd predictions/gravothermal && python predict_gravothermal.py
cd ../rotation_curves && python predict_core_sizes.py
cd ../cluster_offsets && python predict_offsets.py
cd ../delta_neff && python predict_neff.py
cd ../fornax_jeans && python predict_fornax_jeans.py
```

## Key Insight

The MCMC posterior (v38) shows that **without** the relic-density constraint,
astrophysical data prefer higher masses (m_χ ~ 70–90 GeV) with stronger
coupling (α ~ 0.02). This leads to λ ~ 60–330, deep in the resonant regime.
The relic-constrained BPs have λ ~ 2–65.

These predictions discriminate between the two regimes using independent data.
