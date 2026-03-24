# gravothermal/ — Gravothermal Collapse Timescales in dSphs

## Prediction

For each classical dwarf spheroidal (dSph), we compute the SIDM scattering
rate and check whether the halo is in the core-formation, stable-core, or
gravothermal-collapse regime.

**Key test:** The MAP benchmark (λ ≈ 48.6) sits in the core-formation regime
for most classical dSphs — **consistent** with their observed flat density cores.
BP1 (λ~4) predicts stable cores throughout.

## Physics

The scattering rate per particle:
```
Γ = (σ/m) × ρ_s × v_rel
```

The number of scatterings over the halo's lifetime:
```
N_scatter = Γ × t_age
```

| N_scatter | Phase | Profile |
|-----------|-------|---------|
| < 1 | Too few scatterings | NFW cusp survives |
| 1–100 | Core formation | Flat density core (observed) |
| > 100–300 | Gravothermal collapse | Re-cusping begins |

References: Kaplinghat+2016 (PRL 116, 041302), Essig+2019 (PRL 123, 121102),
Nishikawa+2020 (PRD 101, 063009)

## Data

**dsphs_data.csv** — 8 classical Milky Way dSphs:

| Column | Units | Source |
|--------|-------|--------|
| `sigma_v_km_s` | km/s | Walker+2009 (velocity dispersions) |
| `r_half_pc` | pc | McConnachie 2012 |
| `rho_s_Msun_kpc3` | M☉/kpc³ | Read+2019 (Jeans modeling, GAIA proper motions) |
| `r_s_kpc` | kpc | Read+2019 (NFW scale radius) |
| `core_observed` | YES/AMBIGUOUS | Read+2018 (core/cusp classification) |
| `t_age_Gyr` | Gyr | SFH studies |

## Usage

```bash
python predict_gravothermal.py
python predict_gravothermal.py --config custom.json
```

## Expected Output

- Table: N_scatter per dSph × benchmark point, with core/collapse classification
- Plot: `output/gravothermal_prediction.png` — bar chart of N_scatter with thresholds
- Verdict: which benchmarks are consistent with observed cores
