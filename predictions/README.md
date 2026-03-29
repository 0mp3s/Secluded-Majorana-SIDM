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
| 1 | [Gravothermal collapse timescales in dSphs](gravothermal/) | Read+2019, Walker+2009 (GAIA DR3) | λ~49 (MAP) → core formation; λ~4 (BP1) → stable cores |
| 2 | [SPARC rotation curve diversity](rotation_curves/) | Lelli+2016 (SPARC), Oman+2015 | Core size predictions vs measured inner slopes |
| 3 | [DM–galaxy offsets in cluster mergers](cluster_offsets/) | Harvey+2015 (72 clusters) | σ/m(v_cluster) → predicted offset vs upper bound |
| 4 | [Extra radiation ΔN_eff](delta_neff/) | Planck 2018, CMB-S4 forecast | ΔN_eff ≈ 0 (m_φ ≫ T at BBN/CMB → Boltzmann-suppressed) |
| 5 | [Fornax stellar dispersion (Jeans)](fornax_jeans/) | Walker+2009 (2633 stars) | SIDM core flattens σ_los(R); NFW cusp ruled out |
| 6 | **DM mass range: m_χ ∈ [100, 350] GeV preferred** | G8f scan (30 Mar 2026) | Dual-constraint (SIDM + transmutation) peaks at m_χ = 119–143 GeV; best mismatch 0.0% at 119.6 GeV vs MAP 3.9%. Model predicts m_χ > 59 GeV (relic) and prefers ~120–143 GeV (coincidence). Falsifiable by: future SIDM MCMC with extended prior m_χ > 100 GeV. |
| 7 | **Mixing angle θ from A₄ symmetry: tanθ = 1/3** | G8e derivation (29 Mar 2026) | tanθ = g_p/g_s = 1/3 is derived (not assumed) from A₄ CG coefficients + VEV ratio 3/(2√2). Implies sin²θ = 1/10 (or 1/9 with full VEV correction). Falsifiable by: dark SU(2)_d lattice simulation of the scalar potential minimum. |
| 8 | **Λ_d = m_ν universal for all m_χ ∈ [100, 500] GeV** | G8f scan (30 Mar 2026) | Transmutation chain gives Λ_d = 3.031 meV ≡ m_ν(normal hierarchy) for 100% of 6,848 viable SIDM points, independent of m_χ. This is a structural prediction: any SIDM-viable Secluded-Majorana model with dark QCD confinement must satisfy Λ_d = m_ν via the seesaw matching condition. |

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
coupling (α ~ 0.02). The relic-constrained posterior has MAP λ ≈ 48.6,
with the 17 benchmark points spanning λ ~ 2–65.

These predictions discriminate between the two regimes using independent data.

## Theoretical Predictions (G8 investigation, 29–30 Mar 2026)

Predictions 6–8 are **structural predictions** from the G8 gap investigation —
they follow from first principles of the model and are not fits to data.

### P6: DM mass m_χ ∈ [100, 350] GeV preferred (G8f)

From a full parameter scan (200 grains, 6,848 viable SIDM points):
- Best dual-constraint point: **m_χ = 119.6 GeV**, mismatch = **0.0%**
- MAP (98 GeV) has mismatch = 3.9% — not the global minimum
- Distribution peaks at **143 GeV**, declines sharply above ~350 GeV
- Viable range: m_χ ∈ [100, ~500] GeV with soft cutoff at ~350 GeV

**Falsifiable:** MCMC with extended prior m_χ ∈ [50, 500] GeV should find
a new best-fit near 119–143 GeV, NOT at the old MAP of 98 GeV.

### P7: Mixing angle tanθ = 1/3 derived from A₄ (G8e)

The Yukawa mixing angle θ (sin²θ = 1/10) is NOT a free parameter.
It is **derived** from A₄ × SU(2)_d symmetry:
- CG coefficient ratio: g_p/g_s = 1 (bare A₄)
- VEV correction v_p/v_s = 3/(2√2) from A₄ scalar potential minimum
- Combined: **tanθ = 1/3**, sin²θ = 1/10 (or 1/9 with full VEV)

**Falsifiable:** Dark SU(2)_d lattice simulation of (1+1') scalar sector.
If V_min gives v_p/v_s ≠ 3/(2√2), the prediction shifts.

### P8: Λ_d = m_ν universal for all m_χ (G8f + G8a)

The dark QCD confinement scale equals the neutrino mass scale for
**all** m_χ ∈ [100, 500] GeV — not a coincidence at MAP only.

- Λ_d = 3.031 meV ≡ m_ν(normal hierarchy, lightest) for **100%** of 6,848 points
- Independent of m_χ — this is the transmutation chain structure
- The chain: α_d = 2π/(b₀ ln(m_χ M_R/v²)) → Λ_d = m_χ exp(-2π/(b₀ α_d)) = m_ν

**Falsifiable:** If neutrino oscillation experiments establish m_ν ≠ 3 meV
(e.g., inverted hierarchy with m_lightest > 10 meV), the model is ruled out.
Normal hierarchy (m_lightest < 0.1 eV, lightest ≈ 0) is required.
