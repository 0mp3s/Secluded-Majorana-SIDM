# relic_density/ — Relic Density & Benchmark Extraction

Scripts that solve the Boltzmann equation numerically to find parameter
points with Ω_χ h² = 0.120 ± 0.001 (Planck 2018) while satisfying SIDM constraints.

All scripts support `--config custom_config.json` to override default parameters.

## Configuration — config.json

Edit `config.json` to change the scan grid, Boltzmann solver parameters, or SIDM cuts:

```jsonc
{
    "grid": {
        "n_chi": 20,                        // number of m_χ grid points
        "n_phi": 30,                        // number of m_φ grid points
        "m_chi_range_GeV": [10.0, 100.0],   // DM mass range [GeV]
        "m_phi_range_GeV": [1e-3, 50e-3]    // mediator mass range [GeV]
    },
    "boltzmann": {
        "target_omega_h2": 0.1200,          // target relic density
        "bisect_rtol": 1e-4,                // bisection relative tolerance
        "bisect_max_iter": 50,              // maximum bisection iterations
        "alpha_range": [1e-5, 0.05]         // coupling bracket [lo, hi]
    },
    "sidm_cuts": {
        "sigma_m_30_lo": 1.0,               // σ/m(30 km/s) lower bound [cm²/g]
        "sigma_m_30_hi": 10.0,              // σ/m(30 km/s) upper bound [cm²/g]
        "sigma_m_1000_hi": 0.1              // σ/m(1000 km/s) upper bound [cm²/g]
    },
    "output": {
        "all_relic_csv": "../data/v31_all_relic_points.csv",
        "viable_csv": "../data/v31_true_viable_points.csv"
    }
}
```

### Output CSV Formats

**all_relic_csv** — All points that converged (viable or not):

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| m_chi_GeV | float | GeV | Dark matter mass |
| m_phi_GeV | float | GeV | Mediator mass |
| alpha | float | — | Dark coupling |
| omega_h2 | float | — | Relic density |
| lambda | float | — | Dimensionless coupling α·m_χ/m_φ |
| sigma_m_30 | float | cm²/g | σ/m at v=30 km/s |
| sigma_m_1000 | float | cm²/g | σ/m at v=1000 km/s |
| sidm_viable | int | — | 1 if passes SIDM cuts, 0 otherwise |

**viable_csv** — Same columns minus `sidm_viable` (all are viable).

## smart_scan.py (v31)

Full-grid numerical Boltzmann scan:
- For each (m_χ, m_φ) cell, bisects α to find exact Ωh² = target
- Checks SIDM viability at each solution

**Input:** None (scans parameter space directly)
**Output:** CSVs as configured above

## benchmark_extractor.py (v30)

Extracts top benchmark points from the raw VPM scan.

**Input:** `data/all_viable_raw_v8.csv`
**Output:** `data/v30_perfect_benchmarks.csv`

## boltzmann_correction.py (v30)

Verifies benchmark points via exact numerical Boltzmann.

## plot_island.py (v31)

Visualizes the "island of viability" in the (m_χ, m_φ) plane.

**Input:** `data/v31_all_relic_points.csv`
**Output:** `output/v31_island_of_viability.png`, `output/v31_bp1_velocity_profile.png`
