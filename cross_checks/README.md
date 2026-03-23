# cross_checks/ — Validation, Literature & Theory Cross-Checks

Scripts performing independent validation of the VPM solver and 
physics cross-checks against published results.

All scripts support `--config custom_config.json` to override default parameters.

## Configuration — config.json

Edit `config.json` to change input data or velocities:

```jsonc
{
    "benchmark_csv": "../data/v31_true_viable_points.csv",
    "representative_csv": "../data/all_viable_representative_v8.csv",
    "sommerfeld": {
        "velocities_km_s": [30, 200, 1000, 90000]
    },
    "velocity_averaged": {
        "velocities_km_s": [12, 30, 50, 200, 1000, 1500],
        "n_integration_points": 50
    }
}
```

### Input CSV Format (benchmark_csv)

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| m_chi_GeV | float | GeV | Dark matter mass |
| m_phi_MeV | float | MeV | Mediator mass (⚠ MeV, divided by 1000 internally) |
| alpha | float | — | Dark coupling constant |
| lambda | float | — | Dimensionless coupling α·m_χ/m_φ |

### Input CSV Format (representative_csv)

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| m_chi_GeV | float | GeV | Dark matter mass |
| m_phi_GeV | float | GeV | Mediator mass |
| alpha | float | — | Dark coupling constant |
| sigma_m_30 | float | cm²/g | σ/m at v=30 km/s |
| sigma_m_1000 | float | cm²/g | σ/m at v=1000 km/s |

## blind_sanity.py (v28)

Off-grid validation: 10 interpolated viable + 10 off-grid rejects + 
5 monotonicity checks.

**Input:** `representative_csv`
**Output:** `output/v28_out.txt`

## blind_large.py (v29)

Large-scale validation: 200 random points.

**Input:** `representative_csv`
**Output:** `output/v29_out.txt`

## literature_crosscheck.py (v32)

Comparison against Tulin, Yu & Zurek (2013) PRD 87, 115007.

**Input:** Core VPM solver
**Output:** `output/v32_output.txt`

## tyz_comparison.py (v35)

VPM vs Born transfer cross section comparison.

**Input:** Core VPM solver
**Output:** `output/v35_tyz_comparison.png`

## sommerfeld.py (v36)

Sommerfeld enhancement S₀ for s-wave annihilation via numerical RK4 ODE.
Tests all benchmark points from `benchmark_csv`.

**Input:** `benchmark_csv`
**Output:** `output/v36_sommerfeld.png`, `output/v36_output.txt`

## velocity_averaged.py (v37)

Maxwell-Boltzmann velocity-averaged σ/m vs fixed-velocity comparison.
Uses observational data (configurable, same format as observations/config.json).

**Input:** `benchmark_csv`, observations data
**Output:** `output/v37_velocity_averaged.png`

**Input:** `data/v31_true_viable_points.csv`  
**Output:** `output/v37_velocity_averaged.png`, `output/v37_output.txt`

## Example Output

```
output/
├── v28_out.txt               # Blind sanity check
├── v29_out.txt               # Blind large-scale validation
├── v29_output_fresh.txt      # Fresh re-run
├── v32_output.txt            # Literature crosscheck (6/6 PASSED)
├── v35_tyz_comparison.png    # VPM vs Born figure
├── v36_sommerfeld.png        # Sommerfeld enhancement figure
├── v36_output.txt            # Sommerfeld numerical results
├── v37_velocity_averaged.png # MB averaging figure
└── v37_output.txt            # Velocity averaging results
```
