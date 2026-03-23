# cross_checks/ — Validation, Literature & Theory Cross-Checks

Scripts performing independent validation of the VPM solver and 
physics cross-checks against published results.

## blind_sanity.py (v28)

Off-grid validation: 10 interpolated viable + 10 off-grid rejects + 
5 monotonicity checks. Verifies the VPM scan boundaries are correct.

**Input:** `data/all_viable_representative_v8.csv`  
**Output:** `output/v28_out.txt`

## blind_large.py (v29)

Large-scale validation: 200 random points (100 in viable region, 
50 full space, 50 perturbed from viable edges).

**Input:** `data/all_viable_representative_v8.csv`  
**Output:** `output/v29_out.txt`, `output/v29_output_fresh.txt`

## literature_crosscheck.py (v32)

Comparison against Tulin, Yu & Zurek (2013) PRD 87, 115007:
- 6 published benchmark tests
- VPM vs TYZ analytical approximations
- All 6/6 PASSED

**Input:** Core VPM solver  
**Output:** `output/v32_output.txt`

## tyz_comparison.py (v35)

VPM vs Born transfer cross section comparison:
- Born regime: VPM/Born ≈ 0.8–0.9 (expected ≤ 1)
- Beyond-Born: Born fails at λ > 1 (as expected)
- Demonstrates necessity of full VPM computation

**Input:** Core VPM solver  
**Output:** `output/v35_tyz_comparison.png`

## sommerfeld.py (v36)

Sommerfeld enhancement S₀ for s-wave annihilation:
- Full numerical RK4 solution of l=0 radial Schrödinger equation
- S(freeze-out) = 1.003–1.025 → tree-level annihilation valid
- S(30 km/s) = 5.7–417 → irrelevant (late-time, not at freeze-out)

Tests all 17 relic-viable benchmark points.

**Input:** `data/v31_true_viable_points.csv`  
**Output:** `output/v36_sommerfeld.png`, `output/v36_output.txt`

## velocity_averaged.py (v37)

Maxwell-Boltzmann velocity-averaged σ/m vs fixed-velocity comparison:
- Per-system shift ≤ 17%, χ² shift 2–9% (mean 6%)
- All 17/17 relic BPs remain viable after averaging
- Validates the fixed-velocity approximation used in the χ² fit

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
