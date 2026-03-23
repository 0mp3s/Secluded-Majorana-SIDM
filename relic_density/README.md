# relic_density/ — Relic Density & Benchmark Extraction

Scripts that solve the Boltzmann equation numerically to find parameter
points with Ω_χ h² = 0.120 ± 0.001 (Planck 2018) while satisfying SIDM constraints.

## benchmark_extractor.py (v30)

Extracts top benchmark points from the raw VPM scan:
1. Loads `data/all_viable_raw_v8.csv` (80k SIDM-passing points)
2. Applies Kolb-Turner relic filter (0.115 ≤ Ωh² ≤ 0.125)
3. Ranks by proximity to Planck central value
4. Saves top results

**Input:** `data/all_viable_raw_v8.csv`  
**Output:** `data/v30_perfect_benchmarks.csv`, stdout

## boltzmann_correction.py (v30)

Verifies benchmark points via exact numerical Boltzmann (not KT approximation).
Uses bisection on α to find exact Ωh² = 0.120.

**Input:** Benchmark points (hardcoded from v30 results)  
**Output:** stdout

## smart_scan.py (v31)

Full-grid numerical Boltzmann scan:
- For each (m_χ, m_φ) cell, bisects α to find exact Ωh² = 0.120
- Checks SIDM viability at each solution
- Produces the final 17 viable benchmark points

**Input:** None (scans parameter space directly)  
**Output:** `data/v31_true_viable_points.csv`, `data/v31_all_relic_points.csv`

## plot_island.py (v31)

Visualizes the "island of viability" in the (m_χ, m_φ) plane:
- σ/m heatmap with SIDM-viable region overlay
- BP1 velocity-dependent cross section profile

**Input:** `data/v31_all_relic_points.csv`  
**Output:** `output/v31_island_of_viability.png`, `output/v31_bp1_velocity_profile.png`

## Example Output

```
output/
├── v27_output.txt              # Boltzmann solver verification
├── v30_output.txt              # Benchmark extraction summary
├── v31_output.txt              # Smart scan results
├── v31_island_of_viability.png # Figure: parameter island
└── v31_bp1_velocity_profile.png # Figure: BP1 σ/m(v) curve
```
