# observations/ — Astrophysical Comparison & χ² Fitting

Scripts that compare model predictions against astrophysical observations
of self-interacting dark matter.

All scripts support `--config custom_config.json` to override default parameters.

## Configuration — config.json

Edit `config.json` to use different observational data or input CSVs:

```jsonc
{
    "observations": [
        // Each entry: [name, v_km_s, sigma_m_central, sigma_m_lo, sigma_m_hi, reference]
        ["Draco dSph", 12, 0.6, 0.1, 2.0, "KTY16"],
        // ... add/remove systems as needed
    ],
    "scan_data_csv":  "../data/all_viable_raw_v8.csv",   // chi2_fit input
    "relic_bp_csv":   "../data/v31_true_viable_points.csv", // relic BPs input
    "n_workers": 12   // parallel workers for chi2_fit.py
}
```

### Observations Format

Each observation entry is a 6-element list:

| Index | Field | Type | Unit | Description |
|-------|-------|------|------|-------------|
| 0 | name | string | — | System identifier |
| 1 | v_km_s | number | km/s | Characteristic velocity dispersion |
| 2 | sigma_m_central | number | cm²/g | Central σ/m measurement |
| 3 | sigma_m_lo | number | cm²/g | Lower 1σ bound |
| 4 | sigma_m_hi | number | cm²/g | Upper 1σ bound |
| 5 | reference | string | — | Literature reference key |

### Input CSV Formats

**scan_data_csv** — Raw VPM scan results. Required columns:
- `m_chi_GeV` (float): Dark matter mass in GeV
- `m_phi_GeV` (float): Mediator mass in GeV
- `alpha` (float): Dark coupling constant

**relic_bp_csv** — Relic-viable benchmark points. Required columns:
- `m_chi_GeV` (float): Dark matter mass in GeV
- `m_phi_MeV` (float): Mediator mass in **MeV** (⚠ divided by 1000 internally)
- `alpha` (float): Dark coupling constant
- `omega_h2` (float): Relic density

## observational_comparison.py (v33)

Overlays benchmark σ/m(v) curves on 13 real astrophysical constraints.

**Input:** Core VPM solver (via `from v22_raw_scan import sigma_T_vpm`)
**Output:** `output/v33_observational_comparison.png`, `output/v33_output.txt`

Benchmarks are also configurable via config.json:
```json
"benchmarks": [["BP1", 20.69, 0.01134, 0.001048], ...]
```
Format: `[name, m_chi_GeV, m_phi_GeV, alpha]`

## chi2_fit.py (v34)

Full χ² fit of the Secluded Majorana SIDM model to observational systems:
- Uses multiprocessing for parallel VPM evaluation
- Scans parameter points × velocities
- Computes χ²/dof for free best-fit and relic-constrained benchmark points

**Input:** `scan_data_csv`, `relic_bp_csv` (see config.json)
**Output:** `output/v34_chi2_fit.png`, `output/v34_output.txt`, `output/v34_results.csv`

## Example Output

```
output/
├── v33_observational_comparison.png  # Figure: BP curves vs observations
├── v33_output.txt                     # Numerical comparison table
├── v34_chi2_fit.png                   # Figure: χ² landscape
├── v34_output.txt                     # Fit results summary
```
