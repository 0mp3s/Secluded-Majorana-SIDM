# observations/ — Astrophysical Comparison & χ² Fitting

Scripts that compare model predictions against astrophysical observations
of self-interacting dark matter.

## observational_comparison.py (v33)

Overlays BP1 σ/m(v) curve on 13 real astrophysical constraints:
- Dwarf galaxies (Kaplinghat+16): v ~ 30–50 km/s, σ/m ~ 1–10 cm²/g
- Galaxy clusters (Kaplinghat+16): v ~ 1000–1500 km/s, σ/m ~ 0.1 cm²/g
- Bullet Cluster (Markevitch+04): σ/m < 1 cm²/g at v ~ 4700 km/s
- MW satellites (Kamada+17): v ~ 20–100 km/s

**Input:** Core VPM solver (via `from v22_raw_scan import sigma_T_vpm`)  
**Output:** `output/v33_observational_comparison.png`, `output/v33_output.txt`

## chi2_fit.py (v34)

Full χ² fit of the Secluded Majorana SIDM model to 13 observational systems:
- Uses multiprocessing (12 workers) for parallel VPM evaluation
- Scans 5,026 parameter points × 12 velocities = 60,312 VPM calls
- Computes χ²/dof for free best-fit and relic-constrained benchmark points

**Input:** `data/all_viable_raw_v8.csv`, `data/v31_true_viable_points.csv`  
**Output:** `output/v34_chi2_fit.png`, `output/v34_output.txt`, `data/v34_results.csv`

**Key results:**
- Free best-fit: χ²/dof = 0.26
- Relic BP1: χ²/dof = 0.54
- All 17 relic BPs: χ²/dof < 0.85

## Example Output

```
output/
├── v33_observational_comparison.png  # Figure: BP1 vs observations
├── v33_output.txt                     # Numerical comparison table
├── v34_chi2_fit.png                   # Figure: χ² landscape
├── v34_output.txt                     # Fit results summary
```
