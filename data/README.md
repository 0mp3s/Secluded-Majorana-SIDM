# data/ — Shared Data Files

## Raw Scan Output

| File | Rows | Description |
|------|------|-------------|
| `all_viable_raw_v8.csv` | 80,142 | Every SIDM-viable (m_χ, m_φ, α) from VPM scan |
| `all_viable_representative_v8.csv` | ~600 | 1 representative per (m_φ, resonance) cell |

**Columns:** `m_chi_GeV, m_phi_GeV, alpha, sigma_m_30, sigma_m_1000, resonance_idx`

## Relic-Constrained Points

| File | Rows | Description |
|------|------|-------------|
| `v31_true_viable_points.csv` | 17 | Relic + SIDM viable benchmark points |
| `v31_all_relic_points.csv` | ~600 | All relic-correct points (SIDM viable + non-viable) |
| `v30_perfect_benchmarks.csv` | 5 | Top-5 KT-approximation benchmarks |

**Columns (v31):** `m_chi_GeV, m_phi_MeV, alpha, omega_h2, lambda, sigma_m_30, sigma_m_1000`

⚠ **Unit note:** `m_phi_MeV` in v31 files is in MeV. The VPM solver `sigma_T_vpm()` expects m_φ in **GeV**. Divide by 1000.

## χ² Fit Results

| File | Rows | Description |
|------|------|-------------|
| `v34_results.csv` | 5,026 | χ² fit results for all scanned parameter points |
