# data/ — Shared Data Files

All data files are CSV with headers. Scripts resolve paths via `DATA_DIR` or `config.json`.

## Raw Scan Output

| File | Rows | Description |
|------|------|-------------|
| `all_viable_raw_v8.csv` | 80,142 | Every SIDM-viable (m_χ, m_φ, α) from VPM scan |
| `all_viable_representative_v8.csv` | ~600 | 1 representative per (m_φ, resonance) cell |

**Columns (raw scan):**

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| m_chi_GeV | float | GeV | Majorana DM mass |
| m_phi_GeV | float | GeV | Scalar mediator mass |
| alpha | float | — | Dark-sector coupling α_D = g²/(4π) |
| sigma_m_30 | float | cm²/g | Transfer cross section σ_T/m_χ at v=30 km/s |
| sigma_m_1000 | float | cm²/g | Transfer cross section σ_T/m_χ at v=1000 km/s |
| resonance_idx | int | — | Resonance band index (0–3, mapping to λ_crit) |

## Relic-Constrained Points

| File | Rows | Description |
|------|------|-------------|
| `v31_true_viable_points.csv` | 17 | Relic + SIDM viable benchmark points |
| `v31_all_relic_points.csv` | ~600 | All relic-correct points (SIDM viable + non-viable) |
| `v30_perfect_benchmarks.csv` | 5 | Top-5 KT-approximation benchmarks (early scan, before bisection) |

**Columns (v31 relic CSVs):**

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| m_chi_GeV | float | GeV | DM mass |
| m_phi_MeV | float | **MeV** | Mediator mass ⚠ **in MeV** |
| alpha | float | — | Coupling (from bisection on Ωh²) |
| omega_h2 | float | — | Numerical relic density |
| lambda | float | — | Dimensionless coupling α·m_χ/m_φ |
| sigma_m_30 | float | cm²/g | σ_T/m at v=30 km/s |
| sigma_m_1000 | float | cm²/g | σ_T/m at v=1000 km/s |

⚠ **Unit warning:** `m_phi_MeV` in `v31_true_viable_points.csv` is in **MeV**. The VPM solver `sigma_T_vpm()` expects m_φ in **GeV**. Divide by 1000 before passing to the solver.

⚠ **Header difference:** `v31_all_relic_points.csv` uses `m_phi_GeV` (values in **GeV**), while `v31_true_viable_points.csv` uses `m_phi_MeV` (values in **MeV**). Code must handle both unit conventions.

**Extra columns in `v31_all_relic_points.csv`:**

| Column | Type | Description |
|--------|------|-------------|
| sidm_viable | int | 1 if σ/m falls within SIDM window at v = 30 km/s, 0 otherwise |

**Columns in `v30_perfect_benchmarks.csv`:**

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| m_chi_GeV | float | GeV | DM mass |
| m_phi_GeV | float | GeV | Mediator mass (in GeV, unlike v31 files) |
| alpha | float | — | Coupling |
| omega_h2 | float | — | KT-approximation relic density |
| sigma_m_5 | float | cm²/g | σ_T/m at v = 5 km/s |
| sigma_m_10 | float | cm²/g | σ_T/m at v = 10 km/s |
| sigma_m_30 | float | cm²/g | σ_T/m at v = 30 km/s |
| sigma_m_1000 | float | cm²/g | σ_T/m at v = 1000 km/s |
| resonance_idx | int | — | Resonance band index |
| delta_omega | float | — | |Ωh² − 0.120| deviation |

## χ² Fit Results

| File | Rows | Description |
|------|------|-------------|
| `v34_results.csv` | 5,026 | χ² fit results for all scanned parameter points |

**Columns (v34):**

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| chi2 | float | — | Total χ² |
| chi2_dof | float | — | χ² per degree of freedom |
| m_chi_GeV | float | GeV | DM mass |
| m_phi_MeV | float | MeV | Mediator mass |
| alpha | float | — | Coupling |
| lambda | float | — | Dimensionless coupling |
| is_relic | int | — | 1 if relic-constrained BP |
| omega_h2 | float | — | Relic density (relic BPs only) |
| sigma_m_V | float | cm²/g | σ_T/m at each observed velocity V |
