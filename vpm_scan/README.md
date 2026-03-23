# vpm_scan/ — VPM Solver Validation

Scripts that validate the core VPM (Variable Phase Method) phase-shift solver.

## Configuration — config.json

The VPM scan grid is configurable via `config.json`:

```jsonc
{
    "grid": {
        "n_chi": 50,
        "n_phi": 70,
        "n_res": 4,
        "n_alpha": 200,
        "m_chi_range_GeV": [0.1, 100.0],
        "m_phi_range_GeV": [0.1e-3, 200e-3],
        "lambda_crits": [1.68, 6.45, 14.7, 26.0]
    },
    "output": {
        "raw_csv": "../data/all_viable_raw_v8.csv",
        "representative_csv": "../data/all_viable_representative_v8.csv"
    }
}
```

## born_validation.py (v21)

Born-limit validation of the VPM solver:
- Verifies s-wave accuracy against analytic Born approximation
- Tests beyond-Born scaling behavior

**Output:** Text summary to stdout

## error_budget.py (v23)

Systematic error budget for the numerical VPM solver:
- Integrator comparison (RK4 vs RK45)
- Partial-wave truncation (varying l_max)
- Short-range prescription (varying x_min)

**Output:** Text summary to stdout
