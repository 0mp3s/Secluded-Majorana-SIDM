# cosmology/ — Mediator Cosmology & Freeze-out

Scripts analyzing the cosmological viability of the scalar mediator φ.

All scripts support `--config custom_config.json` to override default parameters.

## Configuration — config.json

```jsonc
{
    "benchmark": {
        "m_chi_GeV": 42.919,     // DM mass [GeV]
        "m_phi_GeV": 4.233e-3,   // mediator mass [GeV]
        "alpha": 6.172e-4         // dark coupling
    },
    "higgs_portal": {
        "lambda_h_phi": 1e-6     // Higgs-portal coupling for thermalization
    }
}
```

## mediator_cosmology.py (v25)

Checks whether the scalar mediator φ satisfies:
- BBN constraints (must decay before t ~ 1 s)
- ΔN_eff bounds from Planck
- Thermal equilibrium via Higgs portal coupling λ_hφ

Uses `benchmark.m_phi_GeV` from config.json.

**Output:** `output/v25_output.txt`

## bsf_estimate.py (v26)

Bound-state formation (BSF) rate estimate:
- Compares BSF cross section against perturbative annihilation
- Evaluates at freeze-out and dwarf-galaxy velocities

Uses `benchmark.m_chi_GeV`, `benchmark.m_phi_GeV`, `benchmark.alpha` from config.json.

**Output:** Text summary to stdout
