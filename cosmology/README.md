# cosmology/ — Mediator Cosmology & Freeze-out

Scripts analyzing the cosmological viability of the scalar mediator φ.

## mediator_cosmology.py (v25)

Checks whether the scalar mediator φ satisfies:
- BBN constraints (must decay before t ~ 1 s)
- ΔN_eff bounds from Planck
- Thermal equilibrium via Higgs portal coupling λ_hφ

**Input:** None (self-contained parameter tables)  
**Output:** `output/v25_output.txt`

## bsf_estimate.py (v26)

Bound-state formation (BSF) rate estimate:
- Compares BSF cross section against perturbative annihilation
- Evaluates at freeze-out and dwarf-galaxy velocities
- Determines whether BSF correction is needed

**Input:** None (self-contained benchmarks)  
**Output:** Text summary to stdout

## Example Output

```
output/
├── v25_output.txt    # Mediator cosmology results
```
