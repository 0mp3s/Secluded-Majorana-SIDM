# vpm_scan/ — VPM Solver Validation

Scripts that validate the core VPM (Variable Phase Method) phase-shift solver.

## born_validation.py (v21)

Born-limit validation of the VPM solver:
- Verifies s-wave accuracy against analytic Born approximation
- Tests beyond-Born scaling behavior
- Checks barrier-physics regime

**Input:** None (self-contained benchmarks)  
**Output:** Text summary to stdout

## error_budget.py (v23)

Systematic error budget for the numerical VPM solver:
- Integrator comparison (RK4 vs RK45)
- Partial-wave truncation (varying l_max)
- Short-range prescription (varying x_min)

Tests 3 benchmark points × 2 velocities.

**Input:** None (self-contained benchmarks)  
**Output:** Text summary to stdout

## Example Output

See `output/` directory for saved run outputs (currently empty — 
these scripts output to stdout, redirect with `> output/born_validation.txt`).
