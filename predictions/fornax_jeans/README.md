# Fornax dSph — Jeans Equation Analysis

## Purpose
Solve the spherical Jeans equation with the SIDM-cored DM potential
for Fornax, predict σ_los(R), and compare with Walker+2009 stellar
kinematics. **Blind prediction** — no parameter tuning.

## Method
1. NFW halo (Read+2019: M₂₀₀ = 3.16×10⁹ M☉, c = 18)
2. SIDM core via Kaplinghat+2016 isothermal matching
3. Plummer stellar profile (r_half = 710 pc, M_★ = 2×10⁷ M☉)
4. Solve isotropic (β = 0) Jeans equation analytically
5. Abel projection to σ_los(R)

## Data Sources
- **Kinematics**: Walker et al. 2009, ApJ 704, 1274 (2633 member stars)
- **Halo**: Read et al. 2019, MNRAS 484, 1401
- **Core reference**: Walker & Peñarrubia 2011 (r_core ≈ 500–900 pc)

## Key Test
SIDM core → flat σ_los(R) profile (matches data)  
NFW cusp → rising σ_los towards center (too concentrated)

## Output
- `output/fornax_jeans_prediction.png` — two-panel figure:
  left: DM density profile; right: σ_los(R) vs Walker+2009 data
