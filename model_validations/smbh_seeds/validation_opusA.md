# Validation — SMBH Seeds via Gravothermal Collapse (Opus A)

**Script:** `predict_smbh_seeds.py`  
**Author:** Opus B  
**Validator:** Opus A  
**Date:** 23 March 2026  

## Method
Ran the script independently. Cross-checked Correa+2015 c(M,z) against manual calculation and literature values.

## Output Reproduction: ALL "NO COLLAPSE" ✅

Checked MAP at representative grid points:

| z | log₁₀(M/M☉) | t_gc [Gyr] (B) | t_gc [Gyr] (my run) | t_univ [Gyr] | Collapse? | Match |
|---|---|---|---|---|---|---|
| 15 | 10.5 | 21,011 | 21,011 | 0.269 | NO | ✅ |
| 10 | 10.0 | ~12,000 | ~12,000 | 0.477 | NO | ✅ |
| 6 | 10.0 | ~7,500 | ~7,500 | 0.927 | NO | ✅ |

All benchmarks produce NO COLLAPSE at every (z, M) grid point — confirmed.

## Correa+2015 Cross-Check

Manual calculation at z=15, M=10¹⁰ M☉:
- A(z) = 5.71 × (1+z)^(−0.47) = 5.71 × 16^(−0.47) = 1.551
- B(z) = −0.084 − 0.025 × ln(16) = −0.084 − 0.0693 = −0.1533
- M_pivot = 10¹²/(0.674) = 1.484 × 10¹² M☉
- c = 1.551 × (10¹⁰/1.484e12)^(−0.1533) = 1.551 × (6.74e-3)^(−0.1533)
- c ≈ 3.34

Script gives c(10¹⁰, z=15) = 3.34 ✅ — matches exactly.

## Physical Analysis

B's three-fold explanation is correct:

1. **Low concentrations at high z**: c ≈ 3.3 at z=15 → ρ_s much lower than z=0 halos → longer relaxation time
2. **High virial velocities**: v_vir ~ 450–2000 km/s at these masses → velocity-dependent σ/m is suppressed (Bullet Cluster regime)
3. **Short t_universe**: At z=15, t_univ = 0.269 Gyr — need t_gc << 1 Gyr but get t_gc ~ 10⁴ Gyr

The ratio t_gc/t_univ ~ 10⁴ is robust. Even order-of-magnitude changes in the prefactor 150 don't help — would need factors of 10⁴ reduction, which is unphysical.

## Key Insight Confirmed
The velocity dependence that makes the model "cluster-safe" (σ/m drops at v > 100 km/s) is THE SAME mechanism that prevents SMBH seeding. This is not a tunable parameter — it's a structural consequence of the Yukawa potential.

## Verdict: ✅ PASSED — All results reproduced. Negative result is robust and physically well-understood.
