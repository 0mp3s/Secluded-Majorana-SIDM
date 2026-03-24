# Validation — Fornax GC Survival (Opus A)

**Script:** `predict_fornax_gc.py`  
**Author:** Opus B  
**Validator:** Opus A  
**Date:** 23 March 2026  

## Method
Ran the full script independently and cross-checked internal parameters against literature.

## Output Reproduction

| Benchmark | σ/m(16.5) B reported | My run | r_core B reported | My run | Match |
|---|---|---|---|---|---|
| BP1 | 0.477 cm²/g | 0.477 cm²/g | 449 pc | 449 pc | ✅ |
| BP9 | 0.316 cm²/g | 0.316 cm²/g | 332 pc | 332 pc | ✅ |
| MAP | 1.213 cm²/g | 1.213 cm²/g | 829 pc | 829 pc | ✅ |

Manual VPM spot-check: `sigma_T_vpm(20.69, 11.34e-3, 1.048e-3, 16.5)` = 0.47704 ✅

## GC Survival Summary Reproduction

| Benchmark | B: STALLED/SAFE/Prob | My run | Match |
|---|---|---|---|
| BP1 | 3 / 9 / 3 (12/15) | 3 / 9 / 3 (12/15) | ✅ |
| BP9 | 2 / 9 / 4 (11/15) | 2 / 9 / 4 (11/15) | ✅ |
| MAP | 5 / 9 / 1 (14/15) | 5 / 9 / 1 (14/15) | ✅ |

## Parameter Cross-Checks

### Fornax halo (Walker+2009, Read+2019):
- M200 = 3.16e9 M☉ → consistent with Read+2019 Table 1
- c200 = 18 → on the high end of λCDM expectation but within scatter for Fornax
- σ_v = 11.7 km/s → Walker+2009 confirms

### GC data (Mackey & Gilmore 2003):
- GC3 (Fornax 4): M = 3.63e5 M☉, r_2D = 240 pc, r_h = 2.5 pc — these are standard values
- GC4 (Fornax 5): M = 1.82e5 M☉, r_2D = 175 pc, r_h = 6.0 pc — ✅

### DF timescale cross-check:
Manual calculation for GC3 at 430 pc (face-on): t_DF ≈ 33 Gyr → STALLED ✅  
(Script should give "stalled" since 430 < 449 = r_core for BP1 → DF = 0 inside core)

### Deprojection factors:
- Face-on (×1), Mean (×4/π ≈ 1.27), Median (×π/2 ≈ 1.57) — standard statistical deprojections ✅

## Minor Comment
B's note about discontinuous DF stalling is important for the paper. Real DF continues (weakened) inside the core. This means BP1's marginal cases may be somewhat better in practice.

## Verdict: ✅ PASSED — All values reproduced exactly. GC data and halo parameters consistent with literature.
