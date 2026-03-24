# Validation: CP-Separation MAP Scan — Opus B

**Script:** `cp_separation_MAP.py`  
**Output:** `output/cp_separation_MAP.csv` (500 rows)

## Spot-Checks (independent `sigma_T_vpm()` calls)

| α_s | Expected λ | Got λ | σ/m(30) script | σ/m(30) check | Match |
|---|---|---|---|---|---|
| 5.000e-4 | 3.272 | 3.2722 | 0.120739 | 0.120739 | ✅ |
| 5.734e-3 (MAP) | 48.6 | 48.5944 | 1.714390 | 1.714390 | ✅ |
| 5.000e-3 | 32.72 | 32.7220 | 0.740208 | 0.740208 | ✅ |

## Relic constraint
α_s × α_p = 1.387474e-7 verified for first row: 5.000e-4 × 2.774948e-4 = 1.387474e-7 ✅

## Key numbers
- Total viable: 500/500 ✅
- α_s/α_p range: 1.8 — 11,532 ✅
- Dynamic range: 3.81 decades ✅
- λ range: 3.27 — 261.8 ✅
- All λ > π: YES (3.27 > 3.14) ✅
- MAP in band: YES ✅

## σ/m(1000) range — CORRECTION
A reported: "σ/m(1000) stays safely below 0.11 cm²/g everywhere"  
**Actual max: 0.1935 cm²/g** (at α_s ≈ 0.0123, λ ≈ 81, resonance peak)  
**Still well below 1.0 cm²/g threshold** — all points ARE viable.  
A's table listing max = 0.106 refers to the LAST ROW only, not the global max.  
**Minor text correction needed: 0.11 → 0.19**

## Verdict: PASSED ✅ (with one minor numerical correction)
