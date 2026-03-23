# §7.1 CP-Separation Table — Validation by Opus B

**Date:** 23 March 2026  
**Script:** `model_validations/cp_separation/cp_separation_table.py`  
**Output:** `model_validations/cp_separation/output/cp_separation_table.csv`

---

## 1. α_p = (α_s × α_p)_relic / α_s — ✅ CORRECT

Verified formula: `alpha_p = RELIC_PRODUCT / alpha_s` with `RELIC_PRODUCT = 1.387474e-7`.

| α_s | α_p (script) | α_p (manual) | Match? |
|---|---|---|---|
| 1.3384e-3 | 1.036666e-4 | 1.387474e-7 / 1.3384e-3 = 1.03667e-4 | ✅ |
| 2.6939e-3 | 5.150429e-5 | 1.387474e-7 / 2.6939e-3 = 5.15043e-5 | ✅ |
| 5.4222e-3 | 2.558876e-5 | 1.387474e-7 / 5.4222e-3 = 2.55888e-5 | ✅ |

## 2. λ = α_s m_χ / m_φ (NO factor 2) — ✅ CORRECT

Convention matches our model: λ = α m_χ / m_φ. **No factor-of-2 ambiguity.**

| α_s | λ (script) | λ (manual: α_s × 20.69 / 0.01134) | Match? |
|---|---|---|---|
| 1.3384e-3 | 2.4419 | 2.442 | ✅ |
| 2.6939e-3 | 4.9151 | 4.915 | ✅ |
| 5.4222e-3 | 9.8929 | 9.893 | ✅ |

## 3. σ/m spot-checks via independent VPM calls — ✅ ALL MATCH

Ran `sigma_T_vpm(20.69, 11.34e-3, alpha_s, v)` independently for 3 points:

| α_s | σ/m(30) script | σ/m(30) independent | σ/m(100) script | σ/m(100) independent | σ/m(1000) script | σ/m(1000) independent |
|---|---|---|---|---|---|---|
| 1.3384e-3 | 0.685426 | 0.685426 | 0.589914 | 0.589914 | 0.113488 | 0.113488 |
| 2.6939e-3 | 1.370875 | 1.370875 | 1.173849 | 1.173849 | 0.374306 | 0.374306 |
| 5.4222e-3 | 2.363605 | 2.363605 | 2.044984 | 2.044984 | 0.978176 | 0.978176 |

**All 9 spot-checks match to 6 decimal places.**

## 4. A's Critical λ Correction — ✅ CONFIRMED

| Point | α_s | λ | λ < π? |
|---|---|---|---|
| 1 | 1.3384e-3 | 2.44 | ✅ Yes |
| 2 | 1.5039e-3 | 2.74 | ✅ Yes |
| 3 | 1.6898e-3 | 3.08 | ✅ Yes (barely) |
| 4 | 1.8988e-3 | 3.46 | ❌ No (above π=3.14) |
| ... | ... | ... | ❌ |
| 13 | 5.4222e-3 | 9.89 | ❌ No |

**Result:** Only 3 of 13 points have λ < π. The remaining 10 are in the classical regime.  
This confirms A's correction. Previous paper text claiming "all below first resonance" was wrong.

## 5. Dynamic Range — ✅ CORRECT

- α_s/α_p: 12.9 — 211.9 → ratio 16.4× → **1.22 decades** ✓
- This matches the condition2 viable band width (0.61 decades for α_s, doubled by 1/α_s dependence of α_p)

## Summary

| Check | Status |
|---|---|
| α_p formula | ✅ |
| λ convention | ✅ |
| σ/m values (3 × 3 matrix) | ✅ |
| λ critical correction | ✅ |
| Dynamic range | ✅ |

**VERDICT: Script is CORRECT. No bugs found.**
