# Validation Report — predict_ufd.py (Opus A)

**Date:** 23 March 2026  
**Script:** `model_validations/ufd_crater/predict_ufd.py`  
**Author of script:** Opus B  
**Validator:** Opus A

## 1. Code Review

- Methodology: Kaplinghat+2016 criterion $\rho(r_1) \times (\sigma/m) \times \sigma_v \times t_{\rm age} = 1$, binary search for $r_{\rm core}$ on NFW profile.  
- 15 galaxies: 8 classical dSphs + 6 UFDs + Crater II.  
- 3 benchmarks: BP1, BP9, MAP.  
- Uses `sigma_T_vpm()` from core/v22_raw_scan.py — correct (σ/m depends on α_s only).  
- NFW profile: Correa+2015 c(M,z=0) with h=0.674 — consistent with rest of codebase.  
- $t_{\rm age} = 10$ Gyr — conservative (see discussion on Sextans/Leo I).  
- Core criterion threshold: $N_{\rm scatter} \geq 1$ — standard.

**No bugs found.** Logic, units, and physics all correct.

## 2. Full Script Execution

Ran `predict_ufd.py` independently. All values reproduced exactly:

| Galaxy | Benchmark | σ/m (A) | σ/m (B) | N_scatter (A) | N_scatter (B) | State |
|---|---|---|---|---|---|---|
| Fornax | MAP | 1.354 | 1.354 | 6.93 | 6.93 | CORED ✅ |
| Sextans | MAP | 1.408 | 1.408 | 0.99 | 0.99 | CUSPY ✅ |
| Leo I | MAP | 1.363 | 1.363 | 0.85 | 0.85 | CUSPY ✅ |
| Crater II | MAP | 1.553 | 1.553 | 1.88 | 1.88 | CORED ✅ |
| Tucana III | MAP | 1.567 | 1.567 | — | — | ✅ |
| Fornax | BP1 | 0.476 | 0.476 | — | — | CUSPY ✅ |
| Crater II | BP1 | 0.151 | 0.151 | — | — | CUSPY ✅ |

BP1/BP9: ALL 15 galaxies CUSPY — confirmed.  
MAP: 6/8 classical CORED, 4/6 UFDs CORED — confirmed.  
Crater II MAP: r_core = 197 pc, r_core/r_half = 0.19 — confirmed.

## 3. VPM Cross-Checks (6 Manual Calls)

Independent `sigma_T_vpm()` calls to verify the transfer cross-section:

| Call | Expected | Got | Match |
|---|---|---|---|
| MAP, v=7.6 km/s (Sextans) | 1.4078 | 1.4078 | ✅ |
| MAP, v=2.3 km/s (Crater II) | 1.5529 | 1.5529 | ✅ |
| BP1, v=2.3 km/s (Crater II) | 0.1506 | 0.1506 | ✅ |
| MAP, v=11.5 km/s (Fornax) | 1.3543 | 1.3543 | ✅ |
| MAP, v=9.2 km/s (Leo I) | 1.3628 | 1.3628 | ✅ |
| MAP, v=1.6 km/s (Tuc III) | 1.5671 | 1.5671 | ✅ |

All 6/6 exact matches to 4 decimal places.

## 4. Physics Assessment

- **BP1/BP9 all cuspy:** Consistent with σ/m(v<12) < 0.5 cm²/g — too low for N_scatter > 1 within 10 Gyr.  
- **MAP 6/8 classical CORED:** σ/m ≈ 1.35–1.46 cm²/g in the 7–12 km/s range drives efficient coring.  
- **Sextans N=0.99 borderline:** Within t_age and M_200 uncertainties. Using t_age=12 Gyr → N ≈ 1.19.  
- **Leo I N=0.85:** Genuinely below threshold. Observational evidence for Leo I's core is ambiguous (Mateo+2008). Not a tension.  
- **r_core/r_half NOT universal (67–86% scatter):** Direct consequence of velocity-dependent σ/m. Testable prediction.  
- **Crater II partial:** r_core=197 pc << r_half=1066 pc. Tidal processing required. Consistent with Fattahi+2018, Fu+2019.

## 5. Verdict

**PASSED.** All numerical results independently reproduced. Physics interpretation sound. No corrections needed.
