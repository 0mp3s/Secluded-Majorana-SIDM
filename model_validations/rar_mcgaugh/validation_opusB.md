# §7.3 RAR Analysis — Validation by Opus B

**Date:** 23 March 2026  
**Script:** `model_validations/rar_mcgaugh/rar_analysis.py`  
**Output:** `model_validations/rar_mcgaugh/output/rar_comparison.png` + `rar_data.csv`

---

## 1. g = V²/r Conversion — ✅ CORRECT

```python
def v_to_g(v_km_s, r_kpc):
    v_m_s = v_km_s * 1e3        # km/s → m/s
    r_m = r_kpc * KPC_M         # kpc → m  (KPC_M = 3.086e19)
    return v_m_s**2 / r_m       # centripetal acceleration [m/s²]
```

Manual check: V = 100 km/s, r = 5 kpc  
→ g = (10⁵)² / (5 × 3.086e19) = 10¹⁰ / 1.543e20 = **6.48e-11 m/s²**  
This is 0.54 × g† — reasonable for spiral galaxy at ~5 kpc. ✅

## 2. McGaugh+2016 Formula — ✅ CORRECT

```python
g_obs = g_bar / (1 - exp(-sqrt(g_bar / g†)))
```

Matches published formula: $g_{\rm obs} = \frac{g_{\rm bar}}{1 - e^{-\sqrt{g_{\rm bar}/g_\dagger}}}$

Check: g_bar = 6.48e-11 → x = √(0.54) = 0.735 → denominator = 1 - e^{-0.735} = 0.520  
→ g_obs = 6.48e-11 / 0.520 = **1.25e-10 m/s²** → g_obs/g_bar = 1.92. ✅

## 3. Υ_* Handling — ⚠️ CORRECT CODE, DATA LIMITATION

**Formula:** `g_bar = (V_bar × √Υ_*)² / r` = Υ_* × V_bar² / r ✅  
**Physics:** V_tot² = Υ_* × V_bar² + V_DM² → g_obs = g_bar + g_DM ✅

**Problem:** `fit_galaxy_ac_sidm` uses `minimize_scalar(bounds=(0.01, 3.0))`.  
All 4 dwarfs hit Υ_* = 0.01 (lower bound).

**Root cause:** `sparc_rotation_data.csv` has single `V_bar` column (pre-combined from 3.6μm photometry assuming some default Υ_*). For gas-dominated dwarfs, V_bar ≈ V_gas, and the fit pushes Υ_* → 0 because DM dominates everywhere.

**Fix needed:** Load V_gas, V_disk, V_bulge separately from SPARC. Then:
```
g_bar = (V_gas² + Υ_disk × V_disk² + Υ_bulge × V_bulge²) / r
```
Only Υ_disk and Υ_bulge are free; gas mass is fixed (21cm measurements).

This is a **data preprocessing issue**, not a physics bug. The RAR physics pipeline is correct.

## 4. SIDM Rotation Curve Fitting — ✅ CORRECT

Uses `fit_galaxy_ac_sidm` from `fit_sparc_baryons.py`:
- Iterative adiabatic contraction (Blumenthal+1986)
- SIDM core at r_1 where ρ(r_1) × (σ/m) × v × t_age = 1
- Isothermal interior, AC'd-NFW exterior
- χ² minimization for Υ_*

This is the same validated pipeline from the rotation curve predictions. ✅

## 5. Scatter Computation — ✅ CORRECT

Residuals in log-space relative to McGaugh curve:
```python
resid = log10(g_measured / g_McGaugh(g_bar))
scatter = std(resid)  # in dex
```

This matches the standard RAR scatter analysis (Li+2018, Lelli+2017). ✅

## Summary

| Check | Status | Notes |
|---|---|---|
| g = V²/r units | ✅ | km/s, kpc → m/s² |
| McGaugh formula | ✅ | Matches published |
| Υ_* physics | ✅ | V_bar_scaled = √Υ_* × V_bar |
| Υ_* bounds | ⚠️ | 0.01 lower bound; dwarfs hit it (data issue) |
| SIDM fitting | ✅ | Same validated pipeline |
| Scatter metric | ✅ | Standard dex scatter |
| V_gas/V_disk separation | ❌ | Missing — needs data update |

**VERDICT: Physics pipeline is CORRECT. Dwarf results unreliable due to V_bar pre-combination. Spiral results (Υ_* not at bounds) are trustworthy.**

### Proposed Fix for Dwarfs

1. Download SPARC individual component files (V_gas, V_disk, V_bulge)
2. Modify `load_rotation_data` to return all 3 components
3. In RAR script: `g_bar = g_gas + Υ_disk × g_disk + Υ_bulge × g_bulge`
4. Fit only Υ_disk (fix Υ_bulge = Υ_disk for simplicity)
5. Physical bound: Υ_disk ∈ [0.2, 0.8] (Meidt+2014, Schombert+2019)
