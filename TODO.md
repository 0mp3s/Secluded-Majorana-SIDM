# TODO — Next Steps

משימות פתוחות לאחר השלמת ה-consistency checks (predictions/).

| # | משימה | סוג | Impact | מאמץ | עלול לשבור מודל | זמן עבודה | זמן הרצה | סטטוס |
|---|--------|------|--------|------|----------------|-----------|-----------|--------|
| 1 | **SPARC + baryons rotation curve fit** | ואלידציה | ⭐⭐⭐⭐⭐ | בינוני | לא | — | — | ✅ הושלם |
| 1b | **Sensitivity analysis** (25,200 pts) | ואלידציה | ⭐⭐⭐⭐⭐ | בינוני | לא | — | — | ✅ הושלם |
| 2 | **Fornax core profile** (Jeans solver) | תחזית | ⭐⭐⭐⭐ | בינוני | לא — תחזית עצמאית, לא משנה פרמטרים | ~2 שעות | ~5 דק' | 🔴 פתוח |
| 3 | **Fermi-LAT dSph bounds** (φφ → dark, secondary γ) | constraint | ⭐⭐⭐ | בינוני | לא — secluded → zero SM signal | ~1 שעה | ~1 שניה | ✅ הושלם |
| 4 | **MCMC שרשרת ארוכה יותר** (3000+ צעדים, 50τ) | סטטיסטי | ⭐⭐⭐ | ~50 דק' | לא — אותו מודל, סטטיסטיקה טובה יותר | ~30 דק' | ~50 דק' | 🔴 פתוח |
| 5 | **דיון yₛ ≫ yₚ hierarchy** בפרי-פרינט | טקסט | ⭐⭐ | נמוך | לא | — | — | ✅ הושלם |
| 6 | **ξ evolution** — power law vs log (Farina+16 eq.9) | טקסט | ⭐ | נמוך | לא | — | — | ✅ הושלם |
| 7 | **small-μ₃ inconsistency** — Ω_φ overclosed note | טקסט | ⭐ | נמוך | לא | — | — | ✅ הושלם |

## פירוט

### 1. SPARC + Baryons Fit  ✅
- Fit עקומות סיבוב מלאות: V²_tot = Υ_* V²_bar + V²_SIDM
- 7 גלקסיות מ-SPARC (~80 data points), c(M₂₀₀) Dutton+2014 + adiabatic contraction
- **תוצאה:** 3/4 ננסיות נותנות Υ_* פיזיקלי (0.2–0.8)
- Sensitivity analysis: 25,200 חישובים, סריקת 4 פרמטרים אסטרופיזיקליים
- **פרמטר דומיננטי:** c_factor (concentration scatter); ננסיות לא צריכות AC
- ספירליות דורשות AC מתוחכם יותר (Gnedin+2004) — לא כשל של SIDM

### 2. Fornax Core Profile
- לפתור Jeans equation עם פוטנציאל SIDM
- לנבא r_core [pc] ופרופיל ρ(r) מדויק
- בדיק מול GAIA DR3 + stellar kinematics
- תחזית אמיתית (blind prediction)

### 3. Fermi-LAT dSph Bounds
- χχ → φφ (secluded annihilation)
- φ → dark sector products (no SM tree-level)
- Secondary γ from bremsstrahlung / inverse Compton
- בדיקת upper limits מול Fermi-LAT dSph stacking

### 4. MCMC שרשרת ארוכה יותר
- emcee מזהיר ש-chain length < 50τ (כרגע 2000 צעדים, צריך ≥3000)
- הרצה מחדש של `stats_mcmc/opusB_run_mcmc.py` עם `N_STEPS = 3000`
- ~50 דקות ריצה
- נדרש לגרסה סופית של הפרי-פרינט

### 5. דיון yₛ ≫ yₚ Hierarchy
- Opus B העיר: הסבירו למה yₛ ≫ yₚ (יחס 13–212) טבעי
- אפשרויות: radiative origin ל-yₚ, approximate CP symmetry, loop-suppression
- להוסיף פסקה בדיון (Discussion section) של הפרי-פרינט

### 6. ξ Evolution — Power Law vs Logarithmic
- Opus B ציין: ξ = T_φ/T_γ מתפתח כ-power law (gₛ ratios), לא log-correction
- ההערה נכונה אבל ההשפעה על התוצאות קטנה (qualitative argument)
- להוסיף footnote או הבהרה בפרי-פרינט

### 7. Small-μ₃ Inconsistency Note
- Opus B ציין: באזור μ₃/m_φ < threshold, φ overclosed → inconsistency
- אין השפעה על viable region (שם cannibal עובד)
- להוסיף משפט הבהרה שאזור ה-overclosure הוא excluded by construction

## הושלם

- [x] Gravothermal collapse timescales (predictions/gravothermal/)
- [x] SPARC rotation curve diversity — DM only (predictions/rotation_curves/)
- [x] Cluster merger offsets (predictions/cluster_offsets/)
- [x] ΔN_eff from light mediator (predictions/delta_neff/)
- [x] MCMC posterior sampling (stats_mcmc/)
- [x] χ² fit to 13 astrophysical systems (observations/)
- [x] Relic density — 17 viable BPs (relic_density/)
- [x] Born validation, literature cross-checks (cross_checks/)
- [x] Sommerfeld enhancement (cross_checks/)
- [x] Velocity-averaged cross sections (cross_checks/)
- [x] BBN / mediator cosmology (cosmology/)
- [x] SPARC + baryons rotation curve fit (predictions/rotation_curves/)
- [x] Sensitivity analysis — 25,200 point parameter scan (predictions/rotation_curves/)
- [x] SIDM dwarf upper bound 50→10 cm²/g (Opus B review, condition2 rescan)
- [x] §7 Phenomenological Predictions — all 6 subsections (§7.1–§7.6) integrated into preprint
- [x] דיון yₛ ≫ yₚ hierarchy — §8.4 naturalness paragraph added
- [x] ξ evolution — power law clarification with Farina+16 reference in §5.2
- [x] small-μ₃ overclosure — cannibal annihilation note added in §5.3
- [x] Fermi-LAT dSph bounds — trivially satisfied (secluded → zero SM signal, cosmology/fermi_lat_dsph.py)
