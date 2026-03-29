# דוח בעיות בסקריפטים — ביקורת Opus מקיפה

**תאריך**: 27 מרץ 2026  
**בודק**: סוכן Opus (Claude Opus 4.6), ביקורת עצמאית מלאה

---

## סיכום

| סטטוס | כמות | סקריפטים |
|--------|-------|---------|
| ✅ נקי | 12 | sigma_radiative_stability, fifth_force_constraints, freeze_out_analysis_corrected, verify_cw_bugs, theta_topological, sigma_trapping_ode, verify_a4_cg, test_alpha_convention, dark_force_accumulation, sigma_mass_protection, born_full_amplitude, vev_alignment_stability |
| ⚠️ מיושן חלקית | 1 | dark_axion_full (CW min at θ=0) |
| ❌ Superseded (קונבנציה) | 6 | consistency_check_sidm, consistency_17bp, free_theta_scan, boltzmann_17bp, fornax_gc_check, resonance_bp1 |
| ❌ Superseded (באגים) | 1 | freeze_out_trapping (5 באגים) |
| ✅ חקרני (ללא שגיאות) | 2 | a4_dark_sector_model, dark_qcd_consistency |

---

## ממצאים לפי סקריפט

### 1. `sigma_radiative_stability.py` — Test 1 ✅ CLEAN
- בודק 4 אפשרויות לסקלר σ (A–D). רק אפשרות D (dark axion) שורדת.
- קבועים נכונים: M_Pl = 2.435e18, H_0 = 1.44e-42.
- לא תלוי בקונבנציית α (ניתוח גנרי).

### 2. `dark_axion_full.py` — Test 2 ⚠️ PARTIALLY OUTDATED
- גזירת θ_relic (tan²θ = 1/8 מ-CG) **נכונה**.
- **באג**: טוען CW minimum ב-θ = 0 — תוקן ב-Test 4 ל-θ = π/2.
- ערכי f_needed עקביים פנימית עם journal.

### 3. `fifth_force_constraints.py` — Test 3 ✅ CLEAN
- β = 5 עבור f = 0.2 M_Pl. CW mass screening נכון.
- קבועים נכונים.

### 4. `freeze_out_trapping.py` — Test 4 ❌ SUPERSEDED — 5 באגים מתועדים

| Bug | תיאור | חומרה |
|-----|--------|--------|
| 1 | מקדם CW: `[ln−0.5]` במקום `[ln−1.0]` | גבוהה — factor-of-2 בנגזרת |
| 2 | CW minimum ב-θ=0 (שגוי) — צריך θ=π/2 | גבוהה — כיוון הגלגול הפוך |
| 3 | m²_σ sign שגוי | גבוהה |
| 4 | `t_eval` float במקום array | crash |
| 5 | סיכום hard-coded לא תואם קוד | קוסמטי |

### 5. `verify_cw_bugs.py` — Test 4 Audit ✅ CLEAN
- מזהה ומאמת נומרית את כל 5 הבאגים.

### 6. `freeze_out_analysis_corrected.py` — Test 4 Fix ✅ CLEAN
- כל 5 הבאגים תוקנו. CW minimum ב-π/2. F_th/F_CW ~ 10⁻¹⁰.

### 7. `theta_topological.py` — Test 5 ✅ CLEAN
- מראה θ/(2π) אירציונלי → Z_N לא עובד.
- חיבור טטרהדרלי/A₄/S₄ נחקר נכון.

### 8. `consistency_check_sidm.py` — Test 6 ❌ SUPERSEDED
- **שגיאת קונבנציה**: מכפיל ב-(8/9) את α_CSV שכבר הוא α_s.
- VPM solver עצמו תקין. Majorana spin weights (w_even=1, w_odd=3) נכונים.
- Header SUPERSEDED נוסף.

### 9. `consistency_17bp.py` — Test 7 ❌ SUPERSEDED
- אותה שגיאת קונבנציה כמו Test 6.
- טוען CSV מ-`Secluded-Majorana-SIDM/predictions/output/sweep_17bp_results.csv`.
- Header SUPERSEDED נוסף.

### 10. `free_theta_scan.py` — Test 8 ❌ SUPERSEDED
- שגיאת קונבנציה: α_s = α × cos²θ כאשר α המקורי כבר α_s.
- מגמות איכותיות (מבנה תהודה vs θ) עדיין תקפות; ערכי σ/m מוחלטים ~11% נמוכים.
- Header SUPERSEDED נוסף.

### 11. `a4_dark_sector_model.py` — Test 9 ✅ CLEAN (חקרני)
- חוקר מנגנונים מרובים להפקת θ = arcsin(1/3) מ-A₄.
- ממצא מפתח: עם ψ = (1,1,1)/√3, ξ_s = (1,1,1), ξ_p = (1,0,0) → g_s = 3, g_p = 1 → sin²θ = 1/10, **לא** 1/9.
- **פער מוכר**: סגירת ה-10% דורשת v_p/v_s ≈ 1.061 (תיקון O(1) טבעי).
- גנרטורי S, חוקי מכפלה, ו-VEV alignment מקודדים נכון.

### 12. `boltzmann_17bp.py` — Test 10 ❌ SUPERSEDED
- שגיאת קונבנציה: אותה בעיית (8/9).
- Boltzmann solver מלא (RK4) מיושם היטב.
- Header SUPERSEDED נוסף.

### 13. `sigma_trapping_ode.py` — Test 10 (trapping) ✅ CLEAN
- ODE מצומד (θ, ω, Y) עם כוחות CW + תרמיים.
- משתמש במקדם CW **מתוקן** ([ln − 1.0]).
- מסקנה נכונה: σ לא כלוא אלא אם V_bare קיים או θ בדיד (A₄).
- חלק 5 מנסח בבהירות: θ הוא קבוע תורת חבורות, לא שדה דינמי.

### 14. `verify_a4_cg.py` — Test 10 (CG) ✅ CLEAN עם ממצא חשוב
- A₄ group מאומתת: S² = T³ = (ST)³ = 1, |A₄| = 12. ✅
- CG 3⊗3→1 ו-3⊗3⊗3→1 מאומתות אינווריאנטיות תחת כל 12 האלמנטים. ✅
- תוצאה מפתח: g_s = 3, g_p = 1, **tan²θ = 1/9 → sin²θ = 1/10 → θ = 18.43°** (לא 19.47°).
- **חשוב**: הצימוד 3⊗3⊗3→1 **לא** סימטרי תחת החלפת a↔b → חייבים לסמטרז ל-Majorana. עבור ψ = (1,1,1)/√3, סימטריזציה לא משנה תוצאה.
- יחס VEV: v_p/v_s = 3/(2√2) ≈ 1.061 סוגר פער ל-sin²θ = 1/9. רק 6.1% תיקון.

### 15. `fornax_gc_check.py` — Test 11 ❌ SUPERSEDED
- שגיאת קונבנציה: מכפיל ב-(8/9) את α_CSV.
- עמודת "without decomposition" נותנת תוצאות נכונות.
- Header SUPERSEDED נוסף.

### 16. `test_alpha_convention.py` — Test 12 ✅ SCRIPT מפתח
- **מוכיח**: Scenario A (עם 8/9): 10/17 BPs עוברים. Scenario B (ללא הפחתה): **17/17 עוברים** AND תואמים CSV בפחות מ-2%.
- פותר את שגיאת הקונבנציה ב-6 סקריפטים.

### 17. `dark_force_accumulation.py` — Test 13 ✅ CLEAN
- בודק האם שדות φ/σ מצטברים מ-DM background = ρ_Λ.
- משתמש ב-α = 5.734e-3 (MAP; **מתוקן** מ-3.4e-3 קודם).
- **פיזיקה נכונה**: φ-mediated force מתמסה אקספוננציאלית בטווח קוסמולוגי; שדה σ צריך m_σ ~ H₀ ל-DE, אבל CW נותן m_σ ~ 10⁻¹² eV → פער של 40 סדרי גודל.

### 18. `sigma_mass_protection.py` — Test 13D ✅ CLEAN (חקרני)
- מחשב CW σ mass בשני מקרים: (A) בלי φ VEV → 2-loop sunset → קטן; (B) עם φ VEV אפקטיבי → 1-loop.
- בודק 5 מנגנוני הגנה: Goldstone, dark SUSY, clockwork, extra dimension, radiative seesaw.
- מזהה נכון בעיית f טרנס-פלאנקי.

### 19. `dark_qcd_consistency.py` — Test 14 ✅ CLEAN (מעמיק)
- 6 בדיקות עקביות: תורת חבורות, יחס CP, BBN/N_eff, Λ_d מ-RG, צירוף מקרים, portal.
- תוצאות מפתח: SU(2)_d מועדף (Majorana + BBN). Λ_d ~ 10⁻³ eV דורש α_d(M_Pl) ~ 1/200.
- Λ_d ≈ √(H₀ M_Pl) ~ סקלת מסת נייטרינו — צוין כמרשים.
- Misalignment: Ω_σ = 0.69 עבור θ_i ~ 2 rad.
- **אין שגיאות**. ניתוח BBN ו-RG running נכונים.

### 20. `born_full_amplitude.py` — Test 15 ✅ CLEAN
- אמפליטודת Born מלאה (Dirac algebra) ל-χχ → χχ עם scalar + pseudoscalar vertices.
- α_p = α_s/8 (מ-tan²θ = 1/8). ✅
- **תוצאה מפתח**: תיקון pseudoscalar < 10⁻⁶ מ-scalar. VPM(α_s) מאומת.
- Interference נעלם אחרי spin averaging: Tr[γ⁵ עם < 4 gammas] = 0. ✅

### 21. `vev_alignment_stability.py` — Test 16 ✅ CLEAN
- מאשר: λ₂ > 0 → (1,1,1)/√3; λ₂ < 0 → (1,0,0). ✅
- Two-flavon: רוב מרחב הפרמטרים משמר (1,1,1)×(1,0,0). ✅
- Hessian eigenvalues נבדקים (minima אמיתיים). ✅

### 22. `resonance_bp1.py` — ❌ SUPERSEDED
- שגיאת קונבנציה: מכפיל ב-(8/9).
- m_φ(BP1) **תוקן** מ-11.34 → 10.83 MeV.
- VPM solver ו-resonance scan נכונים; רק ערך α שהוזן היה שגוי.
- Header SUPERSEDED נוסף.

---

## בעיות שנותרו פתוחות

1. **`dark_axion_full.py`** — הפרשנות "CW min at θ=0" בטקסט שגויה (צריך θ=π/2). הסקריפט עצמו מחשב V_CW נכון; רק הטקסט מטעה.
2. **sin²θ = 1/10 vs 1/9** — פער תיאורטי (6.1%) שנסגר ע"י יחס VEV אבל לא מגיע "חינם" מ-A₄.
3. **VPM copies** — 7+ עותקים של אותו VPM solver. לא באג, אבל סיכון לדריפט בין גרסאות.
