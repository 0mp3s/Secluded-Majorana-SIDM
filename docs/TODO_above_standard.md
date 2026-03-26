# TODO — Above Standard: מחקרי · מדעי · אקדמי · פיזיקלי

> **מטרה:** המודל הזה לא רק fit מתמטי לנתונים — הוא טוען משהו על המציאות.
> כדי שזה ילקח ברצינות, אנחנו חייבים להיות **מעל הסטנדרט** בכל ממד.
> המסמך הזה מרכז את כל מה שצריך לעשות כדי להגיע לשם.

---

## A. פיזיקלי — "האם המודל אומר משהו על היקום?"

זו השאלה הכי חשובה. מודל SIDM עם Yukawa mediator שמתאים ל-χ² לא מספיק —
צריך להראות ש:
1. המודל **מנבא** תופעות שלא הוכנסו ל-fit
2. ניבויים אלה **נבדקים** מול נתונים עצמאיים
3. יש **falsifiable predictions** — מה ישלול את המודל

### A1. ✅ כבר קיים

| פריט | סטטוס | מיקום |
|------|--------|-------|
| Relic + SIDM simultaneous (122 BPs) | ✅ | relic_density/ |
| 12 predictions עצמאיות | ✅ | predictions/, model_validations/ |
| Fornax GC survival | ✅ | model_validations/fornax_gc/ |
| SMBH seeds — falsifiable by JWST | ✅ | model_validations/smbh_seeds/ |
| UFD scaling (Crater II, Segue 1) | ✅ | model_validations/ufd_crater/ |
| Majorana vs Dirac fingerprint | ✅ | predictions/majorana_vs_dirac/ |
| σ_SI = 0 prediction (secluded) | ✅ | Built into model |
| ΔN_eff ≈ 0 | ✅ | predictions/delta_neff/ |
| Adiabatic contraction (Blumenthal+86) | ✅ | predictions/rotation_curves/fit_sparc_baryons.py |
| c(M) relation (Dutton & Macciò 2014) | ✅ | predictions/rotation_curves/fit_sparc_baryons.py |

### A2. 📋 לעשות

- [ ] **Gnedin+04 adiabatic contraction** (upgrade מ-Blumenthal)
  - רקע: Blumenthal+86 הוא instantaneous — מניח שהבריונים "נופלים" לפני שה-DM מגיב.
    Gnedin+04 מוסיף modified adiabatic invariant שמתאים טוב יותר ל-N-body.
  - מטרה: האם spiral rotation curves *עובדים* עם AC מתקדם יותר?
  - מה: להוסיף `gnedin_ac()` ב-`fit_sparc_baryons.py`, להריץ על 7 גלקסיות
  - effort: ~2 שעות | impact: גבוה — הופך failure ל-result
  - **priority: HIGH**

- [ ] **Gravothermal fluid equations**
  - רקע: כרגע $t_{\rm grav}/t_{\rm age}$ estimate בלבד (Essig+19 eq. 2)
  - מטרה: core evolution curve $r_c(t)$ — **כמה שנים** עד קריסה, לא רק "כן/לא"
  - מה: סקריפט חדש `predictions/gravothermal/fluid_solver.py`, Balberg+02 4 ODEs
  - effort: ~1 יום | impact: בינוני-גבוה
  - **priority: MEDIUM** — nice to have, לא critical

- [ ] **Diversity problem test מפורש**
  - רקע: Oman+15 הראו שב-ΛCDM rotation curves הם "too similar". SIDM אמור לתת diversity.
  - מטרה: scatter plot of V(2 kpc) vs V_max עבור BPs — האם ה-diversity תואם נתונים?
  - מה: כבר יש predict_core_sizes.py — צריך פלוט נוסף + השוואה ל-Oman+15 Figure 4
  - effort: ~1 שעה | impact: גבוה — diversity problem הוא ה-*motivation* של SIDM
  - **priority: HIGH**

---

## B. מדעי — "האם ניתן לשחזר ולאמת?"

### B1. ✅ כבר קיים

| פריט | סטטוס |
|------|--------|
| 54-stage pipeline | ✅ |
| global_config.json — SSOT | ✅ |
| 8 cross-checks עם סקריפטים נפרדים | ✅ |
| Sensitivity analysis (25,200 pts) | ✅ |
| MCMC with emcee | ✅ |
| Reproducible seeds (42) | ✅ |

### B2. 📋 לעשות

- [ ] **N-body calibration paragraph**
  - מה: פסקה אחת שמצטטת Rocha+13, Robertson+17 — isothermal core approach תואם N-body ב-$0.5 < \sigma/m < 10$
  - מטרה: מנטרל referee question מראש. לא דורש קוד
  - effort: 10 דקות | impact: גבוה (defensive)
  - **priority: HIGH** — להכניס למאמר

- [ ] **Error propagation chain מפורש**
  - מה: מסמך שמראה כיצד שגיאות ב-VPM (§3.3) מתרחשות דרך הפייפליין: VPM → σ/m → χ² → MCMC posteriors
  - מטרה: referee שישאל "how sensitive are your posteriors to the 11% error at 1000 km/s?"
  - effort: ~3 שעות | impact: בינוני-גבוה
  - **priority: MEDIUM**

- [ ] **Convergence test ל-MCMC**
  - מה: Gelman-Rubin $\hat{R}$ diagnostic (כרגע יש רק autocorrelation time)
  - מטרה: סטנדרט אקדמי — emcee paper ממליצים גם Gelman-Rubin
  - effort: ~30 דקות (כבר יש chains) | impact: בינוני
  - **priority: MEDIUM** — קל לעשות, טוב להיות שם

- [ ] **Refactor VPM solver to shared module**
  - רקע: VPM solver (phase shifts, σ_T) מועתק ב-copy-paste ל-6 סקריפטים ב-vpm_scan/.
    כשתוקנה נוסחת ℓ_max ב-v22_raw_scan.py, 5 מתוך 6 העותקים לא עודכנו — גרם לבאג קריטי.
  - מה: להוציא `vpm_phase_shift()`, `sigma_T_vpm()`, ונוסחת ℓ_max לתוך `core/vpm_solver.py`.
    כל סקריפטי vpm_scan/ ייבאו משם במקום לשכפל.
  - מטרה: DRY — למנוע באגים עתידיים כשמשנים את הסולבר
  - effort: ~2 שעות | impact: גבוה (מניעת באגים)
  - **priority: MEDIUM** — אחרי סיום pipeline run

---

## C. אקדמי — "האם זה ברמת פרסום?"

### C1. ✅ כבר קיים

| פריט | סטטוס |
|------|--------|
| LaTeX preprint (main.tex) | ✅ |
| research_journal_v10 (1645 שורות) | ✅ |
| methodology.md מקיף | ✅ |
| CITATION.cff | ✅ |
| Peer review by 3 AI systems | ✅ |

### C2. 📋 לעשות

- [ ] **Limitations section מפורש במאמר**
  - מה: §8 "Limitations & Future Work" — 5 נקודות:
    1. No N-body (calibrated against Rocha+13)
    2. Isothermal core approximation (valid for $N_{\rm scatter} \sim 1$–100)
    3. Spiral rotation curves need higher-fidelity NFW contraction
    4. mχ range limited to 10–100 GeV (TeV deferred)
    5. φ cosmology (cannibal depletion) parametrized, not simulated
  - מטרה: referees *respect* papers that acknowledge limitations. מונע rejection
  - effort: ~1 שעה | impact: **קריטי**
  - **priority: CRITICAL**

- [ ] **Human physicist review**
  - מה: לפני submission — לשלוח ל-physicist (arXiv endorser?) לביקורת
  - מטרה: AI peer review ≠ human peer review. מהנדס תוכנה + AI ← צריך עוד שכבת אישור
  - effort: תלוי בקשר | impact: **קריטי**
  - **priority: CRITICAL** — ללא זה, קשה לפרסם

- [ ] **arXiv endorsement**
  - מה: צריך endorser ב-hep-ph או astro-ph.CO
  - מטרה: ללא endorsement אי-אפשר להעלות. זה prerequisite
  - effort: networking | impact: blocker
  - **priority: CRITICAL** — להתחיל *עכשיו*, לפני שהמאמר מוכן

- [ ] **Comparison table מול מודלי SIDM אחרים**
  - מה: טבלה שמשווה secluded Majorana vs. light mediator (Tulin+13), atomic DM (Kaplan+09), dark photon (Ackerman+09)
  - מטרה: ממקם את המודל שלנו בהקשר הרחב. כל מאמר SIDM רציני עושה את זה בסעיף 1 או 2
  - effort: ~2 שעות (ספרות) | impact: גבוה
  - **priority: HIGH**

---

## D. מחקרי — "האם התהליך מוצק?"

### D1. ✅ כבר קיים

| פריט | סטטוס |
|------|--------|
| Iterative process documented | ✅ (journal) |
| Bug discovery → fix → rerun cycle | ✅ |
| Peer review → FATAL → resolution | ✅ |
| Sensitivity analysis | ✅ |

### D2. 📋 לעשות

- [ ] **Blind analysis element**
  - מה: להוסיף סעיף שמסביר שה-SIDM cuts נלקחו *מהספרות* לפני שראינו תוצאות, לא tuned
  - מטרה: blinding מראה יושרה מדעית. גם אם לא full blinding, ההפרדה בין cuts ל-results חשובה
  - effort: פסקה אחת | impact: בינוני
  - **priority: MEDIUM**

- [ ] **Parameter space coverage analysis**
  - מה: Show that 194K points cover the viable space uniformly — no cherry-picking
  - מטרה: מישהו עלול לשאול "maybe you missed the region where the model fails"
  - effort: ~1 שעה (density plot) | impact: בינוני
  - **priority: LOW**

---

## סיכום לפי Priority

### 🔴 CRITICAL — ללא אלה אין פרסום

| # | פריט | סוג | סעיף |
|---|------|-----|------|
| 1 | Limitations section במאמר | טקסט | C2 |
| 2 | Human physicist review | חיצוני | C2 |
| 3 | arXiv endorsement | חיצוני | C2 |

### 🟠 HIGH — מעלים את העבודה מעל הסטנדרט

| # | פריט | סוג | סעיף |
|---|------|-----|------|
| 4 | Gnedin+04 AC upgrade | קוד | A2 |
| 5 | Diversity problem plot (Oman+15) | קוד + plot | A2 |
| 6 | N-body calibration paragraph | טקסט | B2 |
| 7 | Comparison table vs other SIDM models | טקסט/ספרות | C2 |

### 🟡 MEDIUM — מחזק, לא חובה

| # | פריט | סוג | סעיף |
|---|------|-----|------|
| 8 | Error propagation chain | analysis | B2 |
| 9 | Gelman-Rubin $\hat{R}$ | קוד (קצר) | B2 |
| 10 | Gravothermal fluid solver | קוד | A2 |
| 11 | Blind analysis paragraph | טקסט | D2 |

### 🟢 LOW — bonus

| # | פריט | סוג | סעיף |
|---|------|-----|------|
| 12 | Parameter space coverage plot | plot | D2 |

---

## הערה אישית

> המוטיבציה של הפרויקט הזה לא אקדמית — היא שאלה אמיתית על יקום.
> הרקע שלי (הנדסת תעשייה וניהול, תוכנה) ← **יתרון** ב-reproducibility, pipeline, testing.
> אבל **חובה** להוסיף: physicist review, proper limitations, context vs. other models.
> עם הפריטים הנ"ל, העבודה תהיה ברמה שלא ניתן לפסול בגלל רקע — רק לשפוט לגופו של עניין.
