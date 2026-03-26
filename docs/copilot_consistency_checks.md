# בדיקות עקביות פיזיקלית

נוצר: 2026-03-26 | עדכון: 2026-03-26 | מקור: ניתוח Copilot על בסיס קוד קיים + ידע פיזיקלי כללי

---

## 1. Vacuum Stability — ⚠️ OPEN (אנליטי, לא צריך קוד)

**שאלה:** הפוטנציאל $V(\phi) = \frac{1}{2}m_\phi^2\phi^2 + \mu_3\phi^3 + \lambda_\phi\phi^4$ — יציב? bounded from below?

**למה חשוב:** הכנסנו $\mu_3\phi^3$ לפריפרינט (§4.3) בלי לדון ביציבות. term קוביק בפוטנציאל סקלרי יכול לגרום ל-$V(\phi) \to -\infty$ — שבירת ואקום ספונטנית לא רצויה.

**תוצאה אנליטית (מדיון Copilot-Opus 2026-03-26):** bounded from below דורש $\lambda_\phi > 0$ ושאין מינימום עמוק יותר מ-$\phi=0$. בשילוב עם cannibal threshold $\mu_3/m_\phi \gtrsim 1.7$:
$$\lambda_\phi > \frac{3}{16}(\mu_3/m_\phi)^2 \gtrsim \frac{3}{16}(1.7)^2 \approx 0.54$$
זו הנחה סמויה שצריכה להיכנס למאמר כ-explicit assumption.

**מה קיים:** כלום. אין קוד שבודק את זה.

**מה צריך:** פסקה במאמר שמציינת $\lambda_\phi \gtrsim 0.54$ כתנאי יציבות. לא צריך סקריפט — זה חישוב אנליטי.

---

## 2. Unitarity Bounds — ✅ PASS (בוצע 2026-03-26)

**שאלה:** האם $|a_\ell| \leq 1$ לכל partial wave, לכל ה-BPs, בכל ה-velocities?

**תוצאה:** 122 BPs × 9 velocities (12–4700 km/s) = 1098 checks. **0 failures.** Worst ratio $\tilde\sigma/\tilde\sigma_{\max} = 0.043$ (4.3% of unitarity bound).

**סקריפט:** `_tmp_checks/check2_unitarity.py`

---

## 3. Perturbativity — ✅ PASS (בוצע 2026-03-26)

**שאלה:** כמה מה-BPs בתחום הפרטורבטיבי?

**תוצאה:** כל 122 BPs perturbative. $\alpha_s^{\max} = 4.81 \times 10^{-3}$, $y_s^{\max} = 0.246$, loop parameter $\alpha_s/(4\pi) \leq 0.038\%$.

**סקריפט:** `_tmp_checks/check1_perturbativity.py`

---

## 4. Warm Mediator — Free-Streaming Length

**שאלה:** $\phi$ שנשאר כ-relic לפני cannibal depletion — האם משפיע על free-streaming length ועל matter power spectrum?

**למה חשוב:** `cosmology/mediator_cosmology.py` בודק $\Delta N_{\text{eff}}$ ו-BBN, אבל לא free-streaming. מדיאטור ב-MeV עם מהירויות תרמיות יכול לפעול כ-warm dark radiation ולמחוק מבנים קטנים (כמו WDM).

**מה קיים:** `cosmology/mediator_cosmology.py` — $\Delta N_{\text{eff}}$ בלבד.

**מה צריך:** חישוב $\lambda_{\text{fs}}$ של $\phi$ לפני depletion, והשוואה ל-Lyman-$\alpha$ bound ($\lambda_{\text{fs}} \lesssim 0.1$ Mpc).

---

## 5. Sommerfeld Enhancement/Saturation

**שאלה:** ב-$m_\phi/m_\chi \sim 10^{-3.5}$ וב-cluster velocities (~1000 km/s), האם אנחנו ב-Born regime או שיש הגברת Sommerfeld לא טריוויאלית?

**למה חשוב:** VPM כבר פותר את משוואת שרדינגר מלאה, אז Sommerfeld "כלול" בתשובה. אבל השאלה היא האם ה-Born approximation (שמניח $\lambda = \alpha m_\chi/m_\phi \ll 1$) לגיטימית. עבור BP1: $\lambda \approx 3.4$ — לא קטן.

**מה קיים:** `cross_checks/sommerfeld.py` — בודק enhancement factor. `vpm_scan/vpm_born_ratio.py` — בודק VPM/Born ratio.

**מה צריך:** לוודא ש-`vpm_born_ratio.py` מכסה את כל ה-BPs ולא רק BP1/MAP. אם VPM/Born $\gg 1$, זה אומר שה-VPM תופס פיזיקה שחישוב Born לא היה תופס — וזה בסדר, כי אנחנו משתמשים ב-VPM.

---

---

## בדיקות נוספות שבוצעו (2026-03-26)

### 6. Kaplinghat+2016 Cross-Check — ✅ PASS

14 מערכות (8 dSphs + 6 clusters) מ-KTY16 Table I. כל 5 named BPs: $\chi^2/N < 2$. MAP = 1.02 (excellent fit).

**סקריפט:** `_tmp_checks/check3_kaplinghat.py`

### 7. Fornax GC Timing (Read+2019) — ⚠️ CRITICAL

$\sigma/m < 1.5$ cm²/g at $v \sim 12$ km/s. Named BPs pass (MAP margin +0.18). **64/122 relic-viable points FAIL (52% exclusion).**

Fornas GC הוא mass-space constraint: $m_\phi$ קטן → $\lambda$ גדול → $\sigma/m$ גבוה ב-low-$v$. UV completion (כמו A₄) לא יכול לעזור.

**סקריפט:** `_tmp_checks/check4_fornax_gc.py`

### 8. Oman+2015 Diversity — ✅ PASS

MC scatter ratio = 0.75 (75% of observed diversity). Remaining scatter from baryonic feedback.

**סקריפט:** `_tmp_checks/check5_oman_diversity.py`

### 9. A₄ Compatibility (מדיון Copilot-Opus) — ✅ All 122 Compatible

Factor-of-8 identity: $P_{\text{relic}} = \alpha^2/8$ (verified to 1%). A₄ UV completion ($\alpha_p/\alpha_s = 1/8$) gives $\alpha_s = \alpha_{\text{CSV}}$ — SIDM degenerate. All 122 points match A₄ within 9%.

**סקריפט:** `_tmp_checks/check_a4_compat.py`

---

## סיכום עדיפויות (מעודכן)

| # | בדיקה | סטטוס | חשיבות למאמר |
|---|--------|------|-------------|
| 3 | Perturbativity | ✅ PASS | — |
| 2 | Unitarity | ✅ PASS | — |
| 6 | KTY16 Cross-Check | ✅ PASS | — |
| 7 | Fornax GC | ⚠️ 52% exclusion | **גבוהה** — צריך למאמר |
| 8 | Oman diversity | ✅ PASS (75%) | — |
| 9 | A₄ compatibility | ✅ degenerate | — |
| 1 | Vacuum stability | ⚠️ OPEN ($\lambda_\phi \gtrsim 0.54$) | **גבוהה** — הנחה סמויה |
| 5 | Sommerfeld/Born | 🔲 לא נבדק | בינונית |
| 4 | Warm mediator FSL | 🔲 לא נבדק | נמוכה |
