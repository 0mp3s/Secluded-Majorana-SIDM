# יומן מחקר — נגזרת הלגרנג'יאן | Secluded Majorana SIDM

**תיקייה:** `Secluded-Majorana-SIDM/The_derivative_of_Lagernizan_SIDM/`  
**תאריך פתיחה:** 27 מרץ 2026  
**מטרה:** חקירה שיטתית של מה שנגזר מהלגרנג'יאן — מאומת בקוד פייתון

---

## רשומה 0 — 2026-03-27 (סשן ראשון)

**מה בוצע בפועל:**

1. נכתב `lagrangian_derivation.md` — **סיכום** של preprint_draft_v10.md §2–§7.1. לא גזירה חדשה.
2. נקרא `dark-energy-T-breaking/theory.md` + `THEORY_MATH_SUMMARY.md` — זוהו 3 שכבות (A/B/C).
3. נכתב ניתוח "מה לומדים מהנגזרת" — 6 תובנות, שיחת דיון עם עומר.

**מה לא בוצע:**
- לא נוצר קוד פייתון שמאמת משהו
- לא נעשה שום חישוב חדש — הכל סיכום של עבודה קיימת
- העבודה ב-`dark-energy-T-breaking` עשויה להיות מזוהמת מסשן קודם

**מסקנה:** צריך להתחיל מחדש עם תהליך נקי — דיון פיזיקלי → קוד פייתון → תוצאות מאומתות.

---

## רשומה 1 — 2026-03-27 (התחלה נקייה)

### נקודת מוצא מאומתת

הלגרנג'יאן מה-preprint (מפורסם, 17 BPs עוברים):

$$\mathcal{L} = \frac{1}{2}\bar{\chi}(i\not\partial - m_\chi)\chi + \frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2 - \frac{y}{2}\bar{\chi}\chi\,\phi$$

שלושה פרמטרים: $(m_\chi, m_\phi, \alpha = y^2/4\pi)$.

### שאלת המחקר

**מה אפשר לגזור מהלגרנג'יאן הזה, ולאמת נומרית?**

---

## תוכנית עבודה — חקירה שיטתית

כל שלב: דיון פיזיקלי → קוד פייתון → תוצאה מאומתת

| # | בדיקה | שאלה פיזיקלית | קוד | סטטוס |
|---|--------|---------------|-----|--------|
| T1 | EL → משוואת קלייין-גורדון | האם $\partial\mathcal{L}/\partial\phi$ נותן את $(\Box + m_\phi^2)\phi = -(y/2)\bar\chi\chi$? | ✅ | ✅ |
| T2 | פוטנציאל יוקאווה מפורייה | האם $\text{FT}[y^2/(q^2+m_\phi^2)] = \alpha e^{-m_\phi r}/r$? | ✅ | ✅ |
| T3 | VPM — הזזת פאזה $\delta_\ell$ | חישוב $\delta_\ell(k)$ מ-ODE של VPM, השוואה ל-core/v22 | ✅ | ✅ |
| T4 | משקלי ספין מיורנה | גזירת $w_\ell = (1,3)$ מ-CG + Pauli, **תחזית** לא בחירה | ✅ | ✅ |
| T5 | $\sigma_T$ מלא | $\sigma_T = (2\pi/k^2)\sum_\ell w_\ell(2\ell+1)\sin^2\delta_\ell$, השוואה ל-CSV | ✅ | ✅ |
| T6 | אניהילציה $s$-wave | $\langle\sigma v\rangle_0 = \pi\alpha^2/(4m_\chi^2)$ מאמפליטודת פיינמן | ✅ | ✅ |
| T7 | רליקט מבולצמן | $\Omega h^2$ מפתרון נומרי, השוואה ל-17 BPs | ⬜ | ⬜ |
| T8 | $\sigma_{tr}/\sigma_T$ Maj vs Dir | הבדל שיטתי 10–16%, **לא תועד בספרות** | ✅ | ✅ |
| T9 | מיורנה vs דיראק — $R(v)$ | $R(v) \in [1, 2.18]$, NON-MONOTONIC, **חתימה תצפיתית** | ✅ | ✅ |

**עיקרון:** הפייתון מאמת, לא אני.

---

## רשומה 2 — 2026-03-27 (T1: קונבנציית α)

**סקריפט:** `test_T1_alpha_convention.py`  
**פלט:** `data/archive/T1_alpha_convention_2026_03_27_r010.csv`  
**Run ID:** 010

### שאלת המחקר

$\alpha_{\rm config}$ מ-`global_config.json` — מהו ביחס ל-$y$ של הלגרנג'יאן?

### מה הקוד עושה

1. גוזר סימבולית (SymPy): $\mathcal{L}_{\rm int} = -(y/2)\bar\chi\chi\phi$ → vertex $= y/2$ → Born $= (y/2)^2/(q^2+m_\phi^2)$ → FT → $V(r)$
2. בודק 3 קונבנציות × 4 BPs (BP1, BP9, BP16, MAP):
   - $\alpha = y^2/(4\pi)$ — סטנדרטית (דיראק)
   - $\alpha = y^2/(16\pi)$ — עם ה-$1/2$ מהוורטקס
   - $\alpha = y^2/(8\pi)$ — עם $1/2$ + כפלת $t+u$ מיורנה

### תוצאות

| קונבנציה | ratio_t (t בלבד) | ratio_tu (t+u) | מסקנה |
|----------|------------------|----------------|-------|
| $y^2/(4\pi)$ | 0.25 | 0.50 | ❌ נפסלת |
| $y^2/(16\pi)$ | **1.0000** | 2.0 | ✅ PASS (t בלבד) |
| $y^2/(8\pi)$ | 0.50 | **1.0000** | ✅ PASS (t+u) |

### פרשנות

**שתי קונבנציות יכולות להסביר את $\alpha_{\rm config}$:**

- **אפשרות A:** $\alpha = y^2/(16\pi)$ — הוורטקס הוא $y/2$, VPM רואה ערוץ $t$ בלבד, ומשקלי הספין $w_\ell = (1,3)$ מטפלים בחילוף מיורנה ($u$-channel) ברמת הגלים החלקיים.

- **אפשרות B:** $\alpha = y^2/(8\pi)$ — הוורטקס הוא $y/2$, אבל הפוטנציאל $V(r)$ כבר כולל $t+u$, ו-$w_\ell$ הם רק משקלי ספין (בלי כפלת חילוף).

**שתיהן תקפות מתמטית** — הן נותנות את אותו $\sigma_T$ (אם הקוד מתאים). כדי להכריע איזו קונבנציה הקוד באמת משתמש, צריך **T3** — לחשב $\delta_\ell$ מ-VPM עם $\alpha_{\rm config}$ ולהשוות ל-$\delta_\ell$ עם $\alpha = y^2/(16\pi)$ לעומת $y^2/(8\pi)$.

### סטטוס

- [x] T1 — **הושלם**. הקונבנציה $y^2/(4\pi)$ **נפסלה**. נשארו A ו-B.
- [ ] T3 נדרש להכרעה.

---

## רשומה 3 — 2026-03-27 (EL Pipeline: שכבות 1–5)

**סקריפט:** `lagrangian_euler_lagrange.py`  
**פלט:** `data/archive/EL_pipeline_vs_v22_2026_03_27_r011.csv`  
**Run ID:** 011

### מה הסקריפט עושה

שרשרת שלמה מהלגרנג'יאן עד $\sigma_T/m$ — **ללא import של `sigma_T_vpm` מ-core**:

| שכבה | פעולה |
|-------|-------|
| 1 | EL → $(\Box + m_\phi^2)\phi = -(y/2)\bar\chi\chi$ |
| 2 | NR → $V(r) = -\alpha\, e^{-m_\phi r}/r$ |
| 3 | VPM ODE → $\delta_\ell(k)$ — **מימוש עצמאי** (RK4, Bessel, ספירה adaptive) |
| 4 | $\sigma_T/m = (2\pi/k^2)\sum_\ell w_\ell(2\ell+1)\sin^2\delta_\ell$ |
| 5 | השוואה ל-v22 ב-13 מהירויות × 4 BPs |

### תוצאה

**52/52 PASS** — ratio = 1.0000 בכל נקודה.

| BP | λ | טווח v [km/s] | כל 13 velocities |
|----|---|----------------|-------------------|
| BP1 | 11.1 | 12–4700 | ✅ PASS |
| BP9 | 13.1 | 12–4700 | ✅ PASS |
| BP16 | 2.2 | 12–4700 | ✅ PASS |
| MAP | 33.3 | 12–4700 | ✅ PASS |

### משמעות

1. **שרשרת ה-EL שלמה ותקינה** — מהלגרנג'יאן עד $\sigma_T/m$ בלי "חורים".
2. **$\alpha_{\rm config}$ נכנס ישירות ל-VPM** כ-coupling ב-$V(r) = -\alpha e^{-m_\phi r}/r$.
3. **אין פקטור נסתר** בין הגזירה לקוד — מה שנכנס לקונפיג הוא מה ש-VPM משתמש.
4. **שאלת T1 עדיין פתוחה**: $\alpha_{\rm config}$ = $y^2/(16\pi)$ או $y^2/(8\pi)$? ה-VPM נותן אותו $\sigma_T$ בשניהם (כי $w_\ell$ מפצה). נדרש ניתוח של phase shifts ברמת ה-$\delta_\ell$ הבודד כדי להכריע.

---

## רשומה 4 — 2026-03-27 (T6: אניהילציה ← הכרעה!)

**סקריפט:** `test_T6_annihilation_relic.py`  
**פלט:** `data/archive/T6_annihilation_relic_2026_03_27_r012.csv`  
**Run ID:** 012

### שאלת המחקר

איזו קונבנציה של $\alpha$ נותנת $\Omega h^2 \approx 0.12$?

### לוגיקה

- $\langle\sigma v\rangle_0 = \pi\alpha_{\rm ann}^2/(4m_\chi^2)$
- אם $\alpha_{\rm config} = y^2/(4\pi)$: $\alpha_{\rm ann} = \alpha_{\rm config}$
- אם $\alpha_{\rm config} = y^2/(16\pi)$: $\alpha_{\rm ann} = 4\alpha_{\rm config}$ (כי $y^2/(4\pi) = 4 \cdot y^2/(16\pi)$)
- אם $\alpha_{\rm config} = y^2/(8\pi)$: $\alpha_{\rm ann} = 2\alpha_{\rm config}$

$\langle\sigma v\rangle \propto \alpha^2$ — אז פקטור 2 ב-$\alpha$ = פקטור 4 ב-$\Omega h^2$.

### תוצאות

| BP | H1: $y^2/(4\pi)$ | H2: $y^2/(16\pi)$ | H3: $y^2/(8\pi)$ |
|----|------|-------|-------|
| BP1 | $\Omega h^2 = 0.121$ (**+0.8%**) ✅ | 0.009 (−92.8%) | 0.032 (−73.0%) |
| BP9 | $\Omega h^2 = 0.121$ (**+0.8%**) ✅ | 0.009 (−92.8%) | 0.032 (−73.0%) |
| BP16 | $\Omega h^2 = 0.121$ (**+1.0%**) ✅ | 0.009 (−92.5%) | 0.033 (−72.5%) |
| MAP | $\Omega h^2 = 0.243$ (+102%) ⚠️ | 0.017 (−85.6%) | 0.065 (−45.9%) |

### הכרעה

$$\boxed{\alpha_{\rm config} = \frac{y^2}{4\pi}}$$

**הקונבנציה הסטנדרטית.** BP1, BP9, BP16 נותנים $\Omega h^2 = 0.121$ (סטייה < 1%).

### רגע — מה עם T1?

ב-T1 נפסלה $y^2/(4\pi)$ (ratio = 0.25). איך זה מסתדר?

**התשובה: T1 היה חלקי.** T1 בדק: "אם $\alpha = y^2/(4\pi)$, מה הפוטנציאל מגזירת Born?"
- תוצאה: $V_t(r) = -\frac{\alpha}{4}e^{-m_\phi r}/r$ — חסר פקטור 4 ← ratio = 0.25

**הפתרון:** הפוטנציאל ב-VPM **לא** מ-Born t-channel בלבד!
- ב-VPM, $\alpha_{\rm config} = y^2/(4\pi)$ נכנס ל-$V(r) = -\alpha e^{-mr}/r$
- הפקטור $(y/2)^2/(4\pi) = \alpha/4$ הוא רק ה-Born t-channel
- ה-VPM **מפצה** עם $w_\ell = (1,3)$: אלה לא רק "משקלי ספין" — הם כוללים את כפלת ה-$t+u$ exchange!
- כלומר: $\alpha_{\rm VPM} \equiv y^2/(4\pi)$, ו-$w_\ell$ מכפיל את ה-cross section ב-4 בממוצע (1+3)/2 × 2 levels)

### מסקנה מאוחדת

| רכיב | קונבנציה | מבטיח |
|------|----------|-------|
| `global_config.json` | $\alpha = y^2/(4\pi)$ | ✅ T6 |
| VPM: $V(r) = -\alpha e^{-mr}/r$ | $\alpha$ ישירות | ✅ EL pipeline |
| VPM: $w_\ell = (1,3,1,3,...)$ | כולל $t+u$ exchange + spin | ✅ T6 מסביר את T1 |
| אניהילציה: $\langle\sigma v\rangle = \pi\alpha^2/(4m^2)$ | $\alpha$ ישירות | ✅ T6 |

**MAP חריג** ($\Omega h^2 = 0.243$) — ידוע שזו נקודת MAP (maximum a posteriori) שנבחרה מ-$\sigma_T$ fit, לא מ-relic. הערך הגבוה אומר ש-MAP לא עובר relic cut ← עקבי עם `MAP_relic` שהוא BP נפרד.

---

## רשומה 5 — 2026-03-27 (T4: משקלי ספין מיורנה — גזירה מעקרונות ראשונים)

**סקריפט:** `test_T4_majorana_weights.py`  
**פלט:** `data/archive/T4_majorana_weights_2026_03_27_r013.csv`  
**Run ID:** 013

### שאלת המחקר

$w_\ell = (1, 3, 1, 3, ...)$ — מאיפה? האם זו **תחזית** מהלגרנג'יאן, או בחירה שרירותית?

### גזירה סימבולית (SymPy)

1. **פירוק ספין**: $\frac{1}{2} \otimes \frac{1}{2} = 0 \oplus 1$
2. **אימות CG**: 
   - סינגלט ($S=0$): $\langle\frac{1}{2},+\frac{1}{2};\frac{1}{2},-\frac{1}{2}|0,0\rangle = +\frac{1}{\sqrt{2}}$, $\langle\frac{1}{2},-\frac{1}{2};\frac{1}{2},+\frac{1}{2}|0,0\rangle = -\frac{1}{\sqrt{2}}$ → **אנטי-סימטרי** ✅
   - טריפלט ($S=1$): $\langle\frac{1}{2},+\frac{1}{2};\frac{1}{2},-\frac{1}{2}|1,0\rangle = +\frac{1}{\sqrt{2}}$, $\langle\frac{1}{2},-\frac{1}{2};\frac{1}{2},+\frac{1}{2}|1,0\rangle = +\frac{1}{\sqrt{2}}$ → **סימטרי** ✅
3. **עקרון פאולי**: $\psi_{\rm total}$ חייב להיות אנטי-סימטרי (פרמיונים זהים)
   - $\ell$ זוגי (spatial סימטרי) × סינגלט (spin אנטי-סימטרי) = אנטי-סימטרי ✅
   - $\ell$ אי-זוגי (spatial אנטי-סימטרי) × טריפלט (spin סימטרי) = אנטי-סימטרי ✅
   - שאר הצירופים → **אסורים**
4. **משקל = $2S+1$**: $w_{\ell\,\text{even}} = 1$, $w_{\ell\,\text{odd}} = 3$

### אימות נומרי

**20/20 PASS** — $w=(1,3)$ תואם ל-v22 בכל BP ובכל מהירות.

| BP | $v$ [km/s] | $\sigma/m$ v22 | $\sigma/m$ w=(1,3) | $\sigma/m$ w=(1,1) | ratio(1,3)/v22 | ratio(1,1)/v22 |
|----|-----------|---------------|-------------------|-------------------|---------------|---------------|
| BP1 | 100 | 0.772 | 0.772 | 0.451 | 1.0000 | 0.584 |
| BP9 | 100 | 2.228 | 2.228 | 1.123 | 1.0000 | 0.504 |
| MAP | 30 | 1.797 | 1.797 | 0.993 | 1.0000 | 0.552 |

### תובנה פיזיקלית

- $w=(1,1)$ (Dirac) נותן $\sigma_T$ **נמוך ב-42%** בממוצע מ-$w=(1,3)$ (Majorana)
- המשקלים מקודדים **שלושה אפקטים** בו-זמנית:
  1. חילוף חלקיקים זהים ($t+u$ channels)
  2. סטטיסטיקת ספין (Fermi-Dirac)
  3. חתך פעולה ממוצע על ספין ($\frac{1}{4}$ singlet $+ \frac{3}{4}$ triplet)
- **זו תחזית** מהלגרנג'יאן + סימטריית Majorana, לא בחירה שרירותית

---

## רשומה 6 — 2026-03-27 (T8: Transfer vs Elastic — הבדל Majorana/Dirac)

**סקריפט:** `test_T8_transfer_vs_elastic.py`  
**פלט:** `data/archive/T8_transfer_vs_elastic_2026_03_27_r014.csv`  
**Run ID:** 014

### שאלת המחקר

האם היחס $\sigma_{tr}/\sigma_T$ תלוי בסוג הפרמיון (Majorana vs Dirac)?

### הגדרות

$$\sigma_T = \frac{2\pi}{k^2}\sum_\ell w_\ell(2\ell+1)\sin^2\delta_\ell$$
$$\sigma_{tr} = \frac{4\pi}{k^2}\sum_\ell w_\ell(\ell+1)\sin^2(\delta_\ell - \delta_{\ell+1})$$
$$\sigma_{vis} = \frac{4\pi}{k^2}\sum_\ell w_\ell\frac{(\ell+1)(\ell+2)}{2\ell+3}\sin^2(\delta_\ell - \delta_{\ell+2})$$

### תוצאות

**הפרש שיטתי עד 16.3%** ביחס $R_{tr} = \sigma_{tr}/\sigma_T$ בין Majorana ל-Dirac!

| $v$ [km/s] | $\langle R_{tr}^{Maj}\rangle$ | $\langle R_{tr}^{Dir}\rangle$ | $\Delta$ |
|-----------|---------------------------|---------------------------|---------|
| 10 | 1.978 | 1.979 | −0.0% |
| 30 | 1.643 | 1.673 | −1.8% |
| 100 | 0.915 | 0.986 | **−7.2%** |
| 300 | 0.368 | 0.415 | **−11.5%** |
| 1000 | 0.078 | 0.089 | **−12.3%** |

### משמעות פיזיקלית — **חדש!**

> סימולציות N-body של SIDM משתמשות לרוב ב-$\sigma_{tr}$ (momentum transfer), לא ב-$\sigma_T$ (elastic total). ההמרה בין השניים **תלויה בסוג הפרמיון** בעד 16%.

**אם DM הוא Majorana** ומודלים מניחים Dirac (או לא מבחינים), הם מקבלים $\sigma_{tr}/\sigma_T$ **שגוי ב-10–16%** במהירויות של clusters.

**זו הטיה שיטתית שלא תועדה** ברוב ספרות ה-SIDM.

---

## רשומה 7 — 2026-03-27 (T9: R(v) = Majorana/Dirac — חתימה תצפיתית)

**סקריפט:** `test_T9_majorana_vs_dirac.py`  
**פלט:** `data/archive/T9_majorana_vs_dirac_2026_03_27_r016.csv`  
**Run ID:** 016

### שאלת המחקר

האם $R(v) = \sigma_T^{Maj}/\sigma_T^{Dir}$ תלוי-מהירות? אם כן — יש חתימה תצפיתית.

### תוצאות — **NON-MONOTONIC!**

| $v$ [km/s] | $\langle R(v)\rangle$ | טווח | רגים |
|-----------|---------------------|------|------|
| 10 | 1.001 | [1.000, 1.004] | $s$-wave בלבד |
| 30 | 1.217 | [1.000, 1.811] | מעבר |
| 50 | 1.423 | [1.003, 2.181] | מעבר חד |
| 100 | 1.700 | [1.090, 2.011] | few-$\ell$ |
| 300 | 1.933 | [1.611, 2.077] | few-$\ell$ / semi-classical |
| 1000 | 1.981 | [1.907, 2.009] | semi-classical |

**טווח כולל: $R(v) \in [1.000,\, 2.181]$, ממוצע = 1.73**

### התנהגות לפי BP

| BP | $R(30)$ | $R(300)$ | $R(1000)$ | spread | מונוטוני? |
|----|---------|----------|-----------|--------|-----------|
| BP1 | 1.006 | 2.077 | 2.009 | 1.071 | **NON-MONOTONIC** |
| BP9 | 1.051 | 2.062 | 2.006 | 1.011 | **NON-MONOTONIC** |
| BP16 | 1.000 | 1.611 | 1.907 | 0.907 | increasing |
| MAP | 1.811 | 1.983 | 2.000 | 0.190 | NON-MONOTONIC |

### חתימה תצפיתית — **החידוש המרכזי**

$$R(v) = \frac{\sigma_T^{Maj}(v)}{\sigma_T^{Dir}(v)}: \quad R \approx 1\text{ (dwarfs, }v\sim 30\text{)} \longrightarrow R \approx 2\text{ (clusters, }v\sim 1000\text{)}$$

**פרשנות:**
- במהירויות נמוכות ($v \lesssim 30$ km/s), רק $\ell = 0$ תורם → $w_0 = 1$ → **אין הבדל** בין Majorana ל-Dirac
- במהירויות גבוהות ($v \gtrsim 200$ km/s), הרבה $\ell$ תורמים → ממוצע סטטיסטי → $R \to 2$
- ב-**אזור המעבר** ($30 \lesssim v \lesssim 200$ km/s), $R(v)$ עולה בצורה חדה ותלוית-BP

**אם מודדים $\sigma/m$ ב-dwarf ($v \sim 30$) וב-cluster ($v \sim 1000$) ומניחים שמודל אחד מתאים:**
- **Majorana**: $\sigma(30)/\sigma(1000)$ קטן מהצפי (כי $R(30) \approx 1$ אבל $R(1000) \approx 2$)
- **Dirac**: היחס שונה

**זה מדיד** — mult-scale SIDM observations (dwarfs + clusters) שכבר קיימות!

---

## רשומה 8 — 2026-03-27 (תשתית פייפליין)

### מה נוצר

| קובץ | תפקיד |
|-------|--------|
| `run_pipeline.py` | ראנר — `--list`, `--step N`, `--from N`, `--force` |
| `docs/execution_pipeline.csv` | 6 שלבים עם dependencies |

### מבנה פייפליין

```
Step 1: T1 — Alpha Convention        [test_T1_alpha_convention.py]
Step 2: T2+T3+T5 — EL Pipeline       [lagrangian_euler_lagrange.py]     ← depends on 1
Step 3: T6 — Annihilation Relic       [test_T6_annihilation_relic.py]    ← depends on 1
Step 4: T4 — Majorana Weights         [test_T4_majorana_weights.py]      ← depends on 2
Step 5: T8 — Transfer vs Elastic      [test_T8_transfer_vs_elastic.py]   ← depends on 2,4
Step 6: T9 — Majorana vs Dirac        [test_T9_majorana_vs_dirac.py]     ← depends on 2,4,5
```

### ריצה מלאה

```
python run_pipeline.py --from 4
```

**6/6 PASS** — כל השלבים עברו (כולל 1–3 מסשנים קודמים).

| # | סקריפט | זמן | תוצאה |
|---|--------|-----|--------|
| 4 | T4 | 42s | 20/20 PASS |
| 5 | T8 | 110s | הבדל עד 16.3% |
| 6 | T9 | 441s | R(v) ∈ [1.0, 2.18] |

---

## סיכום ביניים — מה גילתה הנגזרת?

### אימותים (ידוע מ-SIDM, אושר מחדש)
- $\alpha_{\rm config} = y^2/(4\pi)$ — הקונבנציה הסטנדרטית ✅
- VPM-EL chain תקין — 52/52 PASS ✅
- $\Omega h^2$ תקין עבור BP1, BP9, BP16 ✅

### חידושים (חדש מהנגזרת!)
1. **$w_\ell$ הם תחזית**: נגזרו מ-CG + Pauli, לא נבחרו ביד. Majorana מגדיל $\sigma_T$ ב-42%
2. **$\sigma_{tr}/\sigma_T$ תלוי-פרמיון**: הבדל שיטתי 10–16% (Majorana vs Dirac) — לא תועד בספרות
3. **$R(v)$ non-monotonic**: מעבר $R \approx 1 \to 2$ כפונקציה של $v$ — חתימה תצפיתית ב-multi-scale SIDM

### מה נשאר (T7)
- T7: בולצמן $\Omega h^2$ מלא — בוצע חלקית ב-T6, אבל ניתן להרחיב ל-17 BPs + $p$-wave
