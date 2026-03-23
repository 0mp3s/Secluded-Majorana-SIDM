# אפשרויות תיקון למודל Secluded Majorana SIDM

**תאריך:** 23 מרץ 2026  
**מחבר:** Claude Opus 4.6

---

## אבחון: שתי הבעיות המרכזיות

### בעיה 1: האניהילציה היא p-wave, לא s-wave
עבור פרמיון **Majorana** עם צימוד סקלרי $\bar{\chi}\chi\phi$, הביטוי $\bar{v}(p_2)u(p_1)$ בגבול הלא-רלטיביסטי שווה ל:

$$\bar{v}(p_2)u(p_1) = -2\eta^\dagger(\vec{\sigma}\cdot\vec{p})\xi \propto v$$

כלומר **נעלם בסף** ($v \to 0$). התוצאה: s-wave = 0, והגל המוביל הוא p-wave.

### בעיה 2: סקלר יציב ב-10 MeV סוגר את היקום
φ עם מסה ~10 MeV שהיה בשיווי משקל תרמי ($\xi \approx 1$) מייצר:
$$\Omega_\phi h^2 \sim \frac{m_\phi}{93\,\text{eV}} \sim 10^5$$
**overclosure קטסטרופלי** — נדרש מנגנון דלדול.

---

## שלוש אפשרויות תיקון

---

## אפשרות A: החלפת Majorana ב-Dirac (מומלצת)

### למה זה עובד

עבור פרמיון **Dirac**, האניהילציה $\chi\bar{\chi} \to \phi\phi$ דרך אותו צימוד סקלרי $y\bar{\chi}\chi\phi$ **היא s-wave**:

$$\langle\sigma v\rangle = \frac{y^4}{32\pi m_\chi^2} = \frac{\pi\alpha^2}{2m_\chi^2}$$

ההבדל הקריטי: ב-$\chi\bar{\chi}$ (חלקיק-אנטי-חלקיק) אין אנטיסימטריזציה של המצב ההתחלתי כמו ב-Majorana $\chi\chi$. הביטויים $\bar{v}\vec{\gamma}u$ שורדים בסף ונותנים תרומת s-wave סופית. זה תוצאה ידועה בספרות (Kumar & Mardon 2013; Gondolo & Gelmini 1991).

### מה משתנה בקוד

| רכיב | שינוי נדרש | מאמץ |
|------|-----------|------|
| **VPM solver** (`v22_raw_scan.py`) | **ללא שינוי** — אותו פוטנציאל Yukawa | 0 |
| **Phase shifts** (δ_l) | **ללא שינוי** — אותו ODE | 0 |
| **σ_T formula** | נדרש שינוי — ממוצע Dirac (ראה למטה) | בינוני |
| **Boltzmann solver** | שינוי מקדם: $\pi\alpha^2/(2m^2)$ במקום $\pi\alpha^2/(4m^2)$ | קל |
| **Relic scan** | הרצה מחדש עם הנוסחה החדשה | בינוני |
| **Paper text** | "Majorana" → "Dirac" בכל מקום | בינוני |

### שינוי נוסחת הפיזור

ביקום סימטרי (Dirac תרמי), חצי מהמפגשים הם $\chi\chi$ (זהים) וחצי $\chi\bar{\chi}$ (שונים):

**$\chi\chi$ (identical fermions):** 
$$\sigma_{\chi\chi} = \frac{2\pi}{k^2}\left[\sum_{l\,\text{even}}(2l+1)\sin^2\delta_l + 3\sum_{l\,\text{odd}}(2l+1)\sin^2\delta_l\right]$$
(בדיוק כמו Majorana — אותו קוד!)

**$\chi\bar{\chi}$ (distinguishable):**
$$\sigma_{\chi\bar{\chi}} = \frac{4\pi}{k^2}\sum_l(2l+1)\sin^2\delta_l$$
(ללא משקולות סימטריה — נוסחה רגילה)

**ממוצע אפקטיבי:**
$$\sigma_{\rm eff}/m = \frac{1}{2}\sigma_{\chi\chi}/m + \frac{1}{2}\sigma_{\chi\bar{\chi}}/m$$

### שינוי בקוד — `sigma_T_vpm`

```python
# Dirac: effective σ/m for symmetric abundance
# (1) Identical-particle sum (same as before):
sigma_identical = 2π/k² × Σ[w_l (2l+1) sin²δ_l]  # w={1,3,1,3,...}
# (2) Distinguishable sum (new):
sigma_distinguish = 4π/k² × Σ[(2l+1) sin²δ_l]     # no spin weights
# (3) Average:
sigma_eff = 0.5 * sigma_identical + 0.5 * sigma_distinguish
```

Phase shifts δ_l **לא משתנים** — אותו VPM solver בדיוק.

### שינוי Boltzmann

```python
# Before (wrong s-wave for Majorana):
# sv0 = π α² / (4 m²)

# After (correct s-wave for Dirac):
sv0 = π α² / (2 m²)   # factor 2 larger
```

ההשלכה: $\alpha_{\rm relic}$ קטן יותר ב-$\sqrt{2}$, כלומר ~$3.5 \times 10^{-4}$ – $3.4 \times 10^{-3}$ במקום $5.5 \times 10^{-4}$ – $4.8 \times 10^{-3}$. **זה נשאר בתוך התחום ה-SIDM-viable!**

### יתרון מרכזי
**~80% מהקוד הקיים נשמר ללא שינוי.** ה-VPM solver, ה-Born validation, ה-literature cross-checks, ה-chi2 framework, ה-predictions — כולם עובדים עם אותם phase shifts. רק נוסחת ה-σ_T ו-Boltzmann צריכים עדכון.

---

## תיקון בעיית ה-overclosure של φ

### אפשרות A1: קניבליזציה (φ³ self-coupling)

הצימוד $y\bar{\chi}\chi\phi$ **כבר שובר** את סימטריית $Z_2$ ($\phi \to -\phi$). לכן מותר להוסיף $\mu_3 \phi^3/3!$ ללגרנז'יאן — וזה אפילו נוצר רדיאטיבית:

$$\mu_3^{\rm rad} \sim \frac{y^3 m_\chi}{16\pi^2}$$

עם $\mu_3$ (tree-level או loop-induced), התהליך $3\phi \to 2\phi$ מדלדל את שפע ה-φ לפני BBN. זהו מנגנון "cannibal dark matter" (Carlson, Machacek & Hall 1992; Pappadopulo et al. 2016).

**Impact על המודל:**
- פרמטר נוסף אחד ($\mu_3$) — סה"כ 4 פרמטרים
- ה-VPM analysis לא משתנה (φ self-coupling לא משפיע על DM-DM scattering)
- צריך להראות שהדלדול מספיק לפני BBN

### אפשרות A2: φ מתפרק לניוטרינו אפל

הוספת פרמיון קל $\nu_D$ עם $m_{\nu_D} \ll m_\phi$ וצימוד $g_\nu \bar{\nu}_D \nu_D \phi$:

$$\tau_\phi = \frac{8\pi}{g_\nu^2 m_\phi}$$

עבור $\tau < 1$ שנייה (לפני BBN): נדרש $g_\nu > 5 \times 10^{-12}$ — צימוד זעיר.

**Impact:** 2 פרמטרים נוספים ($m_{\nu_D}$, $g_\nu$); ה-$\nu_D$ תורם $\Delta N_{\rm eff} \sim 0.05$ (מתחת ל-Planck).

### אפשרות A3: ξ << 1 (התנתקות מוקדמת)

אם ה-dark sector התנתק מה-SM ב-$T \gg 100$ GeV, ההרבה phase transitions של ה-SM (QCD, EW) מחממים את ה-SM ביחס ל-dark sector:

$$\xi = \left(\frac{g_{*S}(T)}{g_{*S}(T_{\rm dec})}\right)^{1/3} \sim 0.3\text{–}0.5$$

זה מפחית את $n_\phi$ ב-$\xi^3$, אבל **לא מספיק** — עדיין overclosure ב-$\sim 10^3$. צריך לשלב עם A1 או A2.

---

## אפשרות B: להישאר עם Majorana + לקבל p-wave

### התוצאה

עם p-wave: $\langle\sigma v\rangle = b v^2$ כאשר $b = 3\alpha^2/(16m_\chi^2)$.

ב-freeze-out ($v^2 \sim 6/x_f \sim 0.3$): 
$$\langle\sigma v\rangle_{\rm fo} \approx 0.3 \times \frac{3\alpha^2}{16m_\chi^2}$$

לעומת s-wave:
$$\langle\sigma v\rangle_{\rm fo} = \frac{\pi\alpha^2}{4m_\chi^2}$$

היחס: $\frac{p\text{-wave}}{s\text{-wave}} = \frac{0.3 \times 3/16}{\pi/4} \approx 0.072$

כלומר $\alpha_{\rm relic}$ צריך לעלות ב-$\sqrt{1/0.072} \approx 3.7\times$.

### ההשלכה על ה-island

| | s-wave (הנוכחי) | p-wave (מתוקן) |
|---|---|---|
| $\alpha_{\rm relic}$ range | $5\times10^{-4}$ – $5\times10^{-3}$ | $\sim 2\times10^{-3}$ – $2\times10^{-2}$ |
| SIDM viable range | $1.8\times10^{-6}$ – $3.3\times10^{-3}$ | שינוי קטן (פיזור לא תלוי באניהילציה) |
| **חפיפה?** | רחבה (17 BPs) | **צרה מאוד — רק ב-$m_\chi \gtrsim 80$ GeV** |

ההתנגשות: ב-p-wave, הצימוד הנדרש ל-relic ($\alpha \sim 10^{-2}$) גבוה מהתחום ה-SIDM-viable ($\alpha \lesssim 3 \times 10^{-3}$) עבור רוב המסות. **"האי" מצטמצם דרמטית**, אבל ייתכן שנשארת פינה צרה ב-$m_\chi \sim 100$ GeV.

### יתרון
אם עובד — זה כנה ומדויק. לא צריך לשנות מודל.

### חיסרון
הרעיון המרכזי ("exact overlap between relic and SIDM") נהרס. המאמר מאבד את ה-USP שלו.

---

## אפשרות C: Complex Scalar Mediator (Φ מרוכב)

### הרעיון

במקום φ ממשי, להשתמש בסקלר מרוכב Φ עם סימטריית U(1) גלובלית:

$$\mathcal{L} = -\frac{y}{2}\bar{\chi}\chi(Φ + Φ^\dagger)$$

האניהילציה $\chi\chi \to Φ Φ^*$ — **החלקיקים הסופיים שונים** (Φ ≠ Φ*), אז ויכוח ה-parity/הסימטריה לחלקיקים זהים **לא חל**.

### בעיה
- 2 דרגות חופש במקום 1 → $\Delta N_{\rm eff}$ כפול
- ה-overclosure מוכפל
- צריך סימטריית gauge (מוסיף מורכבות)
- **לא ברור שה-s-wave שורד** גם כאן (צריך חישוב מפורש)

**לא מומלץ** — מוסיף מורכבות בלי לפתור את בעיית ה-overclosure.

---

## המלצה סופית

### הנתיב המומלץ: **Dirac + Cannibalization** (אפשרות A + A1)

```
V10 (Majorana, s-wave wrong)  →  V11 (Dirac, s-wave correct, φ³ depletion)
```

**הגדרת המודל V11:**

$$\mathcal{L}_{\rm dark} = \bar{\chi}(i\not{\partial} - m_\chi)\chi - y\bar{\chi}\chi\phi + \frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2 - \frac{\mu_3}{3!}\phi^3$$

**4 פרמטרים:** $(m_\chi, m_\phi, \alpha = y^2/4\pi, \mu_3)$

**כל היתרונות נשמרים:**
- ✅ VPM solver — ללא שינוי (אותו פוטנציאל Yukawa)
- ✅ s-wave annihilation — נכון עבור Dirac
- ✅ Secluded dark sector — אותו רעיון
- ✅ Higgs portal excluded — אותו ארגומנט
- ✅ 3-parameter SIDM scan — אותו scan (+ $\mu_3$ חדש)

**מה שנפתר:**
- ✅ Annihilation s-wave — מוכח בספרות עבור Dirac
- ✅ φ overclosure — $3\phi \to 2\phi$ cannibalization
- ✅ $\Delta N_{\rm eff}$ — φ מתדלדל, לא overclosure

### רשימת משימות

1. **אנליטי:** הוכחה מפורשת ש-$\chi\bar{\chi} \to \phi\phi$ הוא s-wave עבור Dirac (מספיק לצטט Kumar & Mardon 2013)
2. **קוד — σ_T:** להוסיף פונקציית Dirac effective averaging ב-`sigma_T_vpm`
3. **קוד — Boltzmann:** לעדכן $\langle\sigma v\rangle = \pi\alpha^2/(2m_\chi^2)$
4. **קוד — scan:** להריץ מחדש את ה-scan עם הנוסחאות החדשות
5. **אנליטי — φ³:** לחשב את קצב הקניבליזציה ולהראות דלדול לפני BBN
6. **Paper:** לשכתב §§2, 3.2, 5.3, 6.1, 6.2, 8 — Majorana → Dirac, + section על φ depletion
7. **תיקוני באגים:** NFW ρ_s/3, λ convention, CSV column names

---

## סיכום

| אפשרות | s-wave? | Overclosure? | שימור קוד | מורכבות | המלצה |
|---------|---------|-------------|-----------|---------|-------|
| **A: Dirac + φ³** | ✅ | ✅ | ~80% | 4 params | **⭐ מומלצת** |
| B: Majorana p-wave | ❌ (p-wave) | ❌ | ~60% | 3 params | אפשרי אבל חלש |
| C: Complex scalar | ❓ | ❌ | ~70% | 4+ params | לא מומלץ |
