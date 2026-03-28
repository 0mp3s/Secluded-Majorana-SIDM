# מתודולוגיה — חקירת הנגזרת של הלגרנג'יאן

## הבעיה

הלגרנג'יאן של המודל:

$$\mathcal{L} \supset -\frac{y}{2}\bar\chi\chi\,\phi$$

קוד ה-VPM (v22) מחשב פיזור עם פוטנציאל:

$$V(r) = -\alpha \frac{e^{-m_\phi r}}{r}$$

כאשר $\alpha$ נלקח **ישירות** מ-`global_config.json`.

**השאלה:** מהו הקשר המדויק בין $y$ (coupling בלגרנג'יאן) ל-$\alpha$ (פרמטר בקוד)?

---

## למה זה חשוב

### הגזירה מ-Feynman rules

מהלגרנג'יאן, ה-vertex factor הוא $-iy/2$ (הפקטור $1/2$ מיורנה).

**אמפליטודת Born**, ערוץ $t$ בלבד:

$$\mathcal{M}_t = \frac{(y/2)^2}{|\mathbf{q}|^2 + m_\phi^2}$$

**Fourier transform** למרחב ממשי:

$$V_t(r) = -\frac{(y/2)^2}{4\pi}\frac{e^{-m_\phi r}}{r} = -\frac{y^2}{16\pi}\frac{e^{-m_\phi r}}{r}$$

### אבל — מיורנה = חלקיקים זהים

חלקיקי מיורנה הם **זהים** ($\chi = \chi^c$), לכן יש גם ערוץ $u$:

$$\mathcal{M} = \mathcal{M}_t \pm \mathcal{M}_u$$

(סימן ± תלוי בערוץ ספין). ה-**פוטנציאל האפקטיבי** עשוי לקבל פקטור 2 מזה.

### שלוש קונבנציות אפשריות

| קונבנציה | הגדרת $\alpha$ | $V(r)$ | תואם לקוד? |
|-----------|----------------|--------|-------------|
| סטנדרטית | $\alpha = y^2/(4\pi)$ | $-(\alpha/4) e^{-m_\phi r}/r$ | ❌ חסר 1/4 |
| עם vertex $1/2$ | $\alpha = y^2/(16\pi)$ | $-\alpha\, e^{-m_\phi r}/r$ | ✅ (t בלבד) |
| + כפלת Majorana | $\alpha = y^2/(8\pi)$ | $-\alpha\, e^{-m_\phi r}/r$ | ✅ (t+u) |

**לא ניתן לדעת מ"ראש"** — צריך לגזור ולבדוק.

---

## שיטת העבודה

### עיקרון

> **הפייתון מאמת, לא אני.**  
> כל טענה פיזיקלית עוברת דרך קוד שרץ ומוציא תוצאה מספרית.

### תשתית (כמו SIDM)

| רכיב | תפקיד |
|-------|--------|
| `data/config.json` | פרמטרים לבדיקות — שואב מ-global_config דרך `GC` |
| `data/archive/` | כל פלט עם חותמת זמן, לעולם לא נדרס |
| `docs/runs_log.csv` | לוג אוטומטי של כל ריצה |
| `core/` | תשתית (import מ-SIDM core) |

### כלל ברזל

- אין hardcoded values — הכל מקונפיג
- אין דריסת קבצים — timestamped output
- כל ריצה מתועדת ב-runs_log.csv
- סקריפט עוקב קורא עם `get_latest()` — לא שם קובץ

---

## רשימת בדיקות T1–T9

| # | בדיקה | שאלה פיזיקלית |
|---|--------|---------------|
| T1 | EL + Born → V(r) | מהו $\alpha_{\rm derived}$ ביחס ל-$\alpha_{\rm config}$? |
| T2 | Fourier של Yukawa | $\text{FT}[y^2/(q^2+m_\phi^2)] \stackrel{?}{=} \alpha e^{-m_\phi r}/r$ |
| T3 | VPM phase shifts | חישוב $\delta_\ell(k)$ מ-ODE, השוואה ל-v22 |
| T4 | משקלי ספין מיורנה | $w_\ell = (1,3,1,3,...)$ — מאיפה? |
| T5 | $\sigma_T$ מלא | $\sigma_T = (2\pi/k^2)\sum_\ell w_\ell(2\ell+1)\sin^2\delta_\ell$ |
| T6 | אניהילציה $s$-wave | $\langle\sigma v\rangle_0 = \pi\alpha^2/(4m_\chi^2)$ |
| T7 | רליקט בולצמן | $\Omega h^2$ מפתרון נומרי |
| T8 | $\sigma_T^{\rm tr} = \sigma_{\rm el}$ | שוויון עבור מיורנה |
| T9 | Majorana vs Dirac $R(v)$ | $R(v) = \sigma_T^{\rm Maj}/\sigma_T^{\rm Dir}$ |

---

## T1 — פירוט

### קלט
- BP מהקונפיג: $(m_\chi, m_\phi, \alpha)$

### חישוב (SymPy סימבולי)
1. $\mathcal{L}_{\rm int} = -(y/2)\bar\chi\chi\phi$ → vertex factor $= -iy/2$
2. Born amplitude: $\mathcal{M}_t = (y/2)^2 / (q^2 + m_\phi^2)$
3. NR Fourier: $V_t(r) = -\frac{1}{4\pi}\frac{(y/2)^2}{1} \cdot \frac{e^{-m_\phi r}}{r}$
4. מזהים: $\alpha_t = y^2/(16\pi)$
5. עם Majorana t+u: $\alpha_{t+u} = y^2/(8\pi)$

### פלט
- CSV עם: BP, $\alpha_{\rm config}$, $y_{\rm from\_config}$ (בשתי קונבנציות), $\alpha_{\rm derived}$, ratio, PASS/FAIL
- פלט מודפס לקונסול
- נשמר ב-`data/archive/` עם timestamp

### קריטריון הצלחה
אם ratio = 1.0 בדיוק (סימבולי) — הקונבנציה מזוהה.
