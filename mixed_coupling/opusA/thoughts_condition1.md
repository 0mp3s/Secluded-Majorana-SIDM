# מחשבות OpusA — הדרך לנוסחה $a_0 = y_s^2 y_p^2 / (8\pi m_\chi^2)$

*23 מרץ 2026*

---

## נקודת המוצא: ממצא FATAL

ביקורת העמיתים הגיעה עם משפט ברור: "Majorana + pure scalar coupling → p-wave only".
זה נכון. כתמיד בפיזיקה, השאלה היא *למה* זה נכון — כי "למה" פותח דלת.

## התובנה הראשונה: CP selection rule

כשחשבתי על זה סימטרית, ה-argument הוא:

1. שני פרמיוני Majorana זהים ב-s-wave ($L=0$) חייבים להיות ב-antisymmetric spin state → $S=0$ → ¹S₀.
2. Parity: $P = (-1)^{L+1} = -1$.
3. Charge conjugation: $C = (-1)^{L+S} = +1$.
4. לכן **CP = −1** עבור מצב ¹S₀.
5. שני סקלרים אמיתיים זהים ב-s-wave ($L'=0$): CP = +1 תמיד, ללא תלות ב-CP parity של $\phi$.
6. → מעבר מ-CP = −1 ל-CP = +1 **דורש שבירת CP**.

זו הנקודה. זה לא "Majorana can't do s-wave" — זה "CP-conserving Majorana can't do s-wave to identical scalars". ברגע ש-CP נשבר, החסימה מוסרת.

## מה שובר CP?

ה-coupling הכי כללי של פרמיון Majorana לסקלר אמיתי (עם $Z_2$: $\chi \to -\chi$):

$$\Gamma = y_s \mathbb{1} + i y_p \gamma_5$$

- $y_s$ לבד: CP conserving → $a_0 = 0$.
- $y_p$ לבד: גם CP conserving (רק parity flip) → $a_0 = 0$.
- $y_s \neq 0$ **ו-** $y_p \neq 0$ **בו-זמנית**: CP violated! הגורם המתערב $y_s y_p$ הוא CP-odd.

זה לא מסתורי — זו אותה פיזיקה של EDM (Electric Dipole Moment). EDM דורש גם scalar וגם pseudoscalar coupling.

## בניית הסקריפט: מספרים, לא אלגברה

היו לי שתי אפשרויות:
1. אלגברת טרייסים ידנית — 2 עמודי חישוב, סיכוי לשגיאות.
2. מטריצות 4×4 נומריות — computer does the algebra.

בחרתי 2. יותר אמין, ויש cross-check מובנה: אם $a_0(y_s, 0) \neq 0$ מספרית, יש באג בקוד.

### האלגוריתם:

1. הגדר gamma matrices (Dirac rep), γ₅.
2. Self-test: Clifford algebra, $\gamma_5^2 = 1$, $\{\gamma_5, \gamma^\mu\} = 0$.
3. לכל $(y_s, y_p)$: חשב $\Sigma|M|^2$ כפונקציה של $v$ ב-CM frame.
4. Fit $\sigma v(v) = a_0 + a_1 v^2$ → חלץ $a_0$.

### הבאג הראשון: סימן מינוס

בגרסה הראשונה כתבתי $M = M_t - M_u$ (מינוס, כנהוג ל-fermion exchange).
אבל כאן ה-exchange הוא של **בוזונים** ($k_1 \leftrightarrow k_2$, לא $p_1 \leftrightarrow p_2$).
שני $\phi$ זהים → Bose symmetry → $M = M_t + M_u$ (פלוס).

עם מינוס, pure scalar/pseudoscalar **לא** התאפסו — סימן ברור שיש שגיאה.
עם פלוס, pure scalar/pseudoscalar **התאפסו** ל-$\sim 10^{-11}$ (שגיאת עיגול) — כצפוי מ-CP.

זהו self-consistency check חזק: **הפיזיקה של CP מכתיבה את הסימן.**

## הנוסחה

לאחר שכל הבדיקות עברו, הוצאתי את המקדם:

$$a_0(1,1) = 9.2948 \times 10^{-5} \text{ GeV}^{-2}$$
$$\frac{1}{8\pi m_\chi^2} = 9.2948 \times 10^{-5} \text{ GeV}^{-2}$$

**התאמה ל-$5 \times 10^{-7}$ %.**

לכן:
$$\boxed{a_0 = \frac{y_s^2 y_p^2}{8\pi m_\chi^2}}$$

או, במונחי coupling constants $\alpha_i = y_i^2/(4\pi)$:

$$\boxed{a_0 = \frac{2\pi \alpha_s \alpha_p}{m_\chi^2}}$$

## למה $1/(8\pi)$ ולא $1/(16\pi)$?

הנוסחה ה"נאיבית" ל-Majorana scalar (pure $y_s$) היא:
$$\sigma v = \frac{y_s^4}{16\pi m_\chi^2} \quad \text{(p-wave)}$$

ה-$a_0$ שלנו הוא **כפול** מזה (עם $y_s^2 y_p^2$ במקום $y_s^4$). הסיבה: ה-cross term בין הדיאגרמות t ו-u עם mixed coupling מתנהג אחרת. עם פלוס (Bose), הביטויים ה-CP-violating מתחברים constructively.

## Scaling test: 9/9 נקודות

זו הבדיקה הכי חזקה. לא רק ש-$a_0 \propto y_s^2 y_p^2$ — אלא שהיא **בדיוק** $c \cdot y_s^2 y_p^2$ עם אותו $c$ בכל 9 הנקודות:

```
(1,1)    ratio = 1.00000
(2,1)    ratio = 1.00000
(1,2)    ratio = 1.00000
(2,2)    ratio = 1.00000
(3,1)    ratio = 1.00000
(1,3)    ratio = 1.00000
(0.5,0.5) ratio = 1.00000
```

Pure scalar ו-pure pseudoscalar:
```
(1,0)  a₀ = 5.7e-11  →  0 (CF: a₀_ref = 9.3e-5)
(0,1)  a₀ = 1.7e-12  →  0
```

## המשמעות הפיזיקלית

**שבירת CP ב-dark sector היא תכונה, לא באג.** 

ב-Standard Model, CP נשבר ב-CKM matrix. אין סיבה שב-dark sector CP תשמר. ה-"most general Lagrangian" principle אומר שכל coupling שלא נאסר ע"י סימטריה — קיים.

למודל יש 5 פרמטרים: $(m_\chi, m_\phi, \alpha_s, \alpha_p, \mu_3)$ במקום 4 של Dirac. אבל הוא מציע תופעה חדשה — CP violation ב-dark sector — שיכולה להבדל ב-next-generation experiments.

## מה עכשיו?

תנאי 1 **עבר**. הנוסחה מאומתת.

הרליק דנסיטי נמצא על היפרבולה $\alpha_s \alpha_p = \text{const}$ — ולא fine-tuning בכלל, כי כל הנקודות על ההיפרבולה נותנות את אותו $\Omega_\chi h^2$.

תנאי 4 הבא: אשר ש-VPM (SIDM cross section) תלוי רק ב-$\alpha_s$ ו-$y_p$ contribution ל-NR potential הוא negligible כ-$m_\phi^2/m_\chi^2 \sim 10^{-7}$.
