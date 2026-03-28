# מתודולוגיה — אינטגרל המסלול של הלגרנג'יאן

## הבעיה

נתון הלגרנג'יאן של מודל SIDM מג'ורנה מבודד:

$$\mathcal{L} = \frac{1}{2}\bar{\chi}(i\partial\!\!\!/ - m_\chi)\chi + \frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2 - \frac{y}{2}\bar{\chi}\chi\phi$$

**השאלה:** מה ה-path integral עושה עם הלגרנג'יאן הזה? כלומר — מהם כל הצעדים מ-$Z[J,\eta,\bar{\eta}]$ (הפונקציונל המחולל) לכל הנצפים הפיזיקליים?

## העיקרון המנחה

> **הפייתון מאמת, לא אני.**

כל שכבה של אינטגרל המסלול מחושבת נומרית ומושוותה לערכים ידועים.
אין ערכים מקודדים קשיח — הכל נטען מ-`global_config.json` (פרמטרים גלובליים) או `data/config.json` (בחירות פוסט-אינטגרציה).

## שבע השכבות

### שכבה 1 — אינטגרציה גאוסית (שדות חופשיים)

$$Z_0[J,\eta,\bar{\eta}] = \int\!\mathcal{D}\phi\,\mathcal{D}\chi\,\mathcal{D}\bar{\chi}\;\exp\!\left(i\int d^4x\,\mathcal{L}_{\rm free} + \text{sources}\right)$$

השלמת ריבוע בצורה הריבועית ← מפיצי פיינמן:
- סקלרי: $\Delta_F(q^2) = i/(q^2 - m_\phi^2 + i\epsilon)$
- פרמיוני: $S_F(p) = i(\not{p} + m_\chi)/(p^2 - m_\chi^2 + i\epsilon)$
- קודקוד: $-iy/2$ (מג'ורנה יוקאווה)

**בדיקה:** טווח הכוח $r_0 = \hbar c / m_\phi$ ואורך גל דה-ברויי $\lambda_{\rm dB}$.

### שכבה 2 — רמת עץ (Born)

שני קודקודים × מפיץ אחד = פוטנציאל יוקאווה:

$$V(r) = -\frac{\alpha}{r}\,e^{-m_\phi r}, \qquad \alpha = \frac{y^2}{4\pi}$$

**טענה:** קרוב-שדה Born נכשל כי $\lambda = \alpha m_\chi / m_\phi \gg 1$.

### שכבה 3 — דטרמיננט הפלוקטואציות (VPM)

$$\sigma_T = \frac{2\pi}{k^2}\sum_\ell w_\ell(2\ell+1)\sin^2\delta_\ell$$

- $w_\ell = 1$ (זוגי), $w_\ell = 3$ (אי-זוגי) — סימטריית חלקיקים זהים (מג'ורנה)
- הזזות פאזה $\delta_\ell$ הן **בדיוק** ה-log-determinant של אינטגרל המסלול!
- VPM פותר את משוואת שרדינגר הרדיאלית עם פוטנציאל היוקאווה ← $\delta_\ell$ מדויקים

**בדיקה:** Born/VPM $\approx 88\times$–$135\times$ → Born לא שמיש, VPM חובה.

### שכבה 4 — קולמן-ויינברג (T=0)

אינטגרציה על הפרמיון $\chi$ ברקע $\sigma$:

$$V_{\rm CW}(\theta) = -\frac{n_f}{64\pi^2}M^4(\theta)\!\left[\ln\frac{M^2(\theta)}{\mu^2} - \frac{3}{2}\right] + \text{scalar}$$

- $n_f = 2$ (מג'ורנה — שתי דרגות חופש ממשיות)
- $M_{\rm eff}(\theta) = m_\chi\sqrt{\cos^2\theta + \sin^2\theta/K_{A_4}}$
- $K_{A_4} = 9$ מקבוצת $A_4$ (Clebsch-Gordan)

### שכבה 5 — פוטנציאל תרמי (מטסובארה)

$$V_T(\theta,T) = V_{\rm CW}(\theta) + \frac{T^4}{2\pi^2}\left[J_B(m_\phi/T) + n_f\,J_F(M/T)\right]$$

בדיקה: עד $T = 10 m_\chi$ — המינימום לא זז מ-$\theta_{A_4}$.

### שכבה 6 — אינטגרל מסלול אוקלידי (אינסטנטון)

$$\Gamma \sim e^{-S_E}, \qquad S_E \sim (f/m_\sigma)^2 \sim 10^{121}$$

→ **יציבות מוחלטת** — קצב מנהור אפסי.

### שכבה 7 — אנניהילציה + שריד

$$\langle\sigma v\rangle_{s\text{-wave}} = \frac{\pi\alpha^2}{4m_\chi^2}$$

$$\Omega_\chi h^2 \approx 0.107\text{–}0.121$$ (Boltzmann → Kolb-Turner)

## חיבור לאנרגיה אפלה

$$V_{\rm DE}(\sigma) = \Lambda_d^4\,(1 - \cos(\sigma/f))$$

$$m_\sigma = \Lambda_d^2/f$$ (GMOR)

$$H_0 = \sqrt{8\pi V_{\rm eff}/(3M_{\rm Pl}^2)} \approx 67.4\;\text{km/s/Mpc}$$

## מבנה הקבצים

| קובץ | תפקיד |
|-------|--------|
| `lagrangian_path_integral.py` | Pipeline ראשי — כל 7 השכבות |
| `data/config.json` | פרמטרי פוסט-אינטגרציה (לא ב-GC) |
| `data/archive/` | גרסאות CSV עם חותמת זמן |
| `docs/execution_pipeline.csv` | תיעוד שלבי הפייפליין |
| `docs/methodology.md` | מסמך זה |
| `output/` | פלט גרפי (PNG) |
| `research_journal.md` | יומן מחקר מפורט |
