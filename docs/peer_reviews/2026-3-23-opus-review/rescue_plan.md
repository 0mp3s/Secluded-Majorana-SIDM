# תוכנית הצלה — Secluded Majorana SIDM

**תאריך:** 23 מרץ 2026  
**הקשר:** מענה לביקורת עמיתים (Opus, GPT-5.4, Gemini 3.1 Pro)  
**מטרה:** תיקון שתי הבעיות ה-FATAL ושדרוג המודל לגרסה עמידה

---

## מצב נוכחי — שורש הבעיה

שתי בעיות קריטיות שזוהו ע"י שלושה רפרים בלתי-תלויים:

1. **s-wave מתאפס** — לפרמיוני מאיורנה זהים עם צימוד סקלרי $\bar\chi\chi\phi$, ה-s-wave של $\chi\chi \to \phi\phi$ אסור ע"י שימור זוגיות ($0^{-+} \not\to 0^{++}$)
2. **φ overclosure** — סקלר יציב $m_\phi \sim 10$ MeV שהיה תרמי → $\Omega_\phi h^2 \sim 10^4$ (overclosure פי $10^5$)

---

## אסטרטגיה 1: Pure p-wave + Rescan

### הרעיון
קבל שה-annihilation היא p-wave ותסרוק מחדש.

### הנוסחה
$$\langle\sigma v\rangle_{\rm p\text{-}wave} = \frac{3\pi\alpha^2}{2m_\chi^2} \cdot \frac{v^2}{c^2}$$

ב-freeze-out ($x_f \approx 25$): $\langle v^2\rangle = 6/x_f \approx 0.24$, לכן:

$$\langle\sigma v\rangle_{\rm fo} \approx \frac{3\pi\alpha^2}{2m_\chi^2} \times 0.24 \approx \frac{0.36\pi\alpha^2}{m_\chi^2}$$

### השפעה על $\alpha_{\rm relic}$

| | s-wave (ישן) | p-wave (חדש) | יחס |
|---|---|---|---|
| $\langle\sigma v\rangle_{\rm fo}$ | $\frac{\pi\alpha^2}{4m_\chi^2}$ | $\frac{0.36\pi\alpha^2}{m_\chi^2}$ | $\times 1.44$ |
| $\alpha_{\rm relic}$ | $\alpha_0$ | $\alpha_0 / \sqrt{1.44} \approx 0.83\alpha_0$ | **ירד!** |

> **תוצאה מפתיעה**: אם הפרה-פקטור הנכון הוא $3\alpha^2/(2m_\chi^2)$, ה-$\alpha$ הנדרש **קטן יותר** מהמקרה ה-s-wave, לא גדול יותר.

**אבל** — הפרה-פקטור המדויק תלוי בפרטי הדיאגרמה (t-channel + u-channel עם propagator מסיבי $m_\chi$). צריך לחשב את $a_1$ (מקדם ה-p-wave) מפורשות.

### שלבי ביצוע
1. חשב $\mathcal{M}(t) + \mathcal{M}(u)$ עבור $\chi\chi \to \phi\phi$ עם צימוד $y\bar\chi\chi\phi$
2. פתח ל-threshold: $|\mathcal{M}|^2 = a_0 + a_1 v^2 + ...$
3. אשר $a_0 = 0$, חלץ $a_1$ מדויק
4. עדכן v27_boltzmann_relic.py: `sigma_v = a1 * 6.0 / x`
5. הרץ smart_scan מחדש עם p-wave relic density
6. חפש חפיפה חדשה עם SIDM viability

### יתרונות
- שינוי מינימלי — אותו לגרנז'יאן, רק relic density מתעדכן
- ה-VPM solver **לא משתנה** (חתך פיזור אלסטי לא מושפע)
- ה-Higgs portal exclusion **לא משתנה**

### חסרונות
- תלוי ב-$a_1$ המדויק — אם הוא קטן, $\alpha_{\rm relic}$ עולה ויכול לצאת מטווח ה-SIDM
- CMB constraints on p-wave: Planck מגבילה p-wave פחות מ-s-wave (יתרון), אבל indirect detection משתנה

---

## אסטרטגיה 2: Mixed Scalar-Pseudoscalar Coupling (מומלצת)

### הרעיון
שנה את הצימוד ללגרנז'יאן מעורב:

$$\mathcal{L} \supset -\frac{1}{2}\bar\chi(y_s + iy_p\gamma_5)\chi\,\phi$$

כאשר $y_s$ (סקלרי) שולט בפוטנציאל יוקאווה, ו-$y_p$ (פסאודו-סקלרי) מאפשר s-wave annihilation.

### למה זה עובד

**אניהילציה ($\chi\chi \to \phi\phi$):**
- הצימוד $iy_p\gamma_5$ הוא CP-odd → המצב ההתחלתי ${}^1S_0$ עם $J^{PC} = 0^{-+}$ **יכול** לעבור ל-$0^{++}$ דרך שני קודקודים pseudoscalar (מכפלה של שני CP-odd = CP-even)
- ה-s-wave שורד: $\langle\sigma v\rangle \propto y_p^4$

**פיזור אלסטי ($\chi\chi \to \chi\chi$):**
- הפוטנציאל יוקאווה ב-Born limit:
  - $y_s$ → אטרקטיבי (Yukawa רגיל)
  - $y_p$ → רפולסיבי (spin-dependent, מדוכא ב-NR limit)
- ב-$v \to 0$: הרכיב $y_s$ **שולט**, $y_p$ מדוכא פי $v^2$
- ה-VPM solver עם $\alpha = y_s^2/(4\pi)$ **תקף ללא שינוי**

### פרמטריזציה

הגדר:
$$\alpha_s = \frac{y_s^2}{4\pi}, \quad \alpha_p = \frac{y_p^2}{4\pi}, \quad \alpha_{\rm tot} = \alpha_s + \alpha_p, \quad \tan\beta = \frac{y_p}{y_s}$$

| תהליך | שלט ע"י | נוסחה |
|--------|----------|-------|
| SIDM ($\sigma/m$) | $\alpha_s$ | VPM עם $\lambda = \alpha_s m_\chi / m_\phi$ |
| Relic density | $\alpha_p$ (s-wave) + $\alpha_s$ (p-wave) | $\langle\sigma v\rangle \approx \pi\alpha_p^2/(4m_\chi^2)$ |
| Higgs portal | $\alpha_s + \alpha_p$ | ללא שינוי (exclusion עדיין תקף) |

### דוגמה מספרית — BP1 מעודכן

**נוכחי:** $\alpha = 1.048 \times 10^{-3}$, $m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV

**$\alpha_{\rm relic}$ הנדרש עבור s-wave:**
$$\alpha_p = \sqrt{\frac{4m_\chi^2 \langle\sigma v\rangle_{\rm fo}}{\pi}} \approx 1.05 \times 10^{-3}$$

**בחירה לדוגמה:** $\alpha_s = 1.05 \times 10^{-3}$ (SIDM), $\alpha_p = 1.05 \times 10^{-3}$ (relic), $\tan\beta = 1$

- SIDM: $\lambda = \alpha_s m_\chi / m_\phi = 1.91$ → **בדיוק אותו VPM כמו קודם**
- Relic: $\langle\sigma v\rangle_0 = \pi \times (1.05 \times 10^{-3})^2 / (4 \times 20.69^2)$ → מתאים ל-$\Omega h^2 = 0.12$
- הנקודות הקיימות **נשמרות** עם $\alpha_s \approx \alpha_{\rm old}$

### שלבי ביצוע
1. חשב $|\mathcal{M}|^2$ עבור $\chi\chi \to \phi\phi$ עם $y_s + iy_p\gamma_5$ — אשר שה-$y_p^4$ term שורד ב-s-wave
2. חשב NR potential עם mixed coupling — אשר ש-$y_p$ contribution הוא velocity-suppressed
3. עדכן v27: `sigma_v = pi * alpha_p**2 / (4 * m_chi**2)`
4. עדכן smart_scan: סרוק ב-$(\alpha_s, \alpha_p)$ במקום $\alpha$ יחיד
5. הרץ VPM עם $\alpha = \alpha_s$ בלבד
6. חפש חפיפה: $\sigma/m(30)$ viable + $\Omega h^2 = 0.12$

### יתרונות
- **שומר על כל 17 נקודות ה-benchmark** (עם $\alpha_s = \alpha_{\rm old}$)
- פרמטר חופשי אחד נוסף ($\tan\beta$) — טבעי פיזיקלית
- ה-VPM solver לא משתנה כלל
- s-wave → Planck + Fermi-LAT constraints ידועים ומנוהלים

### חסרונות
- צימוד מעורב פחות מינימלי מסקלר טהור
- צריך להצדיק למה $y_s \sim y_p$ (naturalness argument)

---

## פתרון בעיית ה-φ Overclosure

### הבעיה
סקלר יציב $m_\phi \sim 10$ MeV בשיווי משקל תרמי → $\Omega_\phi h^2 \sim 10^4$

### פתרון א': Cannibal Mechanism ($\mu_3\phi^3$)

הצימוד $\bar\chi\chi\phi$ **כבר שובר** $Z_2$ ($\phi \to -\phi$), לכן $\mu_3\phi^3$ **מותר ברנורמליזציה** גם אם לא נכלל ב-tree level. הוסף ללגרנז'יאן:

$$\mathcal{L} \supset -\frac{\mu_3}{3!}\phi^3 - \frac{\lambda_4}{4!}\phi^4$$

**תהליך dominanit:** $3\phi \to \phi$ (cannibal)

$$\Gamma_{3\to 1} \sim \frac{\mu_3^4}{(4\pi)^3 m_\phi^5} n_\phi^2$$

**תנאי לדילול מספיק:** $\Gamma > H$ ב-$T_\phi \sim m_\phi$:

$$\mu_3 \gtrsim m_\phi \left(\frac{16\pi^3 H m_\phi^5}{n_\phi^2}\right)^{1/4}$$

עבור $m_\phi = 10$ MeV, $T \sim 10$ MeV: $\mu_3 \gtrsim O(1)$ MeV — **טבעי לחלוטין**.

**שלבי ביצוע:**
1. הוסף $\mu_3$ ו-$\lambda_4$ ללגרנז'יאן
2. כתוב Boltzmann equation עבור $n_\phi$ עם $3 \to 1$ rate
3. הראה ש-$\mu_3 \sim m_\phi$ מספיק לדלל $\phi$ לפני BBN
4. חשב $\Delta N_{\rm eff}$ שיורי (צריך להיות $\ll 0.3$)

### פתרון ב': Dark Sector קר ($\xi \ll 1$)

אם ה-heavy mediator $\Sigma$ מצמד חלש:

$$\xi \equiv \frac{T_\phi}{T_{\rm SM}} \ll 1$$

אז:
$$\Omega_\phi h^2 \propto \xi^3 \cdot \frac{m_\phi}{93 \text{ eV}}$$

עם $\xi \sim 0.01$: $\Omega_\phi h^2 \sim 10^4 \times 10^{-6} \sim 0.01$ — **מתחת ל-CDM budget**.

**אבל** — $\xi \ll 1$ משפיע גם על relic density ($\chi$ freeze-out ב-dark bath קר), ודורש פתרון מצומד.

### פתרון ג' (מועדף): Boltzmann מצומד

**הפתרון הנכון** שפותר גם §2 וגם §3 של הביקורת:

$$\frac{dn_\chi}{dt} + 3Hn_\chi = -\langle\sigma v\rangle_{\chi\chi\to\phi\phi}(n_\chi^2 - n_{\chi,\rm eq}^2)$$

$$\frac{dn_\phi}{dt} + 3Hn_\phi = +\langle\sigma v\rangle_{\chi\chi\to\phi\phi}(n_\chi^2 - n_{\chi,\rm eq}^2) - \Gamma_{3\to 1}(n_\phi^3 - n_\phi n_{\phi,\rm eq}^2)$$

$$\frac{d\rho_\phi}{dt} + 3H(\rho_\phi + P_\phi) = \text{energy injection from } \chi\chi \to \phi\phi$$

**שלבי ביצוע:**
1. כתוב solver מצומד (3 ODEs)
2. הוסף cannibal rate $3\phi \to \phi$
3. עקוב אחרי $T_\phi(t)$ ו-$n_\phi(t)$ בנפרד
4. אשר $\Omega_\phi h^2 \ll \Omega_\chi h^2$ ב-$T_{\rm BBN}$
5. חלץ $\Delta N_{\rm eff}$ מ-$\rho_\phi(T_{\rm BBN})$

---

## הלגרנז'יאן המעודכן

$$\boxed{\mathcal{L} = \frac{1}{2}\bar\chi(i\partial\!\!\!/ - m_\chi)\chi + \frac{1}{2}(\partial_\mu\phi)^2 - \frac{1}{2}m_\phi^2\phi^2 - \frac{1}{2}\bar\chi(y_s + iy_p\gamma_5)\chi\,\phi - \frac{\mu_3}{3!}\phi^3 - \frac{\lambda_4}{4!}\phi^4}$$

**פרמטרים חופשיים:** $m_\chi,\; m_\phi,\; y_s,\; y_p,\; \mu_3,\; \lambda_4$

| פרמטר | תפקיד | טווח טיפוסי |
|--------|--------|-------------|
| $m_\chi$ | מסת DM | 10–100 GeV |
| $m_\phi$ | מסת mediator | 5–50 MeV |
| $y_s$ ($\alpha_s$) | SIDM cross section | $\alpha_s \sim 10^{-3}$ |
| $y_p$ ($\alpha_p$) | Relic density (s-wave) | $\alpha_p \sim 10^{-3}$ |
| $\mu_3$ | Cannibal depletion of $\phi$ | $\mu_3 \sim m_\phi$ |
| $\lambda_4$ | $\phi$ self-interaction | $\lambda_4 \sim O(1)$ |

**מה לא השתנה:**
- ה-VPM solver (משתמש ב-$\alpha_s$ בלבד)
- ה-Higgs portal exclusion (חזק יותר עם $\alpha_{\rm tot} > \alpha_s$)
- כל ה-SIDM phenomenology (rotation curves, cluster offsets, gravothermal)

---

## סדר עבודה

### שלב 1: אישור אנליטי (חובה לפני כל קוד)
- [ ] חשב $|\mathcal{M}|^2(\chi\chi \to \phi\phi)$ עם $(y_s, y_p)$ — אשר s-wave $\propto y_p^4$
- [ ] חשב NR potential — אשר ש-$y_p$ contribution הוא $O(v^2)$
- [ ] חשב cannibal rate $\Gamma(3\phi \to \phi)$

### שלב 2: Relic density מעודכן
- [ ] עדכן v27_boltzmann_relic.py עם mixed coupling (או pure p-wave כ-fallback)
- [ ] כתוב coupled Boltzmann solver ($n_\chi + n_\phi + T_\phi$)
- [ ] אשר $\Omega_\chi h^2 = 0.120 \pm 0.001$
- [ ] אשר $\Omega_\phi h^2 < 0.001$ (ע"י cannibal depletion)

### שלב 3: סריקה מחדש
- [ ] סרוק ב-$(m_\chi, m_\phi, \alpha_s, \alpha_p)$ — 4D grid
- [ ] SIDM viability: $\sigma/m(30) \in [0.5, 10]$ cm²/g, $\sigma/m(1000) < 0.1$ cm²/g
- [ ] Relic: $\Omega h^2 = 0.12 \pm 0.012$
- [ ] חלץ island of viability חדש

### שלב 4: תיקוני באגים מכניים
- [ ] NFW $\rho_s$: הסר `/3.0` ב-3 קבצים
- [ ] CSV: תאם `m_phi_GeV` ↔ `m_phi_MeV`
- [ ] λ convention: בחר $\lambda = \alpha m_\chi/m_\phi$ (ללא factor 2) בכל מקום
- [ ] Benchmark labeling: BP1 = (20.69, 11.34, $\alpha_s$, $\alpha_p$) בכל המודולים

### שלב 5: עדכון Paper
- [ ] §2 (Lagrangian): עדכן ל-mixed coupling + $\mu_3\phi^3$
- [ ] §5.3 (Cosmology): החלף $\Delta N_{\rm eff}$ naive ← coupled Boltzmann result
- [ ] §6 (Relic density): s-wave עם $\alpha_p$, coupled equations
- [ ] הוסף appendix: derivation of $a_0(\chi\chi \to \phi\phi)$
- [ ] עדכן Table 1 עם $(\alpha_s, \alpha_p, \mu_3)$

---

## סיכון שיורי

| סיכון | הסתברות | השפעה | mitigation |
|-------|----------|-------|------------|
| $a_1$ (p-wave) קטן מהצפוי | 20% | $\alpha_p$ צריך להיות גדול, אולי מחוץ ל-SIDM range | mixed coupling פותר |
| Cannibal לא מספיק יעיל | 10% | $\phi$ overclosure נשאר | הגדל $\mu_3$ או הוסף $\phi \to \nu\nu$ (lepton portal) |
| $y_p$ contribution ל-VPM לא negligible | 15% | צריך לעדכן VPM solver | הוסף $y_p^2$ correction term |
| Island of viability מצטמצם | 30% | פחות benchmark points | מותר — גם 3 BPs מספיקים ל-paper |
| Naturalness של $y_s \sim y_p$ | 5% | Referee שואל "למה?" | RG running מהסקאלה הגבוהה |

---

## שורה תחתונה

> **המודל ניצל.** שינוי מינימלי בלגרנז'יאן — הוספת $y_p\gamma_5$ (pseudoscalar coupling) ו-$\mu_3\phi^3$ (cannibal term) — פותר את שתי הבעיות ה-FATAL בו-זמנית. ה-VPM solver, ה-SIDM phenomenology, וה-Higgs portal exclusion נשארים ללא שינוי. מה שנדרש הוא: (1) חישוב amplitude מפורש, (2) coupled Boltzmann solver, (3) סריקה מחדש ב-4D.
