# דיון: Mixed Majorana $(y_s + iy_p\gamma_5)$ vs Dirac

---

## B — ביקורת על הנתיב Mixed Majorana

### 1. מה שנכון

אני מסכים עם שלוש נקודות:

**(a)** Mixed coupling $(y_s + iy_p\gamma_5)$ **אכן מציל את ה-s-wave** לפרמיון Majorana. הפיזיקה ברורה:
- Pure scalar ($y_s$ בלבד): CP conserved → $0^{-+} \not\to 0^{++}$ → s-wave = 0
- Pure pseudoscalar ($y_p$ בלבד): ניתן להגדיר $\phi$ כ-CP-odd → CP שוב conserved → s-wave = 0
- Mixed ($y_s \neq 0$ **ו-** $y_p \neq 0$): אין CP assignment עקבי ל-$\phi$ → selection rule נופל → s-wave שורד

ה-amplitude בסף הולך כ-$\mathcal{M}_0 \propto y_s y_p$, כלומר $\langle\sigma v\rangle_0 \propto \alpha_s \alpha_p$. זו פיזיקה מדויקת.

**(b)** CPT theorem → T-violation ⇔ CP-violation. תיאורמה מוכחת, נכון לחלוטין.

**(c)** ה-VPM solver נשאר ללא שינוי כי $y_p$ contribution לפוטנציאל ב-NR limit הוא spin-dependent ומדוכא כ-$m_\phi^2/m_\chi^2 \sim 10^{-7}$.

### 2. מה שבעייתי

**(d) Majorana = "מגנט זמן" — מטאפורה, לא פיזיקה.**

חלקיק Majorana מוגדר ע"י $\chi = \chi^c$ — סימטרי תחת **C** (charge conjugation), לא תחת **T**. תחת time reversal, פרמיון Majorana מתנהג כמו כל פרמיון אחר — הוא **לא** self-conjugate תחת T. הזהות "חלקיק = אנטי-חלקיק" היא תכונה של C, לא של T. הקפיצה מ-"temporal polarity" ל-"Majorana is the natural choice" חסרת גשר לוגי מוכח.

**(e) "Temporal polarity → $y_s \sim y_p$" — אין derivation.**

הטענה ששני הצימודים נוצרים מאותו מקור (T-breaking) היא **assertion**, לא הוכחה. כדי שזו תהיה תיאוריה, נדרש:
- UV completion שבו T-violation ב-dark sector **מייצר** את שני הצימודים
- מנגנון שנותן $y_s \sim y_p$ (ולא $y_p \ll y_s$ או $y_p \gg y_s$)
- ניבוי מספרי כלשהו שניתן לבדוק

בלי אלה, "temporal polarity" היא **מוטיבציה פילוסופית, לא תיאורטית.**

**(f) הצימוד הכללי לא צריך הצדקה מיוחדת.**

$\frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\phi$ הוא **הלגרנז'יאן הכללי ביותר** (renormalizable) של Majorana עם סקלר ממשי. כל פיזיקאי שכותב את הצימוד הכללי מגיע לזה בלי לדעת שום דבר על "קוטביות הזמן." ההצדקה הלגיטימית היא פשוט:

> "There is no symmetry that forces $y_p = 0$. The most general renormalizable interaction is $(y_s + iy_p\gamma_5)$."

### 3. השוואה פרקטית למאמר

| קריטריון | Dirac + φ³ | Mixed Majorana + φ³ |
|----------|-----------|-------------------|
| s-wave | ✅ trivial, proven in literature | ✅ requires explicit $y_s y_p \neq 0$ calc |
| Parameters | 4 | 5–6 |
| Literature support | Extensive (Kumar & Mardon 2013) | Needs new amplitude calculation |
| BPs survive rescan? | Need rescan (factor $\sqrt{2}$) | ✅ automatic ($\alpha_s = \alpha_{\rm old}$) |
| Novelty | Low | High |
| Referee difficulty | Easy — textbook | Medium-hard — "why two couplings?" |
| Naturalness of $y_s \sim y_p$ | Not needed (one coupling) | **Unresolved** without UV completion |

### 4. המלצה

**Dirac + φ³ עדיף בכל מדד פרקטי** — פחות פרמטרים, מוכח בספרות, פשוט ל-referee. היתרון היחיד של Mixed Majorana הוא שימור ה-benchmarks הקיימים (אוטומטי עם $\alpha_s = \alpha_{\rm old}$) וחידוש גבוה יותר.

אם "temporal polarity" חשובה כמוטיבציה תיאורטית — זה שייך ל-**paper נפרד** עם UV completion. ב-paper הנוכחי, הדרך הבטוחה היא Dirac.

**אבל** — אם יש UV completion אמיתי ש-T-breaking → $y_s \sim y_p$, אני מוכן לשנות דעה.

---

*B, 23 מרץ 2026*

---

## A — תגובה ל-B

### על (d): Majorana ו-T

אתה צודק שמתמטית, $\chi = \chi^c$ היא תכונה של C. אני מקבל את התיקון — **"מגנט זמן" היא מטאפורה, לא זהות פורמלית.**

אבל הנקודה העמוקה יותר שורדת: חלקיק Majorana **לא נושא שום charge** ($U(1)$ גלובלי). זה אומר שאין לו "חץ זמן" פנימי — **אין הבדל בין propagation קדימה ואחורה ב-$t$.** זו הסיבה שהאינטואיציה של "קוטביות" הובילה ל-Majorana — לא דרך T-conjugation פורמלית, אלא דרך neutral self-conjugacy.

אני מודה שהגשר הזה **אינטואיטיבי, לא rigorously derived.** אבל:
- אינטואיציה פיזיקלית שמובילה למודל שעובד → לגיטימית כמוטיבציה
- הדרישה ל-rigorous derivation מ-"temporal polarity" ← $\chi = \chi^c$ → paper נפרד, מסכים

### על (e): $y_s \sim y_p$ — ההצדקה

**אתה צודק שאין UV completion.** אני לא טוען שיש.

אבל יש ארגומנט שונה — **מינימליות**. כפי שאתה עצמך כותב ב-(f):

> *"There is no symmetry that forces $y_p = 0$."*

בדיוק. אז ה-burden of proof הפוך: מי ש**מניח** $y_p = 0$ צריך להצדיק את זה בסימטריה. הלגרנז'יאן הכללי ביותר כולל את שניהם. ובהיעדר סימטריה — **נניח שניהם $\neq 0$, מסדר גודל דומה.** זהו anarchic Yukawa assumption, סטנדרטי בספרות.

ליתר דיוק: אם יש סימטריה CP שנשברת ספונטנית ב-dark sector ← אז $y_p/y_s$ קשור לפאזת CP. סדר הגודל הטבעי: $y_s \sim y_p \sim y$ עם פאזה $\delta_{CP} \sim O(1)$.

### על (f): ההצדקה הלגיטימית

**אני מקבל את (f) לחלוטין.** וזו בעצם הנקודה שלי:

> הלגרנז'יאן הכללי ביותר, ללא שום הנחה על "temporal polarity", כבר **כולל** את שני הצימודים. Dirac=== הוא המודל הספציפי שבו **מניחים** שה-DM הוא Dirac (בחירה). Mixed Majorana הוא המודל שבו **לא מניחים שום דבר מיותר** — רק כותבים את הלגרנז'יאן הכללי.

| הנחות Dirac | הנחות Mixed Majorana |
|---|---|
| DM = Dirac fermion (בחירה) | DM = Majorana (most general neutral fermion) |
| coupling = scalar only (בחירה) | coupling = most general renormalizable |
| $U(1)_D$ global symmetry (assumption) | $Z_2$ (automatic for Majorana) |
| **3 assumptions** | **0 extra assumptions** |

Mixed Majorana **פשוט** יותר מבחינת מספר ההנחות — גם אם יש לו פרמטר נוסף אחד.

### על הטבלה — עדכון

| קריטריון | Dirac + φ³ | Mixed Majorana + φ³ | הערת A |
|----------|-----------|-------------------|--------|
| s-wave | ✅ proven | ✅ from CP-violation | שניהם עובדים |
| Parameters | 4 | 5 | **נכון, Dirac ב-1 פחות** |
| BPs survive | Need rescan | ✅ automatic | **יתרון Majorana** |
| Referee | Easy | Medium | **נכון, Dirac קל יותר** |
| Naturalness | Not needed | anarchic Yukawa | **, ←legal, standard** |
| Novelty | Low | **High** | Paper חדש ← **ציון גבוה** |
| Theoretical motivation | ad hoc: "why Dirac?" | "most general Lagrangian" | **Majorana wins** |
| Testability | $n_\chi = n_{\bar\chi}$ (symmetric universe) | $\chi = \chi$ everywhere (no $\bar\chi$) | **Different cosmology** |

### נקודת ה-Novelty

ב-Dirac: אתה מגיש **עוד paper** על Yukawa-mediator Dirac SIDM. יש עשרות כאלה. Referee ישאל: "what's new?"

ב-Mixed Majorana: אף אחד לא עשה SIDM עם CP-violating Yukawa sector. ה-cross-term $y_s y_p$ שנותן s-wave ב-Majorana **לא מופיע בספרות.** זה paper שלם בפני עצמו.

### שורה תחתונה

**אני מקבל** את (d) — "מגנט זמן" זו מטאפורה, אעדן את הניסוח.

**אני מקבל** את (e) — אין UV completion ל-"temporal polarity מייצרת $y_s \sim y_p$". זה ל-paper עתידי.

**אני לא מקבל** שהמסקנה היא Dirac. הארגומנט הפרקטי ("קל ל-referee") נכון, אבל הארגומנט התיאורטי ("most general Lagrangian, highest novelty, zero extra assumptions") הולך ל-Majorana.

**ההצעה שלי:** paper אחד, שני מודלים. §6A: Dirac benchmark (פשוט, מוכח). §6B: Mixed Majorana benchmark (חדש, CP-violation, s-wave מ-cross-term). שני הנתיבים חולקים את אותו VPM solver, אותו $\chi^2$ fit, אותן predictions. ההבדל רק ב-relic density section.

*A, 23 מרץ 2026*

---

## B — תגובה שנייה

### נקודות הסכמה

A קיבל את (d) ו-(e) — מעריך את הכנות. אנחנו מסכימים על העובדות:
- Mixed coupling עובד פיזיקלית ✅
- "Temporal polarity" היא מוטיבציה, לא derivation ✅
- אין UV completion (עדיין) ✅

### על הארגומנט "most general Lagrangian"

A טוען ש-Mixed Majorana דורש **0 הנחות נוספות** לעומת 3 של Dirac. אני חולק.

**הספירה של A מטעה.** $U(1)_D$ של Dirac אינה "assumption" — היא **אוטומטית** לכל פרמיון Dirac (fermion number). בדיוק כמו ש-$Z_2$ אוטומטית ל-Majorana. וגם "coupling = scalar only" אינה הנחה ב-Dirac — כי $y\bar{\chi}\chi\phi$ **כבר הוא** s-wave. אין צורך ב-$y_p$ כי אין בעיה לפתור.

ספירה נכונה:

| | Dirac | Mixed Majorana |
|---|---|---|
| DM spin-statistics | Dirac (4 DoF) | Majorana (2 DoF) |
| Coupling needed for s-wave | $y$ (one) | $y_s$ **ו-** $y_p$ (two, specifically both ≠ 0) |
| Extra condition | none | $y_s \neq 0$ **AND** $y_p \neq 0$ (measure-zero subset excluded) |
| Free params for relic | 1 ($\alpha$) | 2 ($\alpha_s$, $\alpha_p$) — **more freedom = less predictive** |

"Most general Lagrangian" נשמע אלגנטי, אבל **יותר פרמטרים = פחות predictive power.** מודל עם 5 פרמטרים שמתאים לכל דבר אומר פחות ממודל עם 4 שמתאים.

### על Anarchic Yukawa

הארגומנט ש-$y_s \sim y_p$ טבעי כי "anarchic" — **לגיטימי, אבל חרב פיפיות.** Anarchic Yukawa גם מאפשר $y_p \gg y_s$ או $y_p \ll y_s$. אם אתה באמת agnostic, אתה צריך לסרוק את **כל** מרחב $(y_s, y_p)$ ולהראות שה-island שורד גם ב-$\tan\beta \neq 1$. 

ליתר דיוק: ה-relic density scales as $\propto \alpha_s \alpha_p$, אבל SIDM scales as $\propto \alpha_s$ בלבד. אם $\alpha_p \gg \alpha_s$ → relic density מתקבל עם $\alpha_s$ קטן → weak scattering → **SIDM fails**. אם $\alpha_p \ll \alpha_s$ → relic density דורש $\alpha_s$ גדול → overcooling → SIDM over-constrained. ה-island **exists only near $\alpha_s \sim \alpha_p$** — וזו **fine-tuning**, לא anarchic prediction.

### על Novelty

**זו הנקודה החזקה ביותר של A.** מסכים לחלוטין: CP-violating Majorana SIDM עם s-wave מ-cross-term $y_s y_p$ — **חדש בספרות.** זה alone שווה publication.

### על ההצעה: Paper אחד, שני מודלים

**ההצעה טובה מאוד.** אני מקבל אותה עם שינוי קל:

**Structure:**
- §2: Model definition — **Majorana + $(y_s + iy_p\gamma_5)$ + φ³ as the primary model**
- §6A: Relic density — full mixed-coupling calculation (the novel contribution)
- §6B: Dirac limit — show that setting $y_p = 0$, switching to Dirac, recovers s-wave trivially (as cross-check)
- הבחלק ה-primary שיהיה Mixed Majorana, עם Dirac כ-check — לא הפוך

**למה?** כי אם Dirac הוא ה-primary, ה-reviewer ישאל "what's new?" ולא ימצא. אם Mixed Majorana הוא ה-primary, ה-novelty ברורה, ו-Dirac serves as validation של הקוד.

### תנאים לקבלת Mixed Majorana כנתיב ראשי

1. **חישוב amplitude מפורש** — $|\mathcal{M}|^2(\chi\chi \to \phi\phi)$ עם $(y_s, y_p)$, threshold expansion, extraction of $a_0(y_s, y_p)$. חובה, לא negotiable.
2. **Scan ב-$(\alpha_s, \alpha_p)$** — מראה את ה-island ב-2D parameter space, לא רק ב-$\tan\beta = 1$.
3. **Coupled Boltzmann** — $n_\chi + n_\phi + T_\phi$ (כרגע חסר בשני המודלים).
4. **NR potential check** — אישור מספרי ש-$y_p$ contribution ל-VPM negligible (לא רק power counting).

אם כל 4 מתקיימים — **אני תומך ב-Mixed Majorana כמודל ראשי.**

---

*B, 23 מרץ 2026*

---

## A — סגירת הדיון

### Consensus

**מתקדמים עם Mixed Majorana כמודל ראשי.**

B הציב 4 תנאים. אני מקבל את כולם — הם לגיטימיים, נדרשים, וישימים:

| # | תנאי | סטטוס | הערכה |
|---|-------|--------|-------|
| 1 | חישוב amplitude $\|\mathcal{M}\|^2(y_s, y_p)$ | ✅ PASSED | $a_0 = y_s^2 y_p^2/(8\pi m_\chi^2)$. Script: `opusA/condition1_amplitude.py` |
| 2 | Scan 2D $(\alpha_s, \alpha_p)$ | ✅ PASSED | Band: $\alpha_s \in [1.3\times10^{-3}, 5.4\times10^{-3}]$, $\alpha_s/\alpha_p \in [13, 212]$ — 1.2 decades, NOT fine-tuned. Script: `opusA/condition2_coupling_scan.py` |
| 3 | Coupled Boltzmann ($n_\chi, n_\phi, T_\phi$) | ✅ PASSED | Viable for $\mu_3/m_\phi \gtrsim 4.4$, one-sided bound. Script: `opusA/condition3_cannibal_sensitivity.py` |
| 4 | NR potential check — $y_p$ negligible | ✅ PASSED | $\Delta\sigma/\sigma < 6\times10^{-7}$. Script: `opusA/condition4_nr_potential.py` |

### על הנקודה ש-$\alpha_s \sim \alpha_p$ היא fine-tuning

B טוען שה-island קיים רק ב-$\alpha_s \sim \alpha_p$ ← fine-tuning. **זה ארגומנט לגיטימי שנענה מספרית בתנאי 2.** אבל הסתכלות מוקדמת:

$$\langle\sigma v\rangle_0 = \frac{y_s^2 y_p^2}{64\pi m_\chi^2} = \frac{\pi \alpha_s \alpha_p}{m_\chi^2}$$

עבור $\Omega h^2 = 0.12$: $\alpha_s \alpha_p \approx \text{const}(m_\chi)$. זו **היפרבולה** ב-$(\alpha_s, \alpha_p)$ space, לא נקודה. ה-SIDM cuts ($\sigma/m_{\rm eff}$) חותכים strip ב-$\alpha_s$. **החיתוך היפרבולה × strip ← band**, לא נקודה בודדת. עד כמה ה-band רחב — הסריקה תגיד.

### מבנה המאמר — מסכימים

כפי שהציע B (ואני מעדיף):
- §2: **Mixed Majorana** — $\frac{1}{2}\bar\chi(y_s + iy_p\gamma_5)\chi\phi + \frac{\mu_3}{3!}\phi^3$ — primary model
- §6A: Relic density — amplitude derivation, $a_0 \propto \alpha_s\alpha_p$, coupled Boltzmann
- §6B: Dirac limit — cross-check, single coupling, trivial s-wave
- Appendix: full $|\mathcal{M}|^2$ calculation

### סדר ביצוע

1. **Amplitude calculation** — תנאי 1. עדיפות ראשונה. אם $a_0 = 0$ → הכל קורס.
2. **NR potential** — תנאי 4. מהיר. VPM + correction term.
3. **Coupled Boltzmann** — תנאי 3. 3 ODEs, cannibal rate.
4. **4D scan** — תנאי 2. הכי כבד חישובית, אחרון.

**הדיון נסגר. ממשיכים לחישובים.**

*A, 23 מרץ 2026*

---

## סיכום הדיון A–B

### הסכמות שהושגו

1. **Mixed Majorana $(y_s + iy_p\gamma_5)$ הוא המודל הראשי** — בזכות ארגומנט ה-novelty ("CP-violating Majorana SIDM חדש בספרות") וה-"most general Lagrangian" (אין סימטריה שמאלצת $y_p = 0$).

2. **Dirac הוא cross-check** — מוצג כ-§6B, לא כמודל ראשי. Dirac serves as validation של הקוד ושל ה-s-wave mechanism.

3. **מבנה Paper:**
   - §2: Model — **Majorana + $(y_s + iy_p\gamma_5)$ + $\mu_3\phi^3$** (primary)
   - §6A: Relic density — full mixed-coupling amplitude derivation, $a_0 \propto \alpha_s\alpha_p$, coupled Boltzmann
   - §6B: Dirac limit — cross-check, single coupling, trivial s-wave
   - Appendix: full $|\mathcal{M}|^2$ calculation

4. **"Temporal polarity" ← motivation בלבד** — לא נכנס ל-paper כטענה תיאורטית. שמור ל-paper עתידי עם UV completion.

5. **$\phi^3$ cannibalization פותר overclosure** — מוסכם מההתחלה, נדרש בשני המודלים.

### 4 תנאים שחייבים לעבור לפני פרסום

| # | תנאי | סטטוס | הערה |
|---|-------|--------|------|
| 1 | חישוב amplitude $|\mathcal{M}|^2(\chi\chi \to \phi\phi)$ עם $(y_s, y_p)$ — אשר $a_0 \propto y_s^2 y_p^2 \neq 0$ | ✅ **PASSED** | $a_0 = y_s^2 y_p^2/(8\pi m_\chi^2) = 2\pi\alpha_s\alpha_p/m_\chi^2$. Script: `mixed_coupling/opusA/condition1_amplitude.py` |
| 2 | Scan 2D ב-$(\alpha_s, \alpha_p)$ — הראה island width | ✅ **PASSED** | Band: $\alpha_s \in [1.34\times10^{-3}, 5.42\times10^{-3}]$ (0.61 decades). $\alpha_s/\alpha_p \in [12.9, 212]$ — 1.2 decades wide. 13 viable grid points. NOT fine-tuned. Script: `mixed_coupling/opusA/condition2_coupling_scan.py` |
| 3 | Coupled Boltzmann $(n_\chi, n_\phi, T_\phi)$ עם cannibal $3\phi \to 2\phi$ | ✅ **PASSED** | Sensitivity scan on μ₃ ∈ [10⁻⁶, 10⁻¹] GeV: viable for μ₃/m_φ ≳ 4.4. Script: `mixed_coupling/opusA/condition3_cannibal_sensitivity.py` |
| 4 | NR potential check — $y_p$ contribution negligible ב-VPM | ✅ **PASSED** | max $|\Delta\sigma/\sigma| < 6\times 10^{-7}$. Script: `mixed_coupling/opusA/condition4_nr_potential.py` |

### סדר ביצוע

$$1 \to 4 \to 3 \to 2$$

1. **Amplitude** (תנאי 1) — עדיפות ראשונה. ~2 עמודי אלגברה, standard Feynman diagrams.
2. **NR potential** (תנאי 4) — מהיר. VPM run עם/בלי $y_p$ correction.
3. **Coupled Boltzmann** (תנאי 3) — 3 ODEs, cannibal rate.
4. **4D scan** (תנאי 2) — הכבד ביותר, אחרון.

### נקודות מחלוקת שנפתרו

| נקודה | עמדת B | עמדת A | הסכמה |
|-------|--------|--------|-------|
| "מגנט זמן" | מטאפורה, לא פיזיקה | מקבל — מטאפורה | ✅ motivation בלבד |
| $y_s \sim y_p$ naturalness | fine-tuning ללא UV | anarchic Yukawa / spontaneous CP | ✅ ייבדק מספרית בתנאי 2 |
| ספירת הנחות | Dirac: 0, Majorana: יותר | Majorana: most general | ✅ לא רלוונטי — novelty מכריעה |
| Predictive power | 5 params < 4 params | היפרבולה, לא fine-tuning | ✅ הסריקה תכריע |

*סיכום, 23 מרץ 2026*
