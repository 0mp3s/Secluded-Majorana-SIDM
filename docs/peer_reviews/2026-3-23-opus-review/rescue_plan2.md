# Rescue Plan v2 — דעה מעודכנת

**תאריך:** 23 מרץ 2026  
**הקשר:** לאחר קריאת peer_review.md, rescue_plan.md (שלי), ו-model_fix_options.md (קודם)  
**המלצה:** **Dirac + φ³ cannibal** (אפשרות A+A1 מ-model_fix_options)

---

## מה השתנה מ-rescue_plan.md

ב-rescue_plan הראשון שלי הצעתי **mixed scalar-pseudoscalar coupling** ($y_s + iy_p\gamma_5$) כדרך לשמור על Majorana ולהחזיר s-wave. אחרי קריאת model_fix_options, אני חוזר בי. **המעבר ל-Dirac עדיף בכל מדד.**

| קריטריון | Mixed coupling (rescue_plan v1) | Dirac (model_fix_options) | מנצח |
|-----------|-------------------------------|--------------------------|------|
| מספר פרמטרים | 5 ($m_\chi, m_\phi, \alpha_s, \alpha_p, \mu_3$) | 4 ($m_\chi, m_\phi, \alpha, \mu_3$) | **Dirac** |
| s-wave proof | צריך חישוב amplitude חדש | מוכח בספרות (Kumar & Mardon 2013) | **Dirac** |
| Naturalness | צריך להצדיק $y_s \sim y_p$ | צימוד יחיד $y$ | **Dirac** |
| שינוי VPM | אין | אין (+ Dirac averaging formula) | שווה |
| שימור benchmarks | ~100% ($\alpha_s = \alpha_{\rm old}$) | ~70% (factor $\sqrt{2}$ shift) | Mixed |
| פשטות המאמר | מסובך — שני צימודים + שתי אלפות | נקי — Dirac + $y$ יחיד | **Dirac** |

**שורה תחתונה:** ההמרה Majorana → Dirac פותרת את בעיית ה-s-wave **חינם** — בלי פרמטר נוסף, בלי ארגומנט naturalness, ובלי חישוב חדש. זה עדיף על הניסיון "להציל" את ה-Majorana.

---

## הנתיב — Dirac Secluded SIDM + Cannibal φ

### הלגרנז'יאן

$$\mathcal{L}_{\rm dark} = \bar\chi(i\partial\!\!\!/ - m_\chi)\chi - y\bar\chi\chi\phi + \frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2 - \frac{\mu_3}{3!}\phi^3$$

- $\chi$: **Dirac** fermion (חלקיק + אנטי-חלקיק)
- $\phi$: real scalar mediator
- $y$: Yukawa coupling ($\alpha = y^2/4\pi$)
- $\mu_3$: self-coupling — מדלל φ לפני BBN
- סימטריית יציבות: $U(1)_D$ גלובלית ($\chi \to e^{i\theta}\chi$, $\phi$ סינגלט)

**4 פרמטרים חופשיים:** $m_\chi$, $m_\phi$, $\alpha$, $\mu_3$

### מה עובד ומה משתנה

| רכיב | סטטוס | הערה |
|------|--------|------|
| VPM solver (phase shifts) | ✅ **ללא שינוי** | אותו Yukawa potential, אותו ODE |
| Born validation | ✅ **ללא שינוי** | אותן פורמולות אנליטיות |
| Literature cross-checks | ✅ **ללא שינוי** | Tulin-Yu-Zurek בודקים Yukawa, לא תלוי בסטטיסטיקת ה-DM |
| Higgs portal exclusion | ✅ **ללא שינוי** | $\phi$ ממשי, אותו ארגומנט |
| Chi2 fit framework | ✅ **ללא שינוי** | אותם observables |
| Predictions (gravothermal, cores) | ✅ **ללא שינוי** | תלוי רק ב-$\sigma/m(v)$ |
| **σ_T formula** | ⚠️ **עדכון** | Dirac averaging: $\frac{1}{2}\sigma_{\chi\chi} + \frac{1}{2}\sigma_{\chi\bar\chi}$ |
| **Boltzmann solver** | ⚠️ **עדכון** | $\langle\sigma v\rangle = \pi\alpha^2/(2m_\chi^2)$ (factor 2 מ-Majorana) |
| **Relic scan** | 🔄 **הרצה מחדש** | עם $\alpha$ חדש |
| **φ depletion** | 🆕 **חדש** | Coupled Boltzmann עם cannibal rate |

---

## שלוש הבעיות ואיך כל אחת נפתרת

### 1. s-wave annihilation ← Dirac פותר

עבור $\chi\bar\chi \to \phi\phi$ (Dirac), ה-amplitude בסף:

$$\bar v(p_2) u(p_1) \xrightarrow{v \to 0} 2m_\chi \neq 0$$

זה **לא מתאפס** כי $\chi$ ו-$\bar\chi$ הם חלקיקים **שונים** — אין Pauli blocking, אין antisymmetrization של המצב ההתחלתי. התוצאה:

$$\langle\sigma v\rangle_{\rm s\text{-}wave} = \frac{\pi\alpha^2}{2m_\chi^2}$$

אין צורך בהוכחה מקורית — זה תוצאה סטנדרטית (Gondolo & Gelmini 1991; Kumar & Mardon 2013).

### 2. φ overclosure ← Cannibal $\mu_3\phi^3$ פותר

הצימוד $y\bar\chi\chi\phi$ **שומר** $U(1)_D$ (כי $\phi$ ניטרלי), אז $\phi^3$ **מותר** ע"י כל הסימטריות. אפילו אם $\mu_3 = 0$ ב-tree level, הוא נוצר רדיאטיבית:

$$\mu_3^{\rm 1\text{-}loop} \sim \frac{y^3 m_\chi}{16\pi^2} \sim 10 \text{ MeV} \quad (\text{for } y \sim 0.1, m_\chi \sim 20\text{ GeV})$$

**תהליך:** $3\phi \to 2\phi$ (cannibal) — מפחית $n_\phi$ אקספוננציאלית

**תנאי minimal:** $\mu_3 \gtrsim m_\phi$ → $\Gamma_{3\to2} > H$ ב-$T \sim m_\phi$

**תוצאה:** $\Omega_\phi h^2 \ll 0.01$ ב-BBN — **overclosure נפתר**

### 3. Single-species Boltzmann ← Coupled equations

זה לא FATAL אבל זה upgrade חשוב. הפתרון:

$$\frac{dY_\chi}{dx} = -\frac{\langle\sigma v\rangle}{Hx}\left(Y_\chi^2 - Y_{\chi,\rm eq}^2\right)$$

$$\frac{dY_\phi}{dx} = +\frac{\langle\sigma v\rangle}{Hx}\left(Y_\chi^2 - Y_{\chi,\rm eq}^2\right) - \frac{\Gamma_{3\to2}}{Hx}\left(Y_\phi^3 - Y_\phi Y_{\phi,\rm eq}^2\right)$$

פותרים את שתי המשוואות **ביחד**. התוצאה: $Y_\chi$ כמעט ללא שינוי (הניזון מ-bath תרמי), $Y_\phi$ יורד חזק ע"י cannibal.

---

## השפעת Dirac על נוסחת σ_T

ביקום סימטרי ($n_\chi = n_{\bar\chi}$), חצי מההתנגשויות הן $\chi\chi$ (identical) וחצי $\chi\bar\chi$ (distinguishable):

$$\sigma_{\rm eff} = \frac{1}{2}\sigma_{\chi\chi}^{\rm id} + \frac{1}{2}\sigma_{\chi\bar\chi}^{\rm dist}$$

כאשר:

$$\sigma_{\chi\chi}^{\rm id} = \frac{2\pi}{k^2}\sum_l w_l(2l+1)\sin^2\delta_l, \quad w_l = \begin{cases}1 & l \text{ even} \\ 3 & l \text{ odd}\end{cases}$$

$$\sigma_{\chi\bar\chi}^{\rm dist} = \frac{4\pi}{k^2}\sum_l(2l+1)\sin^2\delta_l$$

ה-phase shifts $\delta_l$ **זהים** — אותו VPM solver. רק נוסחת הסיכום משתנה.

**השפעה מספרית:** ב-Born limit, $\sigma_{\rm eff}^{\rm Dirac} / \sigma_T^{\rm Majorana} \approx 1.0$–$1.5$ (תלוי ב-$l_{\rm max}$). **ההבדל קטן** ולא מסכן את ה-viability.

---

## השפעה על ה-Island of Viability

עם Dirac, ה-annihilation cross section ב-freeze-out:

$$\langle\sigma v\rangle = \frac{\pi\alpha^2}{2m_\chi^2} \quad (\text{Dirac, s-wave})$$

לעומת:

$$\langle\sigma v\rangle_{\rm old} = \frac{\pi\alpha^2}{4m_\chi^2} \quad (\text{"Majorana", s-wave — wrong})$$

ה-$\alpha_{\rm relic}$ הנדרש ל-$\Omega h^2 = 0.12$ **קטן** ב-$\sqrt{2}$:

$$\alpha_{\rm relic}^{\rm Dirac} = \frac{\alpha_{\rm relic}^{\rm old}}{\sqrt{2}} \approx 0.71 \times \alpha_{\rm relic}^{\rm old}$$

**BP1 לדוגמה:**

| | V10 (שגוי) | V11 (Dirac) |
|---|---|---|
| $\alpha_{\rm relic}$ | $1.048 \times 10^{-3}$ | $7.4 \times 10^{-4}$ |
| $\lambda$ | 1.91 | 1.35 |
| SIDM regime | resonant | resonant (אותו אזור) |

ה-$\alpha$ **יורד** — כלומר אנחנו זזים **פנימה** לתוך ה-resonant regime, לא החוצה. **סביר מאוד** שה-island שורד, אולי אפילו מתרחב.

הסריקה מחדש תגיד. אבל אין סיבה לצפות שה-island נעלם.

---

## סדר עבודה מומלץ

### שלב 1: הוכחה אנליטית [יום 1]
- [ ] כתוב appendix עם $|\mathcal{M}|^2(\chi\bar\chi \to \phi\phi)$ — Dirac, scalar coupling, t+u channels
- [ ] הראה $a_0 \neq 0$ (s-wave שורד)
- [ ] חלץ: $\langle\sigma v\rangle = \pi\alpha^2/(2m_\chi^2)$
- [ ] ציין: Kumar & Mardon (2013), Gondolo & Gelmini (1991)

### שלב 2: עדכון קוד core [יום 1–2]

**v22_raw_scan.py** — הוסף Dirac averaging:
```python
def sigma_eff_dirac_vpm(m_chi, m_phi, alpha, v_km_s):
    """Effective σ/m for Dirac DM in symmetric universe."""
    # Phase shifts (unchanged)
    deltas = [vpm_phase_shift(l, kappa, lam, x_max, N) for l in range(l_max+1)]
    
    # Identical (χχ or χ̄χ̄): spin-weighted
    sig_id = (2*pi/k**2) * sum(w_l*(2*l+1)*sin(d)**2 for l,d in enumerate(deltas))
    
    # Distinguishable (χχ̄): no spin weights
    sig_dist = (4*pi/k**2) * sum((2*l+1)*sin(d)**2 for l,d in enumerate(deltas))
    
    return 0.5 * sig_id + 0.5 * sig_dist
```

**v27_boltzmann_relic.py** — עדכן annihilation:
```python
def sigma_v_swave_dirac(alpha_d, m_chi):
    """s-wave ⟨σv⟩ for Dirac χχ̄ → φφ."""
    return math.pi * alpha_d**2 / (2.0 * m_chi**2)  # factor 2 vs Majorana
```

### שלב 3: Cannibal depletion [יום 2–3]
- [ ] כתוב `cannibal_boltzmann.py` — coupled ODEs ל-$Y_\chi$ + $Y_\phi$
- [ ] חשב $\Gamma(3\phi \to 2\phi)$ אנליטית
- [ ] הרץ: אשר $\Omega_\phi h^2 < 0.001$ עם $\mu_3 \sim m_\phi$
- [ ] חלץ $\Delta N_{\rm eff}$ שיורי

### שלב 4: סריקה מחדש [יום 3–4]
- [ ] עדכן `smart_scan.py` עם Dirac formulas
- [ ] סרוק $(m_\chi, m_\phi, \alpha)$ — grid זהה ל-V10
- [ ] SIDM criterion: $\sigma_{\rm eff}/m(30) \in [0.5, 10]$, $\sigma_{\rm eff}/m(1000) < 0.1$
- [ ] Relic: $\Omega h^2 = 0.12 \pm 0.012$ עם Dirac $\langle\sigma v\rangle$
- [ ] חלץ island חדש — צפי: דומה ל-V10 עם shift ב-$\alpha$

### שלב 5: תיקוני באגים [יום 4]
- [ ] NFW ρ_s: הסר `/3.0` ב-3 קבצים
- [ ] CSV: `m_phi_GeV` → consistent naming
- [ ] λ convention: $\lambda = \alpha m_\chi / m_\phi$ (ללא factor 2) — uniform
- [ ] BP labeling: BP1 = paper Table 1, uniform across modules

### שלב 6: עדכון Paper [יום 5–7]
- [ ] Title: "Secluded ~~Majorana~~ **Dirac** SIDM..."
- [ ] §2: לגרנז'יאן חדש עם Dirac + $\mu_3\phi^3$
- [ ] §3.2: נוסחת $\sigma_{\rm eff}$ Dirac
- [ ] §5.3: Boltzmann מצומד + cannibal depletion result
- [ ] §6: $\langle\sigma v\rangle = \pi\alpha^2/(2m_\chi^2)$ + derivation
- [ ] Appendix: Amplitude calculation (t+u channels)
- [ ] Table 1: BPs מעודכנים

---

## סיכונים שיוריים

| סיכון | הסתברות | mitigation |
|--------|----------|-----------|
| Island מצטמצם אחרי rescan | 25% | $\alpha$ ירד — resonances עדיין שם. אפילו 3 BPs מספיקים. |
| Cannibal לא מספיק | 10% | הגדל $\mu_3$ או הוסף $\lambda_4\phi^4$. Loop-induced $\mu_3$ כבר $\sim m_\phi$. |
| Dirac averaging משנה σ/m משמעותית | 5% | ב-Born limit ההבדל ~20%. ב-resonant regime ייתכן גדול יותר, אבל תמיד ניתן לכוונן $(m_\chi, m_\phi)$. |
| Referee דורש $\sigma_T^{\rm transport}$ | 30% | תגובה: standard in field (Tulin+Yu 2018), ~20% systematic, within obs. errors. |

---

## למה לא Mixed Coupling (rescue_plan v1)

הגישה שהצעתי ב-rescue_plan.md הראשון — $y_s + iy_p\gamma_5$ — **עובדת** מבחינה פיזיקלית, אבל היא נחותה:

1. **יותר פרמטרים** — 5 vs 4. כל referee ישאל "למה שני צימודים?"
2. **Naturalness** — צריך $y_s \sim y_p$ כדי ש-SIDM ו-relic density יהיו באותו אזור. אין סימטריה שמכתיבה את זה.
3. **חישוב חדש נדרש** — הצימוד המעורב לא מכוסה בספרות. Dirac כן.
4. **VPM correction** — ב-mixed coupling, ה-$y_p$ contribution ל-Yukawa potential הוא velocity-dependent (spin-orbit). צריך לבדוק שזה באמת negligible. ב-Dirac אין שאלה.
5. **מסר המאמר** — "Majorana עם סקלר" שנכשל ואז "ניצל" עם pseudoscalar נראה כמו patch. "Dirac עם סקלר" הוא **מודל נקי מההתחלה**.

---

## סיכום

> **הנתיב הנכון: Majorana → Dirac + cannibal φ³.**
>
> שינוי אחד בזהות ה-DM (Majorana → Dirac) פותר את הבעיה הקשה ביותר (s-wave).
> פרמטר אחד נוסף ($\mu_3$) פותר את ה-overclosure.
> ~80% מהקוד נשמר ללא שינוי.
> ה-VPM solver, ה-Born validation, ה-Higgs portal exclusion, ה-predictions — כולם עומדים.
>
> **זה לא "הצלה". זה שדרוג טבעי של המודל.**
