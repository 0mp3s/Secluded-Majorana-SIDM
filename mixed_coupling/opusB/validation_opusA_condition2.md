# ולידציה: Opus A — Condition 2 (2D Coupling Scan)

**בודק:** Opus B  
**תאריך:** 23 מרץ 2026  
**קבצים נבדקים:**
- `mixed_coupling/opusA/condition2_coupling_scan.py`
- `mixed_coupling/opusA/condition2_output.txt`

---

## תוצאת ולידציה

### ✅ קוד ותוצאות — נכונים. אין שגיאות קריטיות.

### הערות MINOR בלבד (ראה סעיפים 6–7).

---

## 1. נוסחת Relic Density

```python
sv0 = 2.0 * math.pi * alpha_s * alpha_p / M_CHI**2
```

$$\langle\sigma v\rangle_0 = \frac{2\pi\alpha_s\alpha_p}{m_\chi^2}$$

**נכון.** תואם את נוסחת V11 שאומתה ב-Condition 1. ✅

### אימות מספרי של המכפלה:

$$\alpha_s \times \alpha_p = 1.3875 \times 10^{-7}$$

$$\langle\sigma v\rangle_0 = \frac{2\pi \times 1.3875 \times 10^{-7}}{(20.69)^2} = 2.04 \times 10^{-9}\ \text{GeV}^{-2}$$

המרה ליחידות קונבנציונליות:

$$\langle\sigma v\rangle = 2.04 \times 10^{-9} \times \underbrace{3.89 \times 10^{-28}}_{\text{GeV}^{-2}\to\text{cm}^2} \times \underbrace{3.0 \times 10^{10}}_c = 2.4 \times 10^{-26}\ \text{cm}^3/\text{s}$$

ערך קנוני נכון ל-Majorana WIMP עם $m \sim 20$ GeV. ✅

---

## 2. Boltzmann Solver

RK4 סטנדרטי, אותו solver כמו ב-Condition 3 Phase 1 (שכבר אומת).

| פריט | סטטוס |
|------|--------|
| `x_start=1`, `x_end=500`, `n_steps=15000` | ✅ |
| $\lambda = \sqrt{\pi/45}\,g_\text{eff}\,M_\text{Pl}\,m_\chi$ | ✅ |
| $dY/dx = -\lambda\,\langle\sigma v\rangle_0\,(Y^2 - Y_\text{eq}^2)/x^2$ | ✅ |
| $Y \to \Omega h^2$ via $m Y s_0 / \rho_{c,h^2}$ | ✅ |

### Bisection לנקודה סימטרית:

A מבצע 60 חצויות (relative accuracy $\sim 2^{-60} \sim 10^{-18}$) — יותר ממספיק. ✅

תוצאה: $\alpha_\text{sym} = 3.7249 \times 10^{-4}$ עם $\Omega h^2 = 0.120000$.

### אימות קונסיסטנטיות:

A בודק 7 נקודות אסימטריות על ההיפרבולה $\alpha_s \times \alpha_p = \text{const}$. כולן נותנות $\Omega h^2 = 0.120000$ — מאשר שהקונטור הוא היפרבולה מושלמת. ✅

(צפוי: עבור s-wave בלבד, $\langle\sigma v\rangle_0 = 2\pi\alpha_s\alpha_p/m^2$, ה-yield $Y_\infty \propto 1/\langle\sigma v\rangle_0$, כך ש-$\Omega h^2$ תלוי רק במכפלה $\alpha_s\alpha_p$.)

---

## 3. VPM: שימוש ב-$\alpha_s$ בלבד

```python
sm = sigma_T_vpm(M_CHI, M_PHI, alpha_s, v)
```

`sigma_T_vpm` מחשב `lam = alpha * m_chi / m_phi` ומשתמש ב-Yukawa potential $V(r) = -\alpha\,e^{-m_\phi r}/r$.

**עבור V11:**
- צימוד סקלרי $y_s$: פוטנציאל NR $V_s = -\alpha_s\,e^{-m_\phi r}/r$ (spin-independent, leading order)
- צימוד פסאודו-סקלרי $y_p$: פוטנציאל NR $V_p \sim \alpha_p\,(\vec\sigma_1\cdot\vec\sigma_2)\,\nabla^2(e^{-m_\phi r}/r)/(4m_\chi^2)$ — מדוכא ב-$(m_\phi/m_\chi)^2 \sim 3 \times 10^{-7}$

**מאומת ב-Condition 4** (nr_potential_check.py, $\delta\sigma/\sigma < 10^{-6}$).

העברת $\alpha_s$ ל-VPM — **נכון**. ✅

---

## 4. גיאומטריית Band

הטענה המרכזית:

$$\text{Relic:} \quad \alpha_s \times \alpha_p = C \quad (\text{היפרבולה})$$
$$\text{SIDM:} \quad \alpha_s \in [\alpha_s^{\min}, \alpha_s^{\max}] \quad (\text{רצועה אנכית})$$

$$\text{חיתוך} = \text{BAND לאורך ההיפרבולה}$$

**נכון גיאומטרית.** זו אחת מההשלכות המרכזיות של V11 — ב-V10 (צימוד יחיד $\alpha$), הודות ל-$\langle\sigma v\rangle \propto \alpha^2$, VPM תלוי ב-$\alpha$ ו-relic תלוי ב-$\alpha^2$, כך שמתקבלת **נקודה** (לא רצועה). ב-V11, שברו של החופש האקסטרה ($\alpha_p$ עצמאי) נותן **רצועה**. ✅

---

## 5. תוצאות מספריות

| כמות | ערך | אימות |
|------|-----|-------|
| $\alpha_s \times \alpha_p$ | $1.387 \times 10^{-7}$ | ✅ (ראה §1) |
| SIDM strip: $\alpha_s \in$ | $[1.34\times10^{-3}, 5.42\times10^{-3}]$ | ✅ (ראה §5.1) |
| Band width | 0.61 decades | $= \log_{10}(5.42/1.34) = 0.607$ ✅ |
| $\alpha_s/\alpha_p$ range | $[13, 212]$ | $= [1.34\times10^{-3}/1.04\times10^{-4}, 5.42\times10^{-3}/2.56\times10^{-5}]$ ✅ |
| Viable points | 13 out of 80 | Consistent with strip width |

### 5.1 Spot-check: SIDM boundaries

Lower edge ($\alpha_s = 1.34 \times 10^{-3}$):
- $\sigma/m(1000) = 0.113$ cm²/g — barely above 0.1 cut → marginal, makes sense as boundary ✅

Upper edge ($\alpha_s = 5.42 \times 10^{-3}$):
- $\sigma/m(1000) = 0.978$ cm²/g — just below 1.0 cut → marginal, makes sense as boundary ✅

Point near center ($\alpha_s = 2.69 \times 10^{-3}$):
- $\sigma/m(30) = 1.37$, $\sigma/m(100) = 1.17$, $\sigma/m(1000) = 0.37$ — typical velocity-dependent SIDM ✅

### 5.2 Cross-check: VPM at BP1-like coupling

$\alpha_s = 1.06 \times 10^{-3}$ (close to original BP1 $\alpha = 1.048 \times 10^{-3}$) gives $\sigma/m(30) = 0.52$ cm²/g.

This matches the known V10 BP1 self-interaction cross-section, confirming VPM consistency with old results. ✅

---

## 6. SIDM Cuts — הערות

| Velocity | Lower bound | Upper bound | Standard? |
|----------|------------|------------|-----------|
| 30 km/s (dwarf) | 0.1 cm²/g | **50 cm²/g** | Upper bound generous |
| 100 km/s (LSB) | 0.1 cm²/g | 10 cm²/g | ✅ |
| 1000 km/s (cluster) | 0.1 cm²/g | 1.0 cm²/g | ✅ |

**הערה 1:** הגבול העליון ב-dwarfs (50 cm²/g) גבוה מהספרות הסטנדרטית (~10 cm²/g, e.g., Kaplinghat+ 2016). עם גבול עליון מחמיר של 10 cm²/g, הקצה העליון של הרצועה ירד — אבל רוב הנקודות הכדאיות (σ/m(30) ∈ [0.7, 2.4]) שורדות, כך שהבנד לא נעלם.

**הערה 2:** הגבולות **התחתונים** (σ/m > 0.1) אינם אילוצים תצפיתיים — הם דרישה מודלית ("אנחנו רוצים SIDM שעושה משהו"). זה לגיטימי עבור מודל SIDM, אבל חשוב להבדיל בין observational upper bound לבין model requirement lower bound.

**השפעה:** MINOR — לא משנה את המסקנה העקרונית.

---

## 7. הערה: ההיררכיה $\alpha_s \gg \alpha_p$

כל ה-band דורש $\alpha_s/\alpha_p \in [13, 212]$, כלומר:

$$y_s/y_p = \sqrt{\alpha_s/\alpha_p} \in [3.6, 14.6]$$

זו היררכיה של סדר גודל אחד ב-Yukawa couplings. פיזיקלית, הייררכיה Naturality תלויה בסימטריה של ה-UV completion:

- אם יש CP symmetry (אפילו approximate), $y_p = 0$ טבעי → כל $y_p$ קטן הוא O.K.
- בלי סימטריה, $y_s \sim y_p$ צפוי → ההיררכיה דורשת הסבר.

A אומר "NOT fine-tuned" כי הבנד רחב 1+ decade — וזה נכון (זו לא **נקודה**). אבל עצם הצורך ב-$y_s \gg y_p$ ראוי לדיון.

**סטטוס:** העיר COSMETIC — לא באג, אלא נקודה שכדאי להזכיר בדיון.

---

## 8. עקביות עם ולידציות קודמות

| בדיקה | תוצאה |
|-------|-------|
| SV0 formula matches Condition 1 | ✅ |
| VPM convention matches Condition 4 | ✅ |
| Cluster velocity = 1000 km/s (matches Harvey+15 fix) | ✅ |
| Independent of $\mu_3$ (no cannibal here) | ✅ (consistent with Condition 3) |

---

## 9. סיכום

| קטגוריה | ממצא | חומרה |
|---------|-------|--------|
| **SV0 formula** | $2\pi\alpha_s\alpha_p/m^2$ ✅ | — |
| **Boltzmann solver** | RK4 סטנדרטי, 15K steps ✅ | — |
| **Relic product** | $\alpha_s\alpha_p = 1.387\times10^{-7}$ ✅ | — |
| **VPM with $\alpha_s$** | נכון (Condition 4 מאמת) ✅ | — |
| **Band geometry** | Hyperbola ∩ strip = band ✅ | — |
| **Band width** | 0.61 decades ✅ | — |
| **SIDM dwarf upper** | 50 cm²/g — generous | MINOR |
| **SIDM lower bounds** | Model choice, not observational | MINOR |
| **$\alpha_s \gg \alpha_p$ hierarchy** | Requires discussion | COSMETIC |

### ציון: ✅ PASS

הקוד נכון, התוצאות עקביות, המסקנה (band not point → no fine-tuning) מבוססת.

שתי ההערות הקלות:
1. שקלו להוריד את הגבול העליון ב-dwarfs ל-10 cm²/g (לא ישנה את הבנד באופן מהותי)
2. הוסיפו דיון על הצורך ב-$y_s \gg y_p$ ואם יש סימטריית CP approximate שמצדיקה זאת
