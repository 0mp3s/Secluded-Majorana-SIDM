# ולידציה: Opus A — Condition 3 (Cannibal φ³ Sensitivity)

**בודק:** Opus B  
**תאריך:** 23 מרץ 2026  
**קבצים נבדקים:**
- `mixed_coupling/opusA/condition3_cannibal_sensitivity.py`
- `mixed_coupling/opusA/condition3_output.txt`

---

## תוצאת ולידציה

### ❌ CRITICAL: נוסחת הקרוס-סקשן $\langle\sigma v^2\rangle_{3\to2}$ שגויה ממדית ומספרית

### ✅ מבנה הקוד, Phase 1, והמסקנות האיכותיות — נכונים

---

## 1. שגיאה קריטית: נוסחת הקניבל

### מה A כתב (שורה 194 + Docstring שורה 37):

$$\langle\sigma v^2\rangle_{3\to2} = \frac{25\sqrt{5}}{512\pi}\frac{\mu_3^6}{m_\phi^9}$$

### מה נכון (Farina, Pappadopulo, Ruderman, Trevisan 2016, eq. 7):

$$\langle\sigma v^2\rangle_{3\to2} = \frac{25\sqrt{5}}{512\pi(4\pi)^6}\frac{\mu_3^6}{m_\phi^{11}}$$

### שתי שגיאות:

| פריט | A | נכון | סיבה |
|------|---|------|-------|
| חזקה של $m_\phi$ | $m_\phi^{-9}$ | $m_\phi^{-11}$ | חסרים 2 propagators × $1/m^2$ |
| מקדם מספרי | $\frac{25\sqrt{5}}{512\pi}$ | $\frac{25\sqrt{5}}{512\pi(4\pi)^6}$ | חסר $(4\pi)^{-6}$ מ-phase space |

### אימות ממדים:

ממשוואת Boltzmann: $\dot{n} + 3Hn = -\langle\sigma v^2\rangle(n^3 - n^2 n_\text{eq})$

$$[\dot{n}] = E^4, \quad [n^3] = E^9 \implies [\langle\sigma v^2\rangle] = E^{-5}$$

| נוסחה | ממד | תקין? |
|-------|-----|-------|
| A: $\mu_3^6/m^9$ | $E^6/E^9 = E^{-3}$ | ❌ |
| נכון: $\mu_3^6/m^{11}$ | $E^6/E^{11} = E^{-5}$ | ✅ |

### גודל השגיאה:

$$\frac{\text{A's rate}}{\text{correct rate}} = (4\pi)^6 \times m_\phi^2 = 3.94 \times 10^6 \times 1.29 \times 10^{-4} \approx 508$$

A מפריז את קצב הקניבל פי ~500. זה אומר שהקניבל פחות יעיל ממה ש-A חישב.

### השפעה על חלון הכדאיות:

$\langle\sigma v^2\rangle \propto \mu_3^6$ → כדי לפצות על הפקטור 508, צריך:
$$\mu_3^{(\text{correct})} = \mu_3^{(\text{A})} \times 508^{1/6} \approx \mu_3^{(\text{A})} \times 2.83$$

A מצא $\mu_3/m_\phi \gtrsim 4.4$ → **Corrected: $\mu_3/m_\phi \gtrsim 12$**

עדיין $O(10)$ — המסקנה האיכותית ($\mu_3 \gtrsim$ few × $m_\phi$) שורדת.

### מקור מלא:

Farina, Pappadopulo, Ruderman, Trevisan, "Phases of Cannibal Dark Matter" (2016):

$$V = \frac{m\lambda_3}{3!}\phi^3 \implies \mu_3 = m\lambda_3$$

$$\langle\sigma_{3\to2} v^2\rangle = \frac{25\sqrt{5}}{512\pi m^5}\left(\frac{\lambda_3}{4\pi}\right)^6 = \frac{25\sqrt{5}}{512\pi(4\pi)^6}\frac{\mu_3^6}{m^{11}}$$

### אימות ממדי עצמאי מה-amplitude:

עבור $\frac{\mu_3}{3!}\phi^3$ ברמת עץ:
- $\mathcal{M}_{3\to2}$ דורש 3 vertices ($\mu_3^3$) ו-2 propagators ($\sim 1/m_\phi^4$ ב-NR):
  $$[\mathcal{M}] = E^3/E^4 = E^{-1} \quad \checkmark \text{ (dim = 4-3-2=-1 for 3→2)}$$
- $\sigma v^2 \sim |\mathcal{M}|^2/m^3 \sim \mu_3^6/(m^8 \cdot m^3) = \mu_3^6/m^{11}$ → $E^{-5}$ ✓

---

## 2. Phase 1: χ Freeze-out

| פריט | סטטוס | הערה |
|------|--------|------|
| RK4 solver | ✅ | Standard, 20000 steps |
| SV0 = $2\pi\alpha_s\alpha_p/m_\chi^2$ | ✅ | נוסחה נכונה (V11) |
| $\lambda = \sqrt{\pi/45}\,g_\text{eff}\,M_\text{Pl}\,m_\chi$ | ✅ | |
| $g_*(T)$ interpolation | ✅ | 15-point table, reasonable |
| $x_\text{fo} \approx 23.6$ | ✅ | Standard for WIMP |
| $Y_\chi(\infty) = 3.00 \times 10^{-12}$ | ✅ | |
| $\Omega_\chi h^2 = 0.0171$ | ✅ | Low because α_s×α_p not re-tuned |

$\Omega = 0.017$ vs target 0.120: צפוי — BP1 כויל עבור V10 ($\sigma v = \pi\alpha^2/m^2$), V11 נותן 8× יותר.

---

## 3. Reference Case (μ₃ = 0)

אימות עצמאי:

$$n_\text{eq}(T=m_\phi) = \left(\frac{m_\phi^2}{2\pi}\right)^{3/2} e^{-1} = 3.40 \times 10^{-8} \text{ GeV}^3$$

$$s(m_\phi) = \frac{2\pi^2}{45} \times 10.75 \times m_\phi^3 = 6.87 \times 10^{-6} \text{ GeV}^3$$

$$Y_\phi = 4.95 \times 10^{-3} \quad \checkmark$$

$$\Omega_\phi h^2 = m_\phi Y_\phi s_0 / \rho_{c,h^2} = 1.54 \times 10^4 \quad \checkmark$$

**Overclosure factor 128,416× — correct.** This is the core of the φ problem.

---

## 4. Operator-Splitting Approach

**Valid.** χ freezes out at $T \sim 1$ GeV, cannibal acts at $T \sim 10$ MeV. Timescale separation $\sim 100×$ justifies decoupling.

Phase 2 uses "instant cannibal freeze-out" ($\Gamma$ vs $H$ threshold), analogous to Lee-Weinberg approximation. Qualitatively correct but misses $O(1)$ factors. Acceptable for a sensitivity scan.

---

## 5. ξ Evolution

A's implementation:
```python
xi *= math.exp(-dlnT / 3.0)  # d(ln ξ)/d(ln T) = -1/3
```

This gives $\xi \propto T^{-1/3}$, i.e., $T_\text{dark} \propto T_\text{vis}^{2/3}$.

**Actual physics** (Carlson, Hall, Machacek 1992): During cannibal regime, $T_\text{dark} \propto 1/\ln(a)$ (logarithmic cooling), which is much more complex than a power law.

**Assessment: Approximate.** A's parametrization captures the qualitative behavior: ξ grows during cannibal phase (3→2 heats dark sector), then freezes after cannibal decoupling. The precise values of ξ are not reliable, but the conclusion $\xi \ll 1$ is robust since even the approximate values give $\xi \in [0.07, 0.17]$.

---

## 6. ΔN_eff = 0

```python
z_bbn = M_PHI / T_dark_bbn  # ≈ 11.34/0.17 ≈ 67
if z_bbn < 3: dneff = (4/7) * xi**4 * g_phi  # never triggered
else: dneff = 0.0
```

**Correct.** φ is deeply NR at BBN ($z \approx 67 \gg 3$), so it contributes to $\Omega_m$, not to radiation ($N_\text{eff}$). The tiny kinetic energy contributes $\Delta N_\text{eff} \sim (T_d/m_\phi)^{5/2} e^{-m/T_d}$, which is negligible.

---

## 7. באג קל: אי-עקביות ב-Ω_φ for small μ₃

| כמות | Reference case | Scan (small μ₃) |
|------|---------------|-----------------|
| $\Omega_\phi h^2$ | 15,400 | 19.0 |

**הסיבה:** עבור $\mu_3$ קטן (cannibal never active), הקוד מעדכן $Y_\phi = Y_\text{eq}(T)$ בכל צעד עד $T_\text{BBN}$. אבל פיזיקלית, φ מתנתק מה-bath ב-$T \sim m_\phi$, ו-$Y_\phi$ קופא שם (ערך גבוה הרבה יותר).

**Verification:**
$$Y_\text{eq}(T_\text{BBN}) = n_\text{eq}(1 \text{ MeV})/s(1 \text{ MeV})$$
$$z = m_\phi/T_\text{BBN} = 11.34 \implies Y \approx 6 \times 10^{-6}$$
$$\Omega = m_\phi \times 6 \times 10^{-6} \times s_0/\rho_{c,h^2} \approx 19 \quad \checkmark$$

**Impact:** אפס — כל הנקודות עם $\mu_3$ קטן הן OVERCLOSED בין כה וכה.

---

## 8. Discretization Artifact

T_cann_fo = 55.79 MeV appears for 13 consecutive μ₃ values (6.87e-5 to 4.50e-4):

$$T_\text{arr} = \text{geomspace}(5m_\phi, T_\text{BBN}, 500) \implies \Delta \log T \approx 0.008$$

The cannibal "activates" at the same grid point for a range of μ₃. This is expected from the instant–freeze-out method with finite resolution. Not a bug — just a limitation.

---

## 9. Physical Conclusions

| מסקנה | סטטוס | הערה |
|-------|--------|------|
| $\mu_3 \ll m_\phi$ → overclosure | ✅ | Robust |
| $\mu_3 \gtrsim$ few × $m_\phi$ → φ depleted | ✅ qualitative | Threshold quantitatively wrong (see §1) |
| $\Omega_\chi$ independent of $\mu_3$ | ✅ | χ freezes out at $T \gg m_\phi$ |
| SIDM independent of $\mu_3$ | ✅ | σ/m depends on α_s, m_χ, m_φ only |
| One-sided bound, no fine-tuning | ✅ | Any $\mu_3$ above threshold works |
| $\Delta N_\text{eff} \approx 0$ | ✅ | φ deeply NR at BBN |

---

## 10. סיכום

| קטגוריה | ממצא | חומרה |
|---------|-------|--------|
| **נוסחת $\langle\sigma v^2\rangle$** | $m^{-9}$ במקום $m^{-11}$ + חסר $(4\pi)^{-6}$ | **CRITICAL** — rate overestimated ×500 |
| **Phase 1 Boltzmann** | נכון | — |
| **Reference overclosure** | $1.54 \times 10^4$ — מאומת | — |
| **SV0 formula** | $2\pi\alpha_s\alpha_p/m^2$ ✅ | — |
| **ξ evolution** | Approximate power-law vs exact log | MINOR — qualitative OK |
| **Small-μ₃ consistency** | Ω_φ = 19 vs 15400 in reference | MINOR — both overclosed |
| **Conclusions** | Qualitatively correct | Threshold shifts to $\mu_3/m_\phi \gtrsim 12$ |

### המלצה:
1. **תקן** את נוסחת $\langle\sigma v^2\rangle$ ל:
   $$\langle\sigma v^2\rangle = \frac{25\sqrt{5}}{512\pi(4\pi)^6}\frac{\mu_3^6}{m_\phi^{11}}$$
2. **הרץ מחדש** את הסריקה עם הנוסחה המתוקנת
3. **עדכן** את הסף $\mu_3/m_\phi$ בהתאם (צפי: ~12 במקום ~4.4)
4. המסקנה העקרונית — "מספיק $\mu_3 \gtrsim \text{few} \times m_\phi$" — **שורדת**
