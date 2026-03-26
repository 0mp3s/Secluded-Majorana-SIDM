# ולידציה: Opus A — Condition 1 (Amplitude $\chi\chi\to\phi\phi$)

**בודק:** Opus B  
**תאריך:** 23 מרץ 2026  
**קבצים נבדקים:**
- `mixed_coupling/opusA/condition1_amplitude.py`
- `mixed_coupling/opusA/condition1_output.txt`
- `mixed_coupling/opusA/thoughts_condition1.md`

---

## תוצאת ולידציה

### ✅ קוד A — נכון. התוצאה המספרית $a_0 = 2\pi\alpha_s\alpha_p/m_\chi^2$ מאומתת.

### ❌ נמצאה שגיאה בגזירה של B עצמו (פירוט בסוף המסמך)

---

## 1. בדיקת מטריצות $\gamma$ (שורות 14–48)

| פריט | סטטוס |
|------|--------|
| Dirac representation, מטריקה $(+,-,-,-)$ | ✅ |
| $\gamma^0 = \text{diag}(1,1,-1,-1)$ | ✅ |
| $\gamma^i$ — ייצוג $\begin{pmatrix}0 & \sigma^i \\ -\sigma^i & 0\end{pmatrix}$ | ✅ |
| $\gamma_5 = i\gamma^0\gamma^1\gamma^2\gamma^3$ | ✅ |
| Self-test: $\{\gamma^\mu,\gamma^\nu\} = 2\eta^{\mu\nu}$ | ✅ (4×4 = 16 cases) |
| Self-test: $\gamma_5^2 = I$ | ✅ |
| Self-test: $\{\gamma_5, \gamma^\mu\} = 0$ | ✅ (4 cases) |

## 2. Slash convention (שורה 53)

```python
def slash(p):
    return p[0]*g0 - p[1]*g1 - p[2]*g2 - p[3]*g3
```

עבור וקטור קונטרה-וריאנטי $p = (E, p_x, p_y, p_z)$:
$$\not{p} = p_\mu\gamma^\mu = p^0\gamma^0 - p^1\gamma^1 - p^2\gamma^2 - p^3\gamma^3 \quad \checkmark$$

## 3. Vertex factor (שורה 77)

```python
G = ys * I4 + 1j * yp * g5
```

$$\Gamma = y_s\mathbf{1} + iy_p\gamma_5 \quad \checkmark$$

תואם ל-Lagrangian $\mathcal{L} = -\frac{1}{2}\bar\chi(y_s + iy_p\gamma_5)\chi\phi$.

## 4. נוסחת Trace (שורות 66–96)

$$\sum_{\text{spins}}|M|^2 = \frac{T_{tt}}{D_t^2} + \frac{T_{uu}}{D_u^2} + \frac{2\text{Re}(T_{tu})}{D_t D_u}$$

$$T_{xy} = \text{Tr}\left[(\not{p}_2 - m)\,\Gamma\,({\not{q}_x + m})\,\Gamma\,(\not{p}_1 + m)\,\Gamma\,({\not{q}_y + m})\,\Gamma\right]$$

**אימות:** לצימוד $\Gamma = y_s + iy_p\gamma_5$, מתקיים:
- $\Gamma^\dagger = y_s - iy_p\gamma_5$ (כי $\gamma_5^\dagger = \gamma_5$)
- $\tilde\Gamma \equiv \gamma^0\Gamma^\dagger\gamma^0 = y_s + iy_p\gamma_5 = \Gamma$ (כי $\gamma^0\gamma_5\gamma^0 = -\gamma_5$)
- לכן $({\Gamma N \Gamma})^\sim = \Gamma N \Gamma$ — הנוסחה נכונה ✅

**סימן Bose:** $M = M_t + M_u$ (פלוס) — נכון לבוזונים זהים, כי $k_1 \leftrightarrow k_2$ ✅

## 5. Phase space (שורות 99–127)

```python
return kf * integral / (512.0 * np.pi * E**3)
```

**גזירה:**

$$\sigma v_{\rm Möl} = \frac{1}{64\pi^2 s}\,\frac{k_f}{k_i}\,\frac{1}{S!}\,\int d\Omega\,\frac{1}{4}\sum|M|^2 $$

עם $s = 4E^2$, $v_{\rm Möl} = 2k_i/E$, $S! = 2$ (זהים), ו-$\int d\Omega = 2\pi\int_{-1}^{1}d\cos\theta$:

$$\sigma v = \frac{k_f}{512\pi E^3}\int_{-1}^{1}d\cos\theta\;\sum|M|^2 \quad \checkmark$$

## 6. תוצאות מספריות

### Part 1: CP selection rule — ✅

| $(y_s, y_p)$ | $a_0$ | צפוי | סטטוס |
|---|---|---|---|
| $(1, 0)$ | $5.66\times10^{-11}$ | $\approx 0$ | ✅ |
| $(0, 1)$ | $1.72\times10^{-12}$ | $\approx 0$ | ✅ |
| $(1, 1)$ | $9.295\times10^{-5}$ | $\neq 0$ | ✅ |

שבירת CP עובדת: pure scalar/pseudoscalar → $a_0 = 0$, mixed → $a_0 \neq 0$.

### Part 2: Scaling test $a_0 \propto y_s^2 y_p^2$ — ✅

9/9 נקודות עם ratio = 1.00000 ± 10⁻⁵. דגימות מספיקות ($y = 0.5, 1, 2, 3$). Tests both $y_s \neq y_p$ and $y_s = y_p$.

### Part 3: Analytical comparison — ✅

$c_{\rm LO} = 1/(8\pi m_\chi^2) = 9.2948 \times 10^{-5}$ GeV⁻² — התאמה ב-$5\times 10^{-7}$.

**⚠️ שגיאה קלה ב-`c_exact`:** A כותב:
```python
c_exact = beta_f * m_chi / (2.0 * np.pi * (2*m_chi**2 - m_phi**2)**2)
```
הנוסחה הנכונה היא $c = \beta_f m_\chi^2 / (2\pi D^2)$ — חסר גורם $m_\chi$ (צריך $m_\chi^2$ במונה, לא $m_\chi$).
זה מסביר את ה-deviation של 19.7× שנמצא. **אבל `c_LO` נכון**, וזו הנוסחה שנבדקת.

### Part 4: Physical benchmark — ✅

$\alpha_{\rm total} = 1.048 \times 10^{-3}$, $\alpha_s = \alpha_p = 5.24 \times 10^{-4}$:
$$a_0 = 4.03 \times 10^{-9}\;\text{GeV}^{-2} = 1.55 \times \text{canonical}$$

מאשר ש-mixed Majorana עם symmetric split נותן s-wave בסדר גודל הנכון.

### Part 5: Velocity dependence — ✅

$\sigma v / a_0 = 0.9999$ at $v = 0.001$, falling to $0.935$ at $v = 0.3$. Confirms s-wave dominance at $v \ll 1$.

### Part 6: Direct threshold — ✅

$a_0(v=0) = 4.030153 \times 10^{-9}$, $a_0(\text{fit}) = 4.030154 \times 10^{-9}$. Agreement: $<10^{-4}\%$.

## 7. מסמך thoughts_condition1.md

- ✅ הכרת הסימן Bose (M_t + M_u) — ומציין שפלוס מביא ל-a₀(1,0) = 0 כצפוי
- ✅ פיזיקת CP מוסברת נכון
- ✅ הנוסחה $a_0 = y_s^2y_p^2/(8\pi m_\chi^2)$ מנוסחת נכון
- ✅ השקילות ל-$2\pi\alpha_s\alpha_p/m_\chi^2$ — נכון: $y^2 = 4\pi\alpha$ אז $16\pi^2\alpha_s\alpha_p/(8\pi m^2) = 2\pi\alpha_s\alpha_p/m^2$

---

## 8. ❌ שגיאה שנמצאה בגזירה של B

בעת הוולידציה של A, גיליתי שגיאת אריתמטיקה **בגזירה שלי עצמי** (`opusB/amplitude_derivation.md`).

### מיקום: סעיף 4, משוואת הקופסה

**B כתב:**
$$\mathcal{M}_0 = \frac{\color{red}{-4i}m_\chi^2\,y_s y_p}{m_\phi^2 - 2m_\chi^2}\,\eta^\dagger\xi$$

**הנכון:**
$$\mathcal{M}_0 = \frac{\color{green}{-8i}m_\chi^2\,y_s y_p}{m_\phi^2 - 2m_\chi^2}\,\eta^\dagger\xi$$

### מקור השגיאה

מהגזירה:
$$\mathcal{M}_0 = \frac{2m_\chi}{D_0}\left[2iy_sy_p \cdot (-2m_\chi\,\eta^\dagger\xi)\right] = \frac{2m_\chi \cdot 2i \cdot (-2m_\chi)}{D_0}\,y_sy_p\,\eta^\dagger\xi = \frac{-8im_\chi^2}{D_0}\,y_sy_p\,\eta^\dagger\xi$$

$(2)(2i)(-2) = -8i$, אבל B כתב $-4i$. **חסר גורם 2.**

### השפעה במורד חישוב

| נוסחה | B (שגוי) | נכון |
|---|---|---|
| $\|\mathcal{M}_0\|^2$ | $16m^4 y_s^2 y_p^2 \|\eta^\dagger\xi\|^2$ | $64m^4 y_s^2 y_p^2\|\eta^\dagger\xi\|^2$ |
| $\sum\|\mathcal{M}_0\|^2$ | $32 m^4/D_0^2$ | $128 m^4/D_0^2$ |
| $\overline{\|\mathcal{M}_0\|^2}$ | $8 m^4/D_0^2$ | $32 m^4/D_0^2$ |
| $(\sigma v)_0$ (LO) | $y_s^2y_p^2/(32\pi m^2)$ | $y_s^2y_p^2/(8\pi m^2)$ |
| $a_0$ | $\pi\alpha_s\alpha_p/(2m^2)$ | $2\pi\alpha_s\alpha_p/m^2$ |

**B's formula is 4× too small.**

### למה ה-cross-check בסעיף 5 של B לא תפס?

B כתב:
$$\overline{|\mathcal{M}_0|^2} = -\frac{m_\chi^4}{D_0^2}\text{Tr}[P_-\Gamma^2 P_+\Gamma^2]$$

**הנכון:**
$$\overline{|\mathcal{M}_0|^2} = -\frac{{\color{green}4}m_\chi^4}{D_0^2}\text{Tr}[P_-\Gamma^2 P_+\Gamma^2]$$

הגורם $4m^4$ מגיע מ: $\overline{|M|^2} = (4m^2/D_0^2)(1/4) \text{Tr}[(-2mP_-)\Gamma^2(2mP_+)\Gamma^2]$  
$= (m^2/D_0^2)(-4m^2)\text{Tr}[P_-\Gamma^2P_+\Gamma^2] = -4m^4/D_0^2 \cdot \text{Tr}[...]$

B כתב $-m^4$ במקום $-4m^4$ — **אותה שגיאת גורם-4**, מה שגרם ל-cross-check להיראות עקבי.

### אימות: A מתאים, B לא

A's numerical result: $a_0(1,1) = 9.2948 \times 10^{-5}$ GeV⁻²

| Formula | Value | Match? |
|---|---|---|
| A: $y_s^2y_p^2/(8\pi m^2)$ | $9.2948 \times 10^{-5}$ | ✅ exact |
| B: $y_s^2y_p^2/(32\pi m^2)$ | $2.324 \times 10^{-5}$ | ❌ 4× too small |

---

## 9. סיכום ולידציה

| פריט | A | B |
|------|---|---|
| שיטה | מספרית (4×4 γ-matrices) | אנליטית (spinor algebra) |
| Self-tests | 3 Clifford + 9 scaling | 1 trace cross-check |
| $a_0$ formula | $2\pi\alpha_s\alpha_p/m^2$ ✅ | $\pi\alpha_s\alpha_p/(2m^2)$ ❌ |
| Match numerics | ✅ (5×10⁻⁷) | ❌ (4× off) |
| CP selection rule | ✅ verified | ✅ qualitative argument correct |
| Bug found | c_exact typo (m→m²), cosmetic | Arithmetic error ×4 in $\mathcal{M}_0$ |

### הנוסחה הנכונה (A's):

$$\boxed{(\sigma v)_0 = \frac{y_s^2\,y_p^2}{8\pi\,m_\chi^2} = \frac{2\pi\,\alpha_s\,\alpha_p}{m_\chi^2}}$$

---

*Opus B, ולידציה של Opus A, 23 מרץ 2026*
