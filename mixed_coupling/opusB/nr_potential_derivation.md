# NR Potential for Mixed Yukawa Coupling: $\chi\chi$ Elastic Scattering

**תנאי 4 מתוך 4 — אישור שה-VPM solver תקף ללא שינוי**

---

## 1. שאלה

המודל V11 כולל צימוד מעורב:

$$\mathcal{L}_{\rm int} = -\frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\,\phi$$

ה-VPM solver מחשב פיזור $\chi\chi \to \chi\chi$ בפוטנציאל Yukawa:

$$V_{\rm VPM}(r) = -\alpha\frac{e^{-m_\phi r}}{r}$$

**שאלה:** האם $y_p$ משנה את הפוטנציאל? אם כן — ב-כמה?

---

## 2. Born Amplitude for $\chi(p_1)\chi(p_2) \to \chi(p_3)\chi(p_4)$

### Diagram: t-channel $\phi$ exchange

$$i\mathcal{M}_t = \frac{i}{t - m_\phi^2}\left[\bar{u}(p_3)\Gamma u(p_1)\right]\left[\bar{u}(p_4)\Gamma u(p_2)\right]$$

with $\Gamma = y_s + iy_p\gamma_5$ and $t = (p_1 - p_3)^2 = -|\vec{q}|^2$ in the NR limit.

### NR reduction of vertex $\bar{u}(p')\Gamma u(p)$

In the NR limit ($|\vec{p}|/m_\chi \ll 1$), using Dirac representation:

$$u(p,s) \approx \sqrt{2m_\chi}\begin{pmatrix}\xi_s \\ \frac{\vec{\sigma}\cdot\vec{p}}{2m_\chi}\xi_s\end{pmatrix}$$

**Scalar part:**

$$\bar{u}'\,y_s\,u = y_s\,u'^\dagger\gamma^0 u \approx 2m_\chi\,y_s\left(\xi'^\dagger\xi + O(p^2/m_\chi^2)\right)$$

Leading term: $\sim 2m_\chi y_s$ — **momentum-independent** → Yukawa potential.

**Pseudoscalar part:**

$$\bar{u}'\,iy_p\gamma_5\,u = iy_p\,u'^\dagger\gamma^0\gamma_5\,u$$

$$\gamma^0\gamma_5 = \begin{pmatrix}0 & I \\ -I & 0\end{pmatrix}$$

$$\bar{u}'\gamma_5 u \approx 2m_\chi\left(\xi'^\dagger\frac{\vec{\sigma}\cdot\vec{p}}{2m_\chi}\xi - \frac{(\vec{\sigma}\cdot\vec{p}')^\dagger}{2m_\chi}\xi'^\dagger\cdot\xi\right)$$

$$= \xi'^\dagger\vec{\sigma}\cdot(\vec{p} - \vec{p}')\xi = -\xi'^\dagger(\vec{\sigma}\cdot\vec{q})\xi$$

Leading term: $\sim y_p\,|\vec{q}|/m_\chi \sim y_p\,m_\phi/m_\chi$ — **suppressed by $m_\phi/m_\chi$**.

---

## 3. Potential in Coordinate Space

### Scalar exchange (leading)

$$V_s(\vec{q}) = -\frac{y_s^2}{|\vec{q}|^2 + m_\phi^2} \xrightarrow{\text{FT}} V_s(r) = -\frac{\alpha_s}{r}e^{-m_\phi r}$$

**This is the standard Yukawa potential used by VPM.**

### Pseudoscalar exchange

$$V_p(\vec{q}) \propto \frac{y_p^2(\vec{\sigma}_1\cdot\vec{q})(\vec{\sigma}_2\cdot\vec{q})}{|\vec{q}|^2 + m_\phi^2}$$

Fourier transforming (see e.g. Donoghue, Golowich & Holstein, or any QFT textbook):

$$V_p(r) = -\frac{\alpha_p}{12m_\chi^2}\left[4\pi\delta^3(\vec{r})(\vec{S}_1\cdot\vec{S}_2) + m_\phi^2\frac{e^{-m_\phi r}}{r}\left((\vec{S}_1\cdot\vec{S}_2) + S_{12}\left(1 + \frac{3}{m_\phi r} + \frac{3}{(m_\phi r)^2}\right)\right)\right]$$

where $S_{12} = 3(\vec{S}_1\cdot\hat{r})(\vec{S}_2\cdot\hat{r}) - \vec{S}_1\cdot\vec{S}_2$ is the tensor operator.

### Scalar–pseudoscalar cross-term

$$V_{sp}(\vec{q}) \propto \frac{y_s y_p\,(\vec{\sigma}\cdot\vec{q})}{|\vec{q}|^2 + m_\phi^2}$$

This is **spin-orbit** type: $V_{sp}(r) \propto \frac{\alpha_s^{1/2}\alpha_p^{1/2}}{m_\chi}\frac{1}{r}\frac{d}{dr}\frac{e^{-m_\phi r}}{r}\,\vec{L}\cdot\vec{S}$

Also suppressed by $1/m_\chi$.

---

## 4. Power Counting

| Contribution | Characteristic scale | Ratio to $V_s$ |
|---|---|---|
| $V_s$ (scalar Yukawa) | $\alpha_s/r$ | **1** |
| $V_p$ (pseudoscalar) | $\frac{\alpha_p\,m_\phi^2}{m_\chi^2\,r}$ | $(m_\phi/m_\chi)^2 \sim 3 \times 10^{-7}$ |
| $V_{sp}$ (cross-term) | $\frac{\sqrt{\alpha_s\alpha_p}\,m_\phi}{m_\chi\,r}$ | $m_\phi/m_\chi \sim 5 \times 10^{-4}$ |

For BP1 ($m_\chi = 13.9$ GeV, $m_\phi = 7.3$ MeV):
- $m_\phi/m_\chi = 5.3 \times 10^{-4}$
- $(m_\phi/m_\chi)^2 = 2.8 \times 10^{-7}$

**Even the largest correction (spin-orbit) is $5 \times 10^{-4}$ relative to the leading Yukawa.**

### After spin-averaging

The tensor term $S_{12}$ averages to zero over spin orientations. The $\vec{S}_1\cdot\vec{S}_2$ term gives $+1/4$ (triplet) or $-3/4$ (singlet), but weighted by $(m_\phi/m_\chi)^2$. The spin-orbit $\vec{L}\cdot\vec{S}$ also averages to zero for s-wave.

**Effective NR potential for VPM (spin-averaged, s-wave dominated):**

$$V_{\rm eff}(r) = -\frac{\alpha_s}{r}e^{-m_\phi r}\left[1 + O\left(\frac{m_\phi^2}{m_\chi^2}\right)\right]$$

---

## 5. Conclusion

> **הפוטנציאל ה-NR נשלט לחלוטין ע"י הצימוד הסקלרי $\alpha_s$.**
> תיקוני $y_p$ מדוכאים ב-$(m_\phi/m_\chi)^2 \sim 10^{-7}$ (לאחר ממוצע ספין).
> ה-VPM solver משמש as-is עם $\alpha_{\rm eff} = \alpha_s$.

---

## 5. Numerical Verification

Script: `nr_potential_check.py`

### Part A: Sensitivity — VPM responds to $\alpha$

Doubling $\alpha$ (extreme upper bound — pretending all of $y_p$ adds directly):

| BP | $m_\chi$ (GeV) | $m_\phi$ (MeV) | $\alpha$ | $v$ (km/s) | $\sigma_T/m(\alpha)$ | $\sigma_T/m(2\alpha)$ | $\Delta$ |
|----|------|------|------|------|------|------|------|
| BP1 | 13.9 | 7.3 | $6.4\times10^{-4}$ | 30 | 1.001 | 2.393 | 139% |
| BP2 | 18.4 | 8.2 | $8.2\times10^{-4}$ | 30 | 1.061 | 2.246 | 112% |
| BP3 | 24.4 | 9.9 | $1.2\times10^{-3}$ | 30 | 0.940 | 1.768 | 88% |

VPM is strongly sensitive to $\alpha$ — as expected.

### Part B: True $y_p$ correction — $\alpha \to \alpha(1 + \varepsilon)$, $\varepsilon = (m_\phi/m_\chi)^2$

| BP | $\varepsilon$ | $v$ (km/s) | $\delta(\sigma_T/m)/(\sigma_T/m)$ | Pass? |
|----|------|------|------|------|
| BP1 | $2.8\times10^{-7}$ | 10 | $4.4\times10^{-5}\%$ | ✅ |
| BP1 | | 30 | $3.9\times10^{-5}\%$ | ✅ |
| BP1 | | 100 | $3.9\times10^{-5}\%$ | ✅ |
| BP1 | | 1000 | $5.4\times10^{-5}\%$ | ✅ |
| BP2 | $2.0\times10^{-7}$ | 10 | $2.7\times10^{-5}\%$ | ✅ |
| BP2 | | 30 | $2.4\times10^{-5}\%$ | ✅ |
| BP2 | | 100 | $2.5\times10^{-5}\%$ | ✅ |
| BP2 | | 1000 | $3.8\times10^{-5}\%$ | ✅ |
| BP3 | $1.6\times10^{-7}$ | 10 | $1.9\times10^{-5}\%$ | ✅ |
| BP3 | | 30 | $1.7\times10^{-5}\%$ | ✅ |
| BP3 | | 100 | $1.8\times10^{-5}\%$ | ✅ |
| BP3 | | 1000 | $3.0\times10^{-5}\%$ | ✅ |

**Maximum $\delta\sigma_T/\sigma_T = 5.4 \times 10^{-5}\%$** — six orders of magnitude below 1%.

### Verdict

$$\boxed{\text{CONDITION 4 PASSED} — \delta\sigma_T/\sigma_T < 10^{-6}}$$

---

## 6. השלכה על מרחב הפרמטרים

ב-V10, הפרמטר $\alpha$ שימש **גם** לאניהילציה **וגם** לפיזור. ב-V11:

- **אניהילציה:** $(\sigma v)_0 \propto \alpha_s \alpha_p$ — תלויה בשני הצימודים
- **פיזור:** $\sigma_T \propto \alpha_s$ בלבד (+ תיקונים זניחים מ-$\alpha_p$)

המשמעות: $\alpha_p$ הוא פרמטר "חופשי" שמשפיע רק על relic density, לא על SIDM. זה מרחיב את ה-island מקו אחד $\alpha = \alpha_{\rm relic}(m_\chi)$ ל**רצועה** בתוך ההיפרבולה $\alpha_s\alpha_p = \text{const}$.

---

*B, 23 מרץ 2026*
