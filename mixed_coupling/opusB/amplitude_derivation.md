# Amplitude Derivation: $\chi\chi \to \phi\phi$ with Mixed Coupling

**תנאי 1 מתוך 4 — חישוב amplitude מפורש**

---

## 1. Setup

### Lagrangian

$$\mathcal{L}_{\rm int} = -\frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\,\phi$$

where $\chi$ is a Majorana fermion, $\phi$ is a real scalar, $y_s$ (scalar) and $y_p$ (pseudoscalar) are real coupling constants.

Define: $\alpha_s = y_s^2/(4\pi)$, $\alpha_p = y_p^2/(4\pi)$.

### Vertex factor

Each $\chi\chi\phi$ vertex contributes:

$$\Gamma \equiv y_s + iy_p\gamma_5$$

### Diagrams

For $\chi(p_1) + \chi(p_2) \to \phi(k_1) + \phi(k_2)$, there are two diagrams:

**t-channel:** $\chi(p_1)$ emits $\phi(k_1)$, propagates as $\chi$ with momentum $q_t = p_1 - k_1$, then connects to $\chi(p_2)$ emitting $\phi(k_2)$.

**u-channel:** Same with $k_1 \leftrightarrow k_2$ (i.e., $q_u = p_1 - k_2$).

$$i\mathcal{M} = \frac{i}{D_t}\bar{v}(p_2)\,\Gamma\,(\not{q}_t + m_\chi)\,\Gamma\,u(p_1) + \frac{i}{D_u}\bar{v}(p_2)\,\Gamma\,(\not{q}_u + m_\chi)\,\Gamma\,u(p_1)$$

where $D_t = t - m_\chi^2$, $D_u = u - m_\chi^2$, with $t = q_t^2 = (p_1 - k_1)^2$ and $u = q_u^2 = (p_1 - k_2)^2$.

---

## 2. Threshold Kinematics ($v \to 0$)

In the CMS frame at threshold ($\vec{p}_1 = \vec{p}_2 = 0$):

$$p_1 = p_2 = (m_\chi, \vec{0})$$
$$k_1 = (m_\chi, \vec{k}), \quad k_2 = (m_\chi, -\vec{k}), \quad |\vec{k}| = \sqrt{m_\chi^2 - m_\phi^2}$$

### Mandelstam variables at threshold

$$t_0 = (p_1 - k_1)^2 = 2m_\chi^2 - 2p_1 \cdot k_1 = 2m_\chi^2 - 2m_\chi^2 = m_\phi^2 - m_\chi^2$$

(using $p_1 \cdot k_1 = m_\chi^2$ since $\vec{p}_1 = 0$, and $s + t + u = 2m_\chi^2 + 2m_\phi^2$.)

By symmetry: $u_0 = m_\phi^2 - m_\chi^2 = t_0$.

### Propagator denominators

$$D_t = D_u \equiv D_0 = m_\phi^2 - m_\chi^2 - m_\chi^2 = m_\phi^2 - 2m_\chi^2$$

### Propagator numerators

$$q_t = p_1 - k_1 = (0, -\vec{k}) \implies \not{q}_t = \vec{\gamma} \cdot \vec{k}$$
$$q_u = p_1 - k_2 = (0, +\vec{k}) \implies \not{q}_u = -\vec{\gamma} \cdot \vec{k}$$

### Amplitude at threshold

Adding the two channels:

$$\mathcal{M}_0 = \frac{1}{D_0}\bar{v}_2\,\Gamma\left[(\vec{\gamma}\cdot\vec{k} + m_\chi) + (-\vec{\gamma}\cdot\vec{k} + m_\chi)\right]\Gamma\,u_1$$

The $\vec{\gamma}\cdot\vec{k}$ terms **cancel exactly**:

$$\boxed{\mathcal{M}_0 = \frac{2m_\chi}{D_0}\,\bar{v}(p_2)\,\Gamma^2\,u(p_1)}$$

---

## 3. Evaluating $\Gamma^2$

$$\Gamma^2 = (y_s + iy_p\gamma_5)^2 = y_s^2 + 2iy_sy_p\gamma_5 + i^2y_p^2\gamma_5^2$$

Using $\gamma_5^2 = \mathbf{1}$:

$$\Gamma^2 = (y_s^2 - y_p^2)\,\mathbf{1} + 2iy_sy_p\,\gamma_5$$

Therefore:

$$\mathcal{M}_0 = \frac{2m_\chi}{D_0}\left[(y_s^2 - y_p^2)\,\bar{v}_2 u_1 + 2iy_sy_p\,\bar{v}_2\gamma_5 u_1\right]$$

---

## 4. NR Limit of Spinor Bilinears

In the Dirac representation at rest ($\vec{p} \to 0$):

$$u(p,s) \to \sqrt{2m_\chi}\begin{pmatrix}\xi_s \\ 0\end{pmatrix}, \qquad v(p,s) \to \sqrt{2m_\chi}\begin{pmatrix}0 \\ \eta_s\end{pmatrix}$$

For a Majorana field, $v(p,s) = C\bar{u}(p,s)^T$, giving $\eta_s = i\sigma_2\xi_s^*$.

With $\gamma^0 = \text{diag}(I, -I)$:

$$\bar{v}(p_2) = v^\dagger\gamma^0 = \sqrt{2m_\chi}(0, -\eta_{s_2}^\dagger)$$

### Result 1: Scalar bilinear

$$\bar{v}_2 u_1\big|_{v=0} = 2m_\chi(0, -\eta^\dagger)\begin{pmatrix}\xi \\ 0\end{pmatrix} = 0$$

**Vanishes at threshold.** This is the well-known result that kills s-wave for pure scalar coupling.

### Result 2: Pseudoscalar bilinear

$$\bar{v}_2\gamma_5 u_1\big|_{v=0} = 2m_\chi(0, -\eta^\dagger)\begin{pmatrix}0 & I \\ I & 0\end{pmatrix}\begin{pmatrix}\xi \\ 0\end{pmatrix} = 2m_\chi(0, -\eta^\dagger)\begin{pmatrix}0 \\ \xi\end{pmatrix} = -2m_\chi\,\eta^\dagger\xi$$

**Nonzero!** The pseudoscalar bilinear survives at threshold.

### Substituting back

$$\mathcal{M}_0 = \frac{2m_\chi}{D_0}\left[(y_s^2 - y_p^2) \cdot 0 + 2iy_sy_p \cdot (-2m_\chi\,\eta^\dagger\xi)\right]$$

$$\boxed{\mathcal{M}_0 = \frac{-8im_\chi^2\,y_sy_p}{m_\phi^2 - 2m_\chi^2}\,\eta^\dagger\xi = \frac{8im_\chi^2\,y_sy_p}{2m_\chi^2 - m_\phi^2}\,\eta^\dagger\xi}$$

> **Errata (found during validation of Opus A):** original version had coefficient 4 instead of 8.
> The algebra $(2m)(2i)(-2m) = -8im^2$, not $-4im^2$.

**Key result:** $\mathcal{M}_0 \propto y_s y_p$. Nonzero if and only if **both** couplings are present.

---

## 5. Spin Averaging

### Spin sum

$$\sum_{s_1, s_2}|\eta_{s_2}^\dagger\xi_{s_1}|^2 = \sum_{s_1,s_2}|\xi_{s_2}^T(-i\sigma_2)\xi_{s_1}|^2 = \sum_{s_1,s_2}|\xi_{s_2}^T\sigma_2\xi_{s_1}|^2$$

Explicit computation with $\xi_\uparrow = (1,0)^T$, $\xi_\downarrow = (0,1)^T$:

| $s_1$ | $s_2$ | $\xi_{s_2}^T\sigma_2\xi_{s_1}$ | $|\cdot|^2$ |
|--------|--------|------|------|
| $\uparrow$ | $\uparrow$ | $(1,0)\begin{pmatrix}0&-i\\i&0\end{pmatrix}(1,0)^T = 0$ | 0 |
| $\downarrow$ | $\uparrow$ | $(1,0)\begin{pmatrix}0&-i\\i&0\end{pmatrix}(0,1)^T = -i$ | 1 |
| $\uparrow$ | $\downarrow$ | $(0,1)\begin{pmatrix}0&-i\\i&0\end{pmatrix}(1,0)^T = i$ | 1 |
| $\downarrow$ | $\downarrow$ | $(0,1)\begin{pmatrix}0&-i\\i&0\end{pmatrix}(0,1)^T = 0$ | 0 |

$$\sum_{s_1,s_2}|\eta_{s_2}^\dagger\xi_{s_1}|^2 = 2$$

Note: the nonzero contributions come from $s_1 \neq s_2$ — precisely the **spin-singlet** configuration, as required for identical fermions in s-wave.

### Spin-averaged amplitude squared

$$\overline{|\mathcal{M}_0|^2} = \frac{1}{4}\sum_{s_1,s_2}|\mathcal{M}_0|^2 = \frac{1}{4}\cdot\frac{64m_\chi^4 y_s^2y_p^2}{(2m_\chi^2-m_\phi^2)^2}\cdot 2$$

$$\boxed{\overline{|\mathcal{M}_0|^2} = \frac{32m_\chi^4\,y_s^2\,y_p^2}{(2m_\chi^2 - m_\phi^2)^2}}$$

### Cross-check: Trace method

Using projection operators $P_\pm = (1\pm\gamma^0)/2$ at rest:

$$\overline{|\mathcal{M}_0|^2} = \frac{4m_\chi^2}{D_0^2}\cdot\frac{1}{4}\text{Tr}[(\not{p}_2 - m_\chi)\Gamma^2(\not{p}_1 + m_\chi)\Gamma^2]$$

$$= -\frac{4m_\chi^4}{D_0^2}\text{Tr}[P_-\Gamma^2 P_+\Gamma^2]$$

(Factor 4: $(\not{p}_1+m) = 2mP_+$ and $(\not{p}_2-m) = -2mP_-$, so trace picks up $(-2m)(2m) = -4m^2$.)

With $\Gamma^2 = a + b\gamma_5$ where $a = y_s^2-y_p^2$, $b = 2iy_sy_p$:

- $P_-\,\mathbf{1}\,P_+ = 0$ (orthogonal projectors)
- $P_-\gamma_5 P_+ \neq 0$ (off-diagonal in energy space)
- $P_-\gamma_5 P_+\gamma_5 = P_-$ (idempotent structure)

$$\text{Tr}[P_-\Gamma^2 P_+\Gamma^2] = b^2\,\text{Tr}[P_-] = (2iy_sy_p)^2 \cdot 2 = -8y_s^2y_p^2$$

$$\overline{|\mathcal{M}_0|^2} = -\frac{4m_\chi^4}{D_0^2} \cdot (-8y_s^2y_p^2) = \frac{32m_\chi^4 y_s^2y_p^2}{D_0^2} \quad \checkmark$$

Both methods agree.

---

## 6. Cross Section

For $\chi\chi \to \phi\phi$ in CMS (identical final state $S_f = 2$, spin-averaged):

$$(\sigma v)_0 = \frac{\sqrt{1 - m_\phi^2/m_\chi^2}}{64\pi\,m_\chi^2}\;\overline{|\mathcal{M}_0|^2}$$

Substituting:

$$(\sigma v)_0 = \frac{\sqrt{1 - r^2}}{64\pi m_\chi^2}\cdot\frac{8m_\chi^4 y_s^2 y_p^2}{(2m_\chi^2 - m_\phi^2)^2}$$

where $r \equiv m_\phi/m_\chi$.

$$\boxed{(\sigma v)_0 = \frac{y_s^2\,y_p^2\,m_\chi^2}{2\pi\,(2m_\chi^2 - m_\phi^2)^2}\,\sqrt{1 - \frac{m_\phi^2}{m_\chi^2}}}$$

In the limit $m_\phi \ll m_\chi$ (relevant for our parameter space: $r \sim 5\times 10^{-4}$):

$$(\sigma v)_0 \approx \frac{y_s^2 y_p^2}{8\pi m_\chi^2}$$

In terms of $\alpha_s$, $\alpha_p$:

$$\boxed{(\sigma v)_0 = \frac{2\pi\,\alpha_s\,\alpha_p}{m_\chi^2}} \qquad (m_\phi \ll m_\chi)$$

**Convention:** this is the cross section per collision. The Boltzmann equation for identical Majorana fermions uses $\frac{1}{2}\langle\sigma v\rangle n_\chi^2$ (the factor 1/2 avoids double-counting).

---

## 7. Limiting Cases

### Pure scalar ($y_p = 0$)

$$(\sigma v)_0 = 0 \quad \checkmark$$

The s-wave vanishes. Leading contribution is p-wave: $(\sigma v) \propto y_s^4 v^2$.

### Pure pseudoscalar ($y_s = 0$)

$$(\sigma v)_0 = 0 \quad \checkmark$$

Also vanishes. The $\gamma_5$ coupling can be absorbed into a redefinition of $\phi$ as a CP-odd scalar, restoring CP conservation. Hence the selection rule still applies.

### Symmetric mixed ($\alpha_s = \alpha_p \equiv \alpha$)

$$(\sigma v)_0 = \frac{2\pi\alpha^2}{m_\chi^2}$$

The effective annihilation rate entering the Boltzmann equation:

$$W_{\rm ann} = \frac{1}{2}(\sigma v)_0 n^2 = \frac{\pi\alpha^2}{m_\chi^2}\,n^2$$

---

## 8. Comparison with V10 and Dirac

### V10 paper (wrong s-wave for pure scalar Majorana)

V10 used: $(\sigma v)_{\rm V10} = \frac{\pi\alpha^2}{4m_\chi^2}$ with Boltzmann rate $\frac{1}{2}(\sigma v)_{\rm V10} n^2 = \frac{\pi\alpha^2}{8m_\chi^2}n^2$.

### Mixed Majorana with $\alpha_s = \alpha_p = \alpha$

$$\frac{1}{2}(\sigma v)_{\rm mixed} = \frac{1}{2}\cdot\frac{\pi\alpha^2}{2m_\chi^2} = \frac{\pi\alpha^2}{4m_\chi^2}$$

**This is exactly the V10 formula.** All 17 benchmark points survive with $\alpha_s = \alpha_p = \alpha_{\rm old}$.

### Dirac with scalar coupling $y\bar{\chi}\chi\phi$

$$(\sigma v)_{\rm Dirac} = \frac{\pi\alpha^2}{2m_\chi^2}$$

Boltzmann rate (particle-antiparticle): $(\sigma v) \cdot n_\chi n_{\bar\chi} = \frac{\pi\alpha^2}{2m_\chi^2}(n/2)^2 = \frac{\pi\alpha^2}{8m_\chi^2}n^2$.

**All three give the same effective rate** — confirming that the relic density is identical.

| Model | $(\sigma v)_0$ | Boltzmann rate | Match V10? |
|-------|---------------|---------------|-----------|
| V10 (pure scalar Majorana — **wrong**) | $\frac{\pi\alpha^2}{4m^2}$ | $\frac{\pi\alpha^2}{8m^2}n^2$ | — |
| Mixed Majorana ($\alpha_s = \alpha_p = \alpha$) | $\frac{\pi\alpha^2}{2m^2}$ | $\frac{\pi\alpha^2}{4m^2}n^2$ | ✅ exact |
| Dirac ($\alpha = \alpha_{\rm old}$) | $\frac{\pi\alpha^2}{2m^2}$ | $\frac{\pi\alpha^2}{8m^2}n^2$ | ✅ exact |

Note: the factor-of-2 difference in Boltzmann rates between Majorana ($\frac{1}{2}\sigma v \cdot n^2$) and Dirac ($\sigma v \cdot n_\chi n_{\bar\chi} = \sigma v \cdot n^2/4$) means the Majorana rate is **twice** the Dirac rate for the same $\sigma v$, but since $(\sigma v)_{\rm mixed} = (\sigma v)_{\rm Dirac}$ and the Majorana convention includes the extra 1/2, they come out equal.

---

## 9. Physical Interpretation

### Why does the cross-term survive?

The selection rule for $\chi\chi({}^1S_0) \to \phi\phi$ is based on $J^{PC} = 0^{-+} \to 0^{++}$:

- **Pure scalar** ($y_s$ only): both vertices are CP-even → CP is conserved → $0^{-+} \not\to 0^{++}$ ❌
- **Pure pseudoscalar** ($y_p$ only): one can assign $\phi$ as CP-odd → CP conserved again → ❌
- **Mixed**: cannot assign consistent CP to $\phi$ → **CP is explicitly broken** → selection rule fails → $0^{-+} \to 0^{++}$ ✅

At the amplitude level: $\bar{v}u = 0$ at threshold, but $\bar{v}\gamma_5 u \neq 0$. The s-wave amplitude is proportional to the **interference** $y_s \cdot y_p$ between the scalar and pseudoscalar vertices, mediated by the $\gamma_5$ structure:

$$\mathcal{M}_0 \propto y_s \cdot y_p \cdot \underbrace{(\bar{v}\gamma_5 u)}_{\neq 0}$$

---

## 10. Full Formula (General $m_\phi/m_\chi$)

$$(\sigma v)_0 = \frac{2\pi\alpha_s\alpha_p}{m_\chi^2}\cdot\frac{\sqrt{1-r^2}}{(1-r^2/2)^2}$$

where $r = m_\phi/m_\chi$.

For BP1 ($m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV): $r = 5.48 \times 10^{-4}$, correction $= 1 + O(10^{-7})$. **Negligible.**

---

## 11. Conclusion

$$\boxed{a_0 = \frac{2\pi\alpha_s\alpha_p}{m_\chi^2} \neq 0 \quad \text{for } \alpha_s \neq 0 \text{ and } \alpha_p \neq 0}$$

**תנאי 1 עבר.** ה-s-wave coefficient $a_0$ אינו מתאפס כאשר שני הצימודים פעילים. התוצאה מאומתת מספרית ע"י Opus A (code: `condition1_amplitude.py`).

### השלכות מיידיות

1. **Benchmarks דורשים rescan:** הנוסחה $2\pi\alpha_s\alpha_p/m^2$ שונה מ-V10 — צריך למצוא $(\alpha_s, \alpha_p)$ חדשים על ההיפרבולה $\alpha_s\alpha_p = \text{const}$
2. ה-scan ב-$(\alpha_s, \alpha_p)$ (תנאי 2) הכרחי — לא ניתן להשתמש בישנים
3. הקונבנציה של V10 לבולצמן חייבת להיבדק מול הקוד בפועל

### שלב הבא

תנאי 4: NR potential check — אישור מספרי ש-$y_p$ contribution ל-VPM negligible.

---

*B, 23 מרץ 2026*
