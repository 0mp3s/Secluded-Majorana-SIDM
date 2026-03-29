# Dark Electromagnetic Duality: Dark Energy as the "Magnetic" Component of Dark Matter

## 1. Core Idea

In classical electromagnetism, a static charge produces only an electric field $\vec{E}$. A moving charge reveals a second component — the magnetic field $\vec{B}$ — which is not a new entity but a **different projection of the same underlying tensor** $F_{\mu\nu}$.

We propose an analogous structure in the dark sector:

| Electromagnetism | Dark Sector |
|---|---|
| $\vec{E}$ (electric field) | $y_s$ (scalar Yukawa — T-even) |
| $\vec{B}$ (magnetic field) | $y_p$ (pseudoscalar Yukawa — T-odd) |
| $F_{\mu\nu}$ (unified tensor) | $\mathcal{Y} = y_s + iy_p = y\,e^{i\sigma/f}$ (complex Yukawa) |
| Lorentz boost mixes $E \leftrightarrow B$ | $\sigma$ field rotates $y_s \leftrightarrow y_p$ |
| $\vec{B}$ does no work ($\vec{F} \perp \vec{v}$) | $y_p$ doesn't interfere with mass term ($\gamma^5$ flips sign) |
| $B/E$ ratio set by $v/c$ (dynamics) | $y_p/y_s$ ratio set by $\theta = \sigma/f$ (dynamics) |
| Like charges repel electrically, parallel currents attract magnetically (sign flip) | DM attracts via $\phi$ (SIDM), vacuum energy repels via $\sigma$ (DE) — sign flip |

**Thesis**: Dark energy is not a separate phenomenon. It is the **T-violating component** of the dark matter interaction, made manifest through CP violation in the Majorana sector — exactly as magnetism is the T-violating component of electrostatics, made manifest through motion.

---

## 2. The Lagrangian

### 2.1 Starting Point: Secluded Majorana SIDM

The established model (122 viable benchmark points, consistent with relic density, SIDM, and all direct/indirect constraints):

$$\mathcal{L}_{\text{SIDM}} = \frac{1}{2}\bar{\chi}(i\!\!\not\!\partial - m_\chi)\chi + \frac{1}{2}(\partial\phi)^2 - V_0(\phi) - \frac{1}{2}\bar{\chi}(y_s + iy_p\gamma^5)\chi\,\phi$$

where $\chi$ is a Majorana fermion (dark matter), $\phi$ is a scalar mediator, and:

$$V_0(\phi) = \frac{1}{2}m_\phi^2\phi^2 + \frac{\mu_3}{6}\phi^3 + \frac{\lambda_4}{24}\phi^4$$

### 2.2 The Duality Extension

We promote the CP phase to a **dynamical field** $\sigma$:

$$y_s + iy_p \;\longrightarrow\; y\cos\!\left(\frac{\sigma}{f}\right) + iy\sin\!\left(\frac{\sigma}{f}\right)\gamma^5 \;=\; y\,e^{i\gamma^5 \sigma/f}$$

The full Lagrangian becomes:

$$\boxed{\mathcal{L} = \mathcal{L}_{\text{SIDM}}\big|_{y_s \to y\cos(\sigma/f),\; y_p \to y\sin(\sigma/f)} \;+\; \frac{1}{2}(\partial\sigma)^2 - V_{\text{bare}}(\sigma)}$$

Here:
- $y$ = total Yukawa coupling (fixed by SIDM + relic)
- $f$ = decay constant of $\sigma$ (free parameter, expected $f \sim 10^{17}\text{-}10^{18}$ GeV)
- $\sigma$ = pseudo-Goldstone boson of a $U(1)$ that rotates the CP phase
- $V_{\text{bare}}(\sigma)$ = explicit breaking potential (from UV physics)

$\sigma$ is a **dark axion** — an axion-like particle that lives entirely in the dark sector.

### 2.3 Couplings as Functions of $\theta \equiv \sigma/f$

$$\alpha_s(\theta) = \frac{y^2}{4\pi}\cos^2\theta, \qquad \alpha_p(\theta) = \frac{y^2}{4\pi}\sin^2\theta$$

---

## 3. Constraints Fix the Angle

### 3.1 SIDM Constraint

The self-interaction cross section in Born approximation:

$$\frac{\sigma_T}{m_\chi} = \frac{\pi\alpha_s^2}{m_\chi^3} = \frac{\pi\alpha^2}{m_\chi^3}$$

This requires $\alpha_s(\theta) = \alpha$, i.e., $y^2\cos^2\theta = 4\pi\alpha$.

### 3.2 Relic Density Constraint

The s-wave annihilation cross section for mixed scalar-pseudoscalar Majorana coupling:

$$\langle\sigma v\rangle_0 = \frac{2\pi\alpha_s\alpha_p}{m_\chi^2}$$

Matching to the pure-scalar relic formula gives:

$$\alpha_s\alpha_p = \frac{\alpha^2}{8}$$

### 3.3 The Universal Angle

Combining SIDM ($\alpha_s = \alpha$) and relic ($\alpha_s\alpha_p = \alpha^2/8$):

$$\alpha_p = \frac{\alpha}{8}, \qquad \frac{\alpha_p}{\alpha_s} = \frac{1}{8} = \tan^2\theta$$

$$\boxed{\theta_{\text{relic}} = \arctan\!\left(\frac{1}{\sqrt{8}}\right) = 19.47° \approx 0.3398 \text{ rad}}$$

This angle is **universal** — independent of $m_\chi$, $m_\phi$, or $\alpha$. It is a pure number from the algebra of Majorana fermions.

The total Yukawa is then:

$$y = \frac{\sqrt{4\pi\alpha}}{\cos\theta_{\text{relic}}} = \frac{\sqrt{4\pi\alpha}}{\sqrt{8/9}} = \frac{3\sqrt{4\pi\alpha}}{2\sqrt{2}}$$

### 3.4 Physical Meaning

| $\theta$ | CP structure | SIDM | Relic | DE |
|---|---|---|---|---|
| $0°$ | Pure scalar | ✓ | ✗ (p-wave only) | None |
| $19.47°$ | Mixed (our Universe) | ✓ | ✓ | $V(\sigma) \neq 0$ |
| $45°$ | CP-democratic | ✗ ($\alpha_s$ too small) | ✗ | Maximal |
| $90°$ | Pure pseudoscalar | ✗ (no t-channel) | ✗ | Maximal |

**Only $\theta = 19.47°$ gives a viable Universe with both SIDM and the correct relic density. This is the Universe we live in.**

---

## 4. The Vacuum Energy

### 4.1 Coleman-Weinberg Potential

At 1-loop, $\chi$ running in a loop generates a $\theta$-dependent potential:

$$V_{\text{CW}}(\theta) = -\frac{1}{32\pi^2}M_{\text{eff}}^4(\theta)\left[\ln\frac{M_{\text{eff}}^2(\theta)}{m_\chi^2} - \frac{3}{2}\right]$$

where the effective mass-squared (at mediator field value $\phi$):

$$M_{\text{eff}}^2(\theta, \phi) = m_\chi^2 + m_\chi y\cos\theta\,\phi + \frac{y^2\phi^2}{4}$$

**Critical observation**: the $\cos\theta$ dependence means:
- At $\theta = 0$: $M_{\text{eff}}^2$ is **maximized** → $V_{\text{CW}}$ is **most negative** (deepest well)
- At $\theta = 90°$: $M_{\text{eff}}^2$ loses the linear term → $V_{\text{CW}}$ is **least negative** (shallowest)

The "magnetic" component ($y_p$, via $\gamma^5$) **doesn't interfere with the mass term** — just as $\vec{B}$ does no work. This is what creates the sign structure.

### 4.2 The Sign Structure — Duality at Work

| Component | Effect on $V_{\text{CW}}$ | Macroscopic manifestation |
|---|---|---|
| "Electric" ($y_s = y\cos\theta$) | Deepens the potential (more negative) | **Attraction** → DM clustering (SIDM) |
| "Magnetic" ($y_p = y\sin\theta$) | Raises the potential (less negative) | **Repulsion** → accelerated expansion (DE) |

This mirrors EM exactly:
- Like electric charges **repel** → positive energy
- Parallel currents (moving charges) **attract** → this is the magnetic sign flip

In the dark sector: the "magnetic" (T-breaking) component contributes **positive** vacuum energy, i.e., **repulsion on cosmological scales**.

### 4.3 Scale of the Vacuum Energy

The $\sigma$-dependent part of $V_{\text{CW}}$, with Planck-suppressed coupling ($f \sim M_{\text{Pl}}$):

$$V_\sigma \sim \frac{y^4 m_\chi^6}{32\pi^2 f^2}$$

Numerical results with $f = 0.2\,M_{\text{Pl}}$:

| Benchmark | $m_\chi$ (MeV) | $\alpha$ | $V_\sigma / \rho_\Lambda$ |
|---|---|---|---|
| BP1 | 20.7 | $1.05 \times 10^{-3}$ | $10^{-6.6}$ |
| BP9 | 42.5 | $2.17 \times 10^{-3}$ | $10^{-4.0}$ |
| MAP | 94.1 | $5.73 \times 10^{-3}$ | $10^{-1.1}$ |
| MAP_relic | 85.8 | $5.52 \times 10^{-3}$ | $10^{-1.4}$ |

**For MAP and MAP_relic, $V_\sigma$ is within one order of magnitude of $\rho_\Lambda$ with $f \sim 0.2\,M_{\text{Pl}}$.**

The required $f$ for exact match:

| Benchmark | $f$ for $V_\sigma = \rho_\Lambda$ | $f / M_{\text{Pl}}$ |
|---|---|---|
| BP1 | $1.3 \times 10^{15}$ GeV | $5.3 \times 10^{-4}$ |
| BP9 | $2.3 \times 10^{16}$ GeV | $9.5 \times 10^{-3}$ |
| MAP | $6.6 \times 10^{17}$ GeV | 0.27 |
| MAP_relic | $4.9 \times 10^{17}$ GeV | 0.20 |

These are **GUT/string scales** — physically motivated, not arbitrary.

---

## 5. Radiative Stability

### 5.1 Shift Symmetry Protection

$\sigma$ is a pseudo-Goldstone boson of $U(1)_{\text{CP}}: \sigma \to \sigma + c$.

This symmetry is only broken by $V_{\text{bare}}(\sigma)$ and by the CW potential (which is periodic in $\sigma/f$). The leading 1-loop correction to $m_\sigma$ with $M = M_{\text{Pl}}$:

$$\delta m_\sigma^2 = \frac{3m_\chi^2}{8\pi^2 M_{\text{Pl}}^2} \sim (10^{-12} \text{ eV})^2$$

This is much larger than $H_0^2 \sim (10^{-33} \text{ eV})^2$, which means **$m_\sigma$ is NOT naturally at $H_0$** with Planck-scale suppression alone. However, $m_\sigma$ is radiatively stable in the sense that it doesn't receive quadratically divergent corrections — only logarithmic ones controlled by $m_\chi/f$.

### 5.2 The Mass Hierarchy

For quintessence, we need $m_\sigma \lesssim H_0$. The CW-generated mass is:

$$m_\sigma^2 = \frac{1}{f^2}\frac{d^2 V_{\text{CW}}}{d\theta^2} \sim \frac{m_\chi^4}{16\pi^2 f^2}$$

With $f = 0.2\,M_{\text{Pl}}$ and $m_\chi = 90$ MeV:

$$m_\sigma \sim \frac{m_\chi^2}{4\pi f} \sim \frac{(0.09)^2}{12.6 \times 4.9\times10^{17}} \sim 10^{-21} \text{ GeV} \sim 10^{-12} \text{ eV}$$

This is $\sim 10^{21}$ times larger than $H_0$. **Open problem**: either $f$ must be much larger, or an additional mechanism (e.g., alignment, monodromy) must suppress $m_\sigma$ further.

---

## 6. The Minimum Problem

### 6.1 CW Drives to $\theta = 0$

The Coleman-Weinberg potential has its minimum at $\theta = 0$ (pure scalar coupling, no CP violation). This is because $M_{\text{eff}}^2$ is maximized there, making $V_{\text{CW}}$ most negative.

**But $\theta = 0$ gives a Universe with no relic density matching** (p-wave only annihilation). The physical Universe requires $\theta = 19.47°$.

### 6.2 Two Interpretations

**Interpretation A — UV Stabilization (conventional):**
$V_{\text{bare}}(\sigma)$ from string compactification or extra dimensions provides a potential with a minimum at $\theta_{\text{relic}}$. The location of the minimum is a prediction of the UV theory, not fine-tuning.

**Interpretation B — Dynamical (the duality view):**

In EM, the ratio $B/E$ is not set by a potential — it is set by **dynamics** (the velocity of the charge). Similarly, $\theta$ may not be set by $V_{\text{bare}}$ but by the **cosmological dynamics** of the dark sector:

- In the early Universe, $\sigma$ rolls and oscillates
- Dark matter freeze-out **selects** $\theta_{\text{relic}}$ — only the right angle gives the observed $\Omega_{\text{DM}}$
- $\sigma$ gets trapped at $\theta_{\text{relic}}$ by Hubble friction
- The residual vacuum energy at this angle **is** the dark energy

This is analogous to the **misalignment mechanism** for axion dark matter, but here the misalignment produces dark energy rather than dark matter.

### 6.3 CW Energy Barrier at the Wrong Angle

$$\Delta V_{\text{CW}}(\theta=0 \to \theta_{\text{relic}}) \sim 10^{-10} \text{ GeV}^4 \sim 10^{37} \times \rho_\Lambda$$

This means $V_{\text{bare}}$ must overcome a barrier $10^{37}$ times larger than $\rho_\Lambda$. The **difference** between $V_{\text{bare}}$ and $V_{\text{CW}}$ at $\theta_{\text{relic}}$ must equal $\rho_\Lambda$ — this is the CC fine-tuning problem, reappearing in a new form.

---

## 7. Observable Consequences

### 7.1 Fifth Force in the Dark Sector

$\sigma$ mediates a long-range force between $\chi$ particles:

$$V_5(r) = -\frac{g_{\sigma\chi\chi}^2}{4\pi}\frac{e^{-m_\sigma r}}{r}, \qquad g_{\sigma\chi\chi} \sim \frac{m_\chi}{f}$$

The fifth-force parameter:

$$\beta = \frac{g_{\sigma\chi\chi}\,M_{\text{Pl}}}{m_\chi} = \frac{M_{\text{Pl}}}{f} \approx 5 \quad (f = 0.2\,M_{\text{Pl}})$$

**This is NOT a side-effect of the theory. In the duality picture, this IS dark energy** — expressed as a particle-level interaction rather than a cosmological constant.

Testable via:
- CMB (integrated Sachs-Wolfe effect)
- BAO (modified growth rate of structure)
- Weak lensing (DM clustering on Hubble scales)
- Bullet Cluster (DM-DM interaction beyond SIDM)

### 7.2 Equation of State

$w = -1 + \epsilon$ where $\epsilon$ depends on the slow-roll dynamics. With CW potential alone, $\epsilon$ is negligibly small. A more realistic $V_{\text{bare}}$ could give detectable $w \neq -1$.

### 7.3 Redshift-Dependent SIDM

If $\sigma$ evolves even slightly over cosmic time:

$$\frac{\Delta\alpha_s}{\alpha_s} = -2\tan\theta_{\text{relic}} \cdot \Delta\theta$$

This would manifest as different SIDM cross sections at different redshifts — testable by comparing dwarf galaxies ($z \sim 0$) with cluster mergers ($z \sim 0.5\text{-}1$).

---

## 8. Benchmark Points

| Parameter | BP1 | BP9 | MAP | MAP_relic |
|---|---|---|---|---|
| $m_\chi$ (MeV) | 20.69 | 42.53 | 94.07 | 85.84 |
| $m_\phi$ (MeV) | 11.34 | 10.92 | 11.10 | 15.35 |
| $\alpha$ | $1.048\times10^{-3}$ | $2.165\times10^{-3}$ | $5.734\times10^{-3}$ | $5.523\times10^{-3}$ |
| $y$ (total Yukawa) | 0.1217 | 0.1749 | 0.2847 | 0.2794 |
| $\theta_{\text{relic}}$ | 19.47° | 19.47° | 19.47° | 19.47° |
| $\alpha_s$ | $1.048\times10^{-3}$ | $2.165\times10^{-3}$ | $5.734\times10^{-3}$ | $5.523\times10^{-3}$ |
| $\alpha_p$ | $1.310\times10^{-4}$ | $2.706\times10^{-4}$ | $7.168\times10^{-4}$ | $6.904\times10^{-4}$ |
| $f$ for $V_\sigma = \rho_\Lambda$ | $1.3\times10^{15}$ | $2.3\times10^{16}$ | $6.6\times10^{17}$ | $4.9\times10^{17}$ |
| $\beta$ (fifth force) | 5.0 | 5.0 | 5.0 | 5.0 |

---

## 9. Summary of Status

### What works ✓
- Lagrangian is fully consistent with SIDM (122 BPs) and relic density
- Universal angle $\theta = 19.47°$ emerges from Majorana algebra alone
- The EM duality analogy is structurally exact: $y_s + iy_p\gamma^5 \leftrightarrow E + iB$
- Sign structure correct: "magnetic" (T-breaking) component raises $V$ → repulsion → DE
- Scale $\rho_\Lambda$ appears naturally for MAP-like BPs with $f \sim 0.2\,M_{\text{Pl}}$ (GUT/string scale)
- Shift symmetry protects $m_\sigma$ from quadratic divergences
- Fifth force $\beta \approx 5$ is testable

### Open problems ✗
- CW minimum is at $\theta = 0$, not $\theta_{\text{relic}}$ — needs $V_{\text{bare}}$ or dynamical mechanism
- $m_\sigma$ from CW is $\sim 10^{-12}$ eV, not $H_0 \sim 10^{-33}$ eV — hierarchy of 21 orders
- CC problem re-emerges: $V_{\text{bare}}$ cancellation against $V_{\text{CW}}$ must leave residual $\sim \rho_\Lambda$
- SM sector contributions to vacuum energy untouched

### Next steps
1. **Dynamical trapping**: Can cosmological evolution of $\sigma$ during freeze-out trap it at $\theta_{\text{relic}}$?
2. **Monodromy / alignment**: Can $m_\sigma$ be further suppressed by known mechanisms?
3. **Observable predictions**: Compute $\beta$-dependent modifications to CMB, BAO, lensing
4. **UV completion**: What string/GUT model gives $f \sim 0.2\,M_{\text{Pl}}$ with the right $V_{\text{bare}}$?

---

## 10. Comparison to Existing Frameworks

| Framework | Similarity | Difference |
|---|---|---|
| Quintessence | Slowly rolling scalar drives DE | Ours: $\sigma$ is a dark axion, not generic scalar |
| Axion DE (Frieman+ 1995) | Axion with $f \sim M_{\text{Pl}}$ drives DE | Ours: axion rotates CP phase of DM coupling specifically |
| Coupled DE (Amendola 2000) | DE field coupled to DM | Ours: coupling is not ad hoc — it IS the Yukawa structure |
| MaVaN (Fardon+ 2004) | DE from neutrino mass variation | Ours: from Majorana DM, not neutrinos |
| Duality in DM (novel) | — | EM-like duality structure in Yukawa sector is new |

---

*First draft: March 26, 2026*
*Based on numerical analysis in `sigma_radiative_stability.py` and `dark_axion_full.py`*
