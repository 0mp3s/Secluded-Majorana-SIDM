# Self-Interacting Dark Matter from a Majorana Fermion with Light Scalar Mediator:
# Velocity-Dependent Cross Sections via Yukawa Partial-Wave Analysis

**Omer P.**  
*Independent researcher*

**March 2026 — V10 Draft**

---

## Abstract

We construct a minimal **secluded** dark matter model in which a Majorana fermion $\chi$ couples to a light real scalar $\phi$ through a Yukawa interaction $-(y/2)\bar{\chi}\chi\,\phi$. The dark sector has no tree-level coupling to the Standard Model. The $\phi$-mediated self-interaction produces velocity-dependent scattering governed by the attractive Yukawa potential $V(r) = -\alpha\,e^{-m_\phi r}/r$ with $\alpha = y^2/(4\pi)$.

We compute the self-interaction transfer cross section using the **Variable Phase Method (VPM)** — a full partial-wave analysis that captures quasi-bound-state resonances of the Yukawa potential. A systematic scan over $(m_\chi, m_\phi, \alpha)$ — with $m_\chi \in [0.1, 100]$ GeV, $m_\phi \in [0.1, 200]$ MeV — identifies **80,142 parameter points** satisfying the SIDM sweet-spot criteria $\sigma/m(30\text{ km/s}) \in [1, 10]\text{ cm}^2/\text{g}$ and $\sigma/m(1000\text{ km/s}) < 0.1\text{ cm}^2/\text{g}$.

The annihilation $\chi\chi \to \phi\phi$ is **s-wave**, yielding a thermal relic density $\Omega h^2 = 0.120$ for $\alpha \sim 5 \times 10^{-4}$–$5 \times 10^{-3}$ — squarely within the SIDM-viable region. A dedicated cosmological scan using an exact numerical Boltzmann solver identifies **17 benchmark points** forming an "island of viability" at $m_\chi \in [10, 100]$ GeV, $m_\phi \in [7.6, 14.8]$ MeV that simultaneously satisfy relic density, SIDM, and cluster constraints. The primary benchmark has $m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV, $\alpha = 1.05 \times 10^{-3}$, $\Omega h^2 = 0.120$, $\sigma/m(30) = 0.52$ cm$^2$/g, and $\sigma/m(1000) = 0.072$ cm$^2$/g. A quantitative $\chi^2$ fit to 13 astrophysical systems spanning $v = 12$–$4700$ km/s yields $\chi^2/\nu = 0.38$ for the best relic-constrained point and $\chi^2/\nu = 0.12$ for the unconstrained best fit, with all pulls below $1.2\sigma$. A Bayesian MCMC posterior analysis confirms that all 17 relic benchmarks lie within the 95% credible region, with marginalized constraints $m_\chi = 73^{+43}_{-56}$ GeV, $m_\phi = 10.8^{+4.8}_{-4.3}$ MeV (68% CI). We show that a Higgs portal coupling to the SM is generically **excluded** for light mediators ($m_\phi \lesssim \mathcal{O}(\text{GeV})$) by the tension between direct detection ($\sigma_{\rm SI} \propto 1/m_\phi^4$) and BBN lifetime constraints; a secluded dark sector is the natural resolution. The scalar has only 1 degree of freedom; since $m_\phi \gg T_\nu$ at BBN and CMB, $\phi$ is non-relativistic and Boltzmann-suppressed, contributing $\Delta N_{\rm eff} \approx 0$ — trivially consistent with all $N_{\rm eff}$ bounds.

---

## 1. Introduction

The standard cold dark matter (CDM) paradigm, while remarkably successful at explaining large-scale structure, faces persistent challenges at galactic and sub-galactic scales. These include the core-cusp problem [1, 2], the missing satellites problem [3], the too-big-to-fail problem [4], and the diversity of rotation curves [5]. Self-interacting dark matter (SIDM), first proposed by Spergel & Steinhardt [6], offers a compelling resolution: DM particles scatter off each other at short range, thermalizing the inner halo and producing flat density cores.

The key phenomenological requirement is velocity-dependent self-interaction:
$$\sigma/m \sim 1\text{–}10 \text{ cm}^2/\text{g at } v \sim 30 \text{ km/s (dwarfs)}$$
$$\sigma/m \lesssim 0.1 \text{ cm}^2/\text{g at } v \sim 1000 \text{ km/s (clusters)} \quad [7, 8]$$

This velocity dependence arises naturally from Yukawa-type interactions mediated by a light boson [9, 10].

In this work we consider a **secluded dark sector** containing a Majorana fermion $\chi$ coupled to a **real scalar mediator** $\phi$ via a Yukawa interaction $-(y/2)\bar{\chi}\chi\,\phi$. We show that the commonly-invoked Higgs portal connection to the SM ($\lambda_{H\phi}|H|^2\phi^2$) is incompatible with light mediators due to the $1/m_\phi^4$ enhancement of the direct detection cross section (§5.1), and instead adopt a secluded framework with no direct SM coupling. This minimal setup has several attractive features:

1. The DM–DM self-interaction is an **attractive Yukawa potential** $V(r) = -\alpha\,e^{-m_\phi r}/r$, supporting velocity-dependent resonant scattering.
2. The annihilation $\chi\chi \to \phi\phi$ is **s-wave**, enabling a correct thermal relic density with couplings in the SIDM-viable range.
3. The secluded dark sector has $\sigma_{\rm SI} = 0$ — **no tension with direct detection experiments** — and annihilation products remain dark, evading CMB energy-injection bounds.
4. The model is **minimal and anomaly-free** — only three parameters $(m_\chi, m_\phi, \alpha)$, no gauge symmetry, no anomaly cancellation needed.

We compute the self-interaction cross sections using the **Variable Phase Method (VPM)** — a full partial-wave analysis that correctly captures resonance structure near quasi-bound states of the Yukawa potential. This is essential in the regime $\lambda = \alpha m_\chi / m_\phi \sim 1$–$30$ relevant to our parameter space.

**Related work.** Light-mediator SIDM was introduced by Feng, Kaplinghat, Tu & Yu [17] in the context of hidden charged dark matter, and developed systematically by Tulin, Yu & Zurek [9] who provided fitting formulas for the Yukawa transfer cross section. The comprehensive review by Tulin & Yu [10] covers the full landscape of SIDM models. Our work differs from these precedents in several respects: (i) we use the full VPM partial-wave analysis rather than the approximate Hulthén-potential fitting formulas of [9] (also adopted in recent automated pipelines such as sidmkit [24]), which fail to capture resonance structure at $\lambda \gtrsim 1$ (Appendix C); (ii) we demonstrate quantitatively that the Higgs portal is generically excluded for light mediators (§5.1), motivating a secluded dark sector rather than the visible-decay models assumed in much of the earlier literature; (iii) we perform a combined SIDM + relic density scan with an exact Boltzmann solver (§4.4, §6.2), identifying a well-defined island of viability rather than individual benchmark points; and (iv) we provide a quantitative $\chi^2$ fit to 13 astrophysical systems (§4.6) spanning four decades in velocity.

**Physical picture.** The central insight of this work is that quantum mechanics alone — applied to the simplest possible dark sector — is sufficient to explain the observed pattern of dark matter self-interactions from dwarf galaxies to galaxy clusters. A Majorana fermion exchanging a light scalar is nothing more than a Yukawa scattering problem, the same physics that governs the deuteron. The velocity dependence required by observations is not engineered; it is an inevitable consequence of the transition between the Born and resonant regimes of the Yukawa potential as the de Broglie wavelength of the dark matter particle crosses the mediator range. At low velocities (dwarf galaxies), the scattering probes quasi-bound-state resonances that enhance $\sigma/m$; at high velocities (galaxy clusters), the interaction becomes short-ranged and weak. The thermal relic constraint independently selects coupling values that place the dark matter squarely in this transition region, creating a well-defined "island of viability" where the cosmological abundance and the self-interaction phenomenology are simultaneously determined — yielding a one-parameter family of benchmark points, each with a fully determined $\sigma/m(v)$ curve.

Two structural conclusions follow. First, the Higgs portal to the SM is not merely unnecessary — it is *forbidden* by the combined requirements of light-mediator SIDM and direct detection bounds (§5.1). The dark sector must be secluded. This transforms the persistent null results from direct detection experiments from a puzzle into a prediction: the dark matter was never going to scatter off nuclei, because the required coupling channel does not exist. Second, the Majorana nature of the fermion imprints a quantum-statistical signature on the scattering cross section (§3.2, §7.7): the spin-weighted partial-wave sum produces a non-monotonic Majorana-to-Dirac ratio $R(v)$ with a regime ($v \approx 40$–95 km/s) where identical-particle scattering *exceeds* distinguishable-particle scattering — a direct manifestation of the Pauli principle in the dark sector that no Dirac model can replicate.

---

## 2. The Model (The Scalar Mediator)

To resolve the small-scale structure anomalies without violating theoretical unitarity or parity constraints, we consider a simplified dark sector consisting of a Majorana fermion dark matter candidate, $\chi$, interacting via a real scalar mediator, $\phi$. The interaction Lagrangian is given by:

$$\mathcal{L}_{\text{int}}=-\frac{y}{2}\bar{\chi}\chi\phi$$

In the non-relativistic limit, this scalar exchange generates a universally attractive Yukawa potential between the dark matter particles, irrespective of their spin configuration:

$$V(r)=-\frac{\alpha}{r}e^{-m_\phi r}$$

where $\alpha=y^2/(4\pi)$. Because $\chi$ is a Majorana fermion — identical to its own antiparticle — the DM–DM scattering amplitude must be antisymmetrized, leading to non-trivial spin-weighted partial-wave sums (even-$\ell$ singlet weight 1, odd-$\ell$ triplet weight 3) and an overall identical-particle factor of $1/2$; the full derivation is given in §3.2. A natural connection between $\phi$ and the Standard Model would be the Higgs portal coupling $\lambda_{H\phi}|H|^2\phi^2$, which induces a mixing angle $\sin\theta$ and enables both $\phi$ decay and DM–nucleon scattering. We investigated this possibility and show in \S5.1 that it is **generically excluded** for light mediators ($m_\phi \lesssim \mathcal{O}(\text{GeV})$): the $1/m_\phi^4$ enhancement of $\sigma_{\rm SI}$ closes the window between direct detection and BBN bounds by a factor of $\sim 3 \times 10^4$. The dark sector is therefore **secluded** — $\phi$ has no tree-level coupling to the SM.

### 2.1 Parameters

The model has three parameters:

| Parameter | Range scanned / used | Physical meaning |
|-----------|---------------------|------------------|
| $m_\chi$ | 0.1–100 GeV | DM mass |
| $m_\phi$ | 0.1–200 MeV | Mediator mass |
| $\alpha = y^2/(4\pi)$ | $10^{-6}$–$5 \times 10^{-3}$ | Dark fine-structure constant |

All three enter the SIDM cross section and relic density. The model contains no free parameters related to the SM–dark sector coupling.

---

## 3. Self-Interaction Cross Section

### 3.1 Variable Phase Method

The DM–DM scattering cross section via $\phi$ exchange is governed by the attractive Yukawa potential:

$$V(r) = -\alpha\, \frac{e^{-m_\phi r}}{r}$$

We compute the transfer cross section using the **Variable Phase Method (VPM)**, which solves for the phase shifts $\delta_l$ via the ODE [11]:

$$\frac{d\delta_l}{dx} = +\frac{\lambda\, e^{-x}}{\kappa\, x}\left[\hat{j}_l(\kappa x)\cos\delta_l - \hat{n}_l(\kappa x)\sin\delta_l\right]^2$$

where:
- $x = m_\phi r$ (dimensionless radial coordinate)
- $\kappa = k/m_\phi = \mu v_{\rm rel}/m_\phi$ with $\mu = m_\chi/2$ the reduced mass
- $\lambda = \alpha\, m_\chi / m_\phi$ (coupling strength parameter)
- $\hat{j}_l, \hat{n}_l$ are Riccati–Bessel functions

The positive sign corresponds to an **attractive** potential.

### 3.2 Cross Section for Identical Majorana Fermions

For identical Majorana fermions, quantum statistics modifies the partial-wave sum:

$$\sigma_T = \frac{2\pi}{k^2}\left[\sum_{l\,\text{even}} (2l+1)\sin^2\delta_l + 3\sum_{l\,\text{odd}} (2l+1)\sin^2\delta_l\right]$$

The factor of 3 for odd partial waves arises from spin statistics: the Majorana fermion has spin-1/2, and the triplet spin state (weight 3) combines with odd orbital angular momentum to form the antisymmetric total wave function. The factor $2\pi/k^2$ (rather than $4\pi/k^2$) includes the conventional factor of $1/2$ for identical-particle scattering rates.

*Note:* For identical Majorana fermions this formula gives both the elastic **and** the momentum-transfer cross section $\sigma_T^{\rm tr} \equiv \int(1-\cos\theta)\,d\sigma/d\Omega\,d\Omega$; the two are **identically equal**. The proof is simple: spin averaging separates the scattering into a singlet channel (even $l$ only) and a triplet channel (odd $l$ only). Within each channel, $|f(\theta)|^2$ is symmetric under $\theta \to \pi - \theta$ (because $P_l(\cos(\pi-\theta)) = (-1)^l P_l(\cos\theta)$ and all contributing $l$ share the same parity), so $\cos\theta\,|f|^2$ is antisymmetric and its angular integral vanishes. Hence $\int \cos\theta\,d\sigma/d\Omega\,d\Omega = 0$ and $\sigma_T^{\rm tr} = \sigma_{\rm el}$ exactly. We verify this numerically in Appendix C: the ratio $\sigma_T^{\rm tr}/\sigma_{\rm el}$ equals unity to within $10^{-12}$ (machine precision) across both BP1 ($\lambda = 1.9$) and MAP ($\lambda = 48.6$) at all nine SIDM-relevant velocities ($v = 12$–$4700$ km/s), using 200-point Gauss–Legendre angular integration of the full symmetrized differential cross section.

### 3.3 Numerical Implementation

The VPM ODE is integrated using a 4th-order Runge-Kutta scheme with:
- Adaptive $x_{\rm max}$: 50 (low $\kappa$), 80 (medium), 100 (high $\kappa$)
- Adaptive step count: 4000–12000 steps
- Truncation: $l_{\rm max} = \min(\max(3, \lfloor\kappa\rfloor + 3), 80)$ with early termination when the contribution falls below $10^{-3}$ of the running sum
- Centrifugal barrier: $x_{\rm min} = \max(10^{-5}, l/\kappa)$ for $l > 0$

The implementation is verified against:
1. **Free particle test:** $\delta_l(\lambda = 0) = 0$ for all $l$ — PASS
2. **Born limit:** At $\lambda \ll 1$, VPM agrees with the analytical Born formula to within $< 0.1\%$; at benchmark couplings ($\lambda \sim 2$–$49$), Born overestimates by factors of 3–5 at low velocity (Appendix C.2), confirming that the full VPM is essential for our parameter space.
3. **Scipy ground truth:** Agreement with `scipy.integrate.solve_ivp` (RK45, rtol=$10^{-10}$) to $< 0.001\%$ — PASS

### 3.4 Numerical Error Budget

A systematic error analysis identifies three components, tested on the relic-constrained benchmark BP1 ($\lambda = 1.91$), the marginal benchmark BP17 ($\sigma/m(1000) = 0.099$), and the MCMC best-fit MAP ($\lambda = 48.6$):

| Source | 30 km/s (relic BPs) | 1000 km/s (relic BPs) | 1000 km/s (MAP) | Method |
|--------|--------------------|-----------------------|-----------------|--------|
| **Integrator** (RK4 step-size) | < 0.001% | < 0.001% | < 0.01% | Compare default vs 16 000 steps |
| **Truncation** ($l_{\rm max}$) | < 0.01% | 1.8–3.6% | ~28% | Compare $l_{\rm max}$ vs $l_{\rm max}+20$ |
| **Prescription** ($x_{\rm min}$) | < 0.01% | 10.6–13.1% | ~32% | Vary $x_{\rm min}$ factor by ±20% |
| **Total** (quadrature) | **< 0.01%** | **10.8–13.6%** | **~43%** |  |

At the primary SIDM velocity (30 km/s, dwarf galaxies), the total systematic is negligible (< 0.01%) for all relic-constrained benchmarks. At cluster scales ($v \sim 1000$ km/s), the 10.8–13.6% systematic for relic BPs is dominated by the $x_{\rm min}$ prescription — an inherent sensitivity of the VPM phase-shift integration boundary. For the MAP point ($\lambda = 48.6$, deep resonant regime with many contributing partial waves), the systematic reaches ~43% at 1000 km/s due to large truncation and prescription sensitivities; this does not affect our viability selection since MAP is used only for the observational fit, not the relic cut.

Among the 17 relic benchmarks, two have $\sigma/m(1000)$ near the selection cut: BP16 ($\sigma/m = 0.096$) and BP17 ($\sigma/m = 0.099$). Under the worst-case $x_{\rm min}$ shift (+5.4%), these would reach $\sim$0.10–0.105 cm$^2$/g, marginally above our conservative selection cut of 0.1 cm$^2$/g. However, the observational cluster constraints are substantially weaker: Harvey et al. [21] quote $\sigma/m < 0.47$ cm$^2$/g at 68% CL at 1000 km/s, and our 0.1 cut is a conservative choice rather than a hard observational bound. All 17 benchmarks remain safely below the observational limits even with the full systematic applied.

---

## 4. Parameter Scan

### 4.1 Scan Grid

The scan covers the following grid, centered on the first four Yukawa quasi-bound-state resonances:

| Dimension | Range | Points | Scale |
|-----------|-------|--------|-------|
| $m_\chi$ | 0.1–100 GeV | 50 | log |
| $m_\phi$ | 0.1–200 MeV | 70 | log |

The upper bound $m_\chi \leq 100$ GeV is a practical choice: for heavier DM, the required $\alpha_{\rm relic} \propto m_\chi$ (§6.1) pushes $\lambda = \alpha m_\chi/m_\phi$ well above the fourth Yukawa resonance, where the cross section becomes increasingly oscillatory and the scan grid would need to be significantly finer. Extending to $m_\chi \sim$ TeV is straightforward in principle but computationally expensive and deferred to future work.
| Resonance | $\lambda_{\rm crit} \in \{1.68, 6.45, 14.7, 26.0\}$ | 4 | — |
| $\alpha$ (around $\alpha_{\rm crit}$) | $\pm 30\%$ log | 200 | log |

Total evaluations: $50 \times 70 \times 4 \times 200 = 2{,}800{,}000$.

### 4.2 Selection Criteria

A point $(m_\chi, m_\phi, \alpha)$ is deemed viable if:

$$\sigma/m(30\text{ km/s}) \in [1, 10] \text{ cm}^2/\text{g}$$
$$\sigma/m(1000\text{ km/s}) < 0.1 \text{ cm}^2/\text{g}$$

The initial scan also imposes $m_\phi > 2m_e = 1.022$ MeV, which would ensure $\phi \to e^+e^-$ before nucleosynthesis in a Higgs-portal model. Although this BBN filter is not required in the secluded framework adopted here ($\S$5) — where $\phi$ is stable and contributes negligibly to $\Delta N_{\rm eff}$ ($\S$5.3) — we retain it to demonstrate that viable points exist even under the more stringent assumption. In the cosmological scan (\S4.4), the dwarf-scale lower bound is relaxed to $\sigma/m(30) \geq 0.5$ cm$^2$/g, consistent with the range inferred from recent analyses of diverse rotation curves and halo density profiles [10, 13, 19]. In particular, Kamada et al. [19] demonstrate that cores form at $\sigma/m \gtrsim 0.3$ cm$^2$/g across diverse dwarf galaxy rotation curves, and Kaplinghat et al. [13] fit a wide range of systems with $\sigma/m \sim 0.5$–3 cm$^2$/g.

### 4.3 Results

| Quantity | Count |
|----------|-------|
| Raw viable samples | **80,142** |
| BBN-safe ($m_\phi > 2m_e$) | **50,751** |
| Representative (1 per grid cell) | **646** |
| BBN-safe representative | **394** |

**Parameter ranges of viable region:**

| Parameter | Min | Max |
|-----------|-----|-----|
| $m_\chi$ | 13.9 GeV | 100 GeV |
| $\alpha$ | $1.8 \times 10^{-6}$ | $3.3 \times 10^{-3}$ |
| $\lambda = \alpha m_\chi / m_\phi$ | 1.18 | 28.93 |

### 4.4 Cosmological Benchmark Points

To go beyond the raw SIDM scan (§4.1–4.3), we perform a dedicated cosmological scan over a $20 \times 30$ grid in $(m_\chi, m_\phi)$ with $m_\chi \in [10, 100]$ GeV and $m_\phi \in [1, 50]$ MeV. For **each grid cell**, we bisect in $\alpha$ to find the unique coupling that yields $\Omega h^2 = 0.120$ via an exact numerical Boltzmann solver (4th-order Runge-Kutta; see §6.2), then evaluate the SIDM cross sections with our VPM solver. This "cosmology as compass" approach treats relic density as an **input** rather than a post-hoc filter.

Of 600 converged cells, **17 satisfy both** the relic constraint ($\Omega h^2 = 0.120 \pm 0.001$) and the relaxed SIDM criteria ($\sigma/m(30) \geq 0.5$ cm$^2$/g, $\sigma/m(1000) < 0.1$ cm$^2$/g). The coarse $20 \times 30$ grid may miss narrow viable strips between cells; a finer grid would likely interpolate additional points within the same island but is unlikely to shift its boundaries significantly, since the 17 points already span the full $(m_\chi, m_\phi)$ range continuously. The 0.5 cm$^2$/g threshold represents the conservative end of the observationally-permitted window: it is sufficient to produce $\mathcal{O}$(kpc) cores in dwarfs [19] while remaining consistent with the diversity of rotation curves [13]. The 17 points form a contiguous **"island of viability"** in parameter space (Figure 1):

$$m_\chi \in [10, 100] \text{ GeV}, \quad m_\phi \in [7.6, 14.8] \text{ MeV}, \quad \lambda \equiv \frac{\alpha\,m_\chi}{m_\phi} \in [0.73, 32.4]$$

We select the point with the lowest $\sigma/m(1000)$ (largest safety margin from the cluster bound) as our primary benchmark:

**Table 1: Primary Benchmark Point (BP1)**

| Parameter | Value |
|:----------|:------|
| $m_\chi$ | 20.69 GeV |
| $m_\phi$ | 11.34 MeV |
| $\alpha$ | $1.048 \times 10^{-3}$ |
| $\lambda$ | 1.91 |
| $\Omega h^2$ | 0.1200 (numerical Boltzmann) |
| $\sigma/m(30\text{ km/s})$ | 0.516 cm$^2$/g |
| $\sigma/m(1000\text{ km/s})$ | 0.072 cm$^2$/g |

**Table 2: Range of the 17 viable points**

| Parameter | Min | Max |
|:----------|:----|:----|
| $m_\chi$ [GeV] | 10.0 | 100.0 |
| $m_\phi$ [MeV] | 7.56 | 14.85 |
| $\alpha$ | $5.54 \times 10^{-4}$ | $4.81 \times 10^{-3}$ |
| $\lambda$ | 0.73 | 32.4 |
| $\sigma/m(30)$ [cm$^2$/g] | 0.50 | 0.79 |
| $\sigma/m(1000)$ [cm$^2$/g] | 0.072 | 0.099 |

![Figure 1: Island of viability in the $(m_\chi, m_\phi)$ plane. Left: color scale shows $\sigma/m(30\text{ km/s})$. Right: $\sigma/m(1000\text{ km/s})$. Stars mark the 17 viable points satisfying all constraints.](v31_island_of_viability.png)

![Figure 2: Velocity-dependent cross section $\sigma/m(v)$ for BP1 ($m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV). Shaded bands indicate the dwarf ($v \sim 30$ km/s) and cluster ($v \sim 1000$ km/s) velocity scales.](v31_bp1_velocity_profile.png)

*Note:* All relic densities are obtained from the full numerical Boltzmann equation (§6.2), not the Kolb-Turner analytic approximation — the latter underestimates $\Omega h^2$ by $\sim$25–30% at these masses, which is sufficient to invalidate benchmark selections.

### 4.5 Comparison with Observational Data

To test whether our benchmark points are consistent with astrophysical observations, we compare the predicted $\sigma/m(v)$ curves against data extracted from five published analyses spanning four decades in velocity:

| System | $v$ [km/s] | $\sigma/m$ [cm$^2$/g] | 68% CL range | Reference |
|:-------|:----------:|:---------------------:|:------------:|:---------:|
| Draco dSph | 12 | 0.6 | [0.1, 2.0] | [13] |
| Fornax dSph | 12 | 0.8 | [0.2, 3.0] | [13] |
| NGC 2976 | 60 | 2.0 | [0.5, 5.0] | [13] |
| NGC 1560 | 55 | 3.0 | [1.0, 8.0] | [13] |
| IC 2574 | 50 | 1.5 | [0.3, 5.0] | [13] |
| NGC 720 (group) | 250 | 0.5 | [0.1, 1.5] | [13] |
| NGC 1332 (group) | 280 | 0.3 | [0.05, 1.0] | [13] |
| Abell 611 | 1200 | 0.1 | [0.02, 0.3] | [13] |
| Abell 2537 | 1100 | 0.15 | [0.03, 0.4] | [13] |
| Diverse rotation curves | 40 | 3.0 | [0.5, 10.0] | [19] |
| Bullet Cluster$^*$ | 4700 | — | $< 1.25$ | [20] |
| 72 cluster mergers$^*$ | 1000 | — | $< 0.47$ | [21] |
| TBTF dwarfs | 30 | 1.0 | [0.5, 5.0] | [22] |

**Table 3: Observational constraints on $\sigma/m$ used for comparison.** The nine entries from Kaplinghat et al. [13] (Draco through Abell 2537) are derived from SIDM halo profile fits to observed rotation curves and velocity dispersions — these are model-dependent inferences that assume an SIDM+NFW profile, not direct measurements of the scattering cross section. The four remaining entries (Diverse RC [19], Bullet Cluster [20], Harvey et al. [21], TBTF [22]) are independent of the SIDM halo model. Entries marked $^*$ are one-sided upper limits: they contribute $\chi^2 = 0$ when the theory prediction lies below the limit, and are penalized only when exceeded.

BP1 falls within the allowed range of **11 out of 13** observational systems (Figure 3). The two marginal cases — NGC 2976 ($\sigma/m_{\rm theory} = 0.50$ vs lower bound 0.50 cm$^2$/g, a 0.4% shortfall) and NGC 1560 ($\sigma/m_{\rm theory} = 0.50$ vs lower bound 1.0 cm$^2$/g) — are both rotation-curve systems where the inferred $\sigma/m$ has large uncertainties ($\sim$0.5 dex) from halo-profile modeling. The overall agreement across dwarfs ($v \sim 10$–50 km/s), galaxies ($v \sim 50$–300 km/s), and clusters ($v \sim 1000$–5000 km/s) demonstrates that the velocity dependence of our model naturally reproduces the astrophysical pattern: $\sigma/m \sim 0.5$ cm$^2$/g in dwarfs, declining through the group and cluster scales to $\sigma/m \sim 0.05$ cm$^2$/g at $v \sim 1000$ km/s.

![Figure 3: Predicted $\sigma/m(v)$ for BP1 (solid blue), BP5 (dashed orange), and BP17 (dotted green) compared with observational constraints from Kaplinghat et al. [13], Kamada et al. [19], Randall et al. [20], Harvey et al. [21], and Elbert et al. [22]. Error bars show 68% CL ranges. All three benchmarks trace a similar velocity-dependent envelope consistent with the observed pattern.](v33_observational_comparison.png)

### 4.6 Quantitative $\chi^2$ Fit to Observational Data

To move beyond the qualitative compatibility assessment of §4.5, we perform a full $\chi^2$ fit of the predicted $\sigma/m(v)$ to all 13 observational systems listed in Table 3. For each parameter point we evaluate the VPM transfer cross section at all 12 unique velocities and compute:

$$\chi^2 = \sum_{i=1}^{N_{\rm eff}} \left(\frac{\sigma/m_{\rm theory}(v_i) - \sigma/m_{\rm obs,\,i}}{\sigma_i}\right)^2$$

with asymmetric errors: $\sigma_i^+ = (\sigma/m_{\rm upper} - \sigma/m_{\rm obs})/1.0$ and $\sigma_i^- = (\sigma/m_{\rm obs} - \sigma/m_{\rm lower})/1.0$ at 68% CL, selecting $\sigma_i^+$ or $\sigma_i^-$ according to the sign of the residual. The Bullet Cluster [20] and Harvey et al. [21] constraints are **one-sided upper limits**: they contribute $\chi^2 = 0$ when the theory prediction satisfies $\sigma/m < \sigma/m_{\rm limit}$, and are only penalized when exceeded. Since all our viable points satisfy these bounds, $N_{\rm eff} = 11$ active constraints. With 3 free parameters, the number of degrees of freedom is $\nu = 10$.

We evaluate 5,009 points sampled from the 80,142 raw viable set plus all 17 relic benchmarks, totaling 5,026 points $\times$ 12 velocities = 60,312 full VPM evaluations. A fine-grained scan around the top 50 candidates adds 400 further evaluations.

**Unconstrained best fit.** The global minimum is:

$$\chi^2_{\rm min}/\nu = 1.17/10 = 0.12$$

at $m_\chi = 100$ GeV, $m_\phi = 9.15$ MeV, $\alpha = 2.84 \times 10^{-3}$, $\lambda = 31.0$. All 11 active pulls are below $1\sigma$ (largest: NGC 1560 at $-0.61\sigma$); the two upper-limit systems (Bullet, Harvey) contribute zero. Representative predictions:

| System | $v$ [km/s] | Theory | Obs | Pull |
|:-------|:----------:|:------:|:---:|:----:|
| Draco dSph | 12 | 1.32 | 0.60 | $+0.52$ |
| NGC 1560 | 55 | 1.77 | 3.00 | $-0.61$ |
| Diverse RC | 40 | 1.95 | 3.00 | $-0.42$ |
| Abell 2537 | 1100 | 0.11 | 0.15 | $-0.33$ |
| Bullet Cluster$^*$ | 4700 | 0.008 | $< 1.25$ | — |
| TBTF dwarfs | 30 | 1.91 | 1.00 | $+0.23$ |

The low $\chi^2/\nu \ll 1$ reflects the generous observational uncertainties characteristic of astrophysical SIDM constraints — typical 68% CL ranges span $\sim$0.5 dex (a factor of 3), compared with the $\mathcal{O}(10\%)$ precision common in collider physics. In this regime, $\chi^2/\nu < 1$ does not imply over-fitting but rather that the data do not yet have the resolving power to distinguish between models at the $\sim$factor-of-2 level. A more discriminating test will require tighter observational error bars, e.g., from JWST-era kinematics of ultra-faint dwarfs. Nevertheless, the fit demonstrates that a single set of parameters can **simultaneously** match all 13 systems spanning $v = 12$–$4700$ km/s — a non-trivial consistency check given the four-decade velocity baseline.

**Circularity note.** Nine of the 13 observational entries (Draco through Abell 2537, all from Kaplinghat et al. [13]) are derived from SIDM halo-profile fits — their extraction assumes a cored DM density profile, which is itself an outcome of SIDM. These data therefore provide a consistency check rather than an independent test. The remaining four entries (Kamada et al. [19], Randall et al. [20], Harvey et al. [21], Elbert et al. [22]) are derived from methods independent of the SIDM halo model (rotation-curve diversity, merger kinematics, lensing, and subhalo counts). Restricting to the four independent constraints alone, all benchmarks remain compatible: the Bullet and Harvey upper limits are satisfied, the TBTF and diverse-RC central values lie within the predicted range.

**Relic-constrained best fit.** Restricting to the 17 benchmark points satisfying $\Omega h^2 = 0.120$ (§4.4):

$$\chi^2_{\rm relic}/\nu = 3.82/10 = 0.38$$

at the $(m_\chi, m_\phi)$ grid cell nearest to BP1: $m_\chi = 20.69$ GeV, $m_\phi = 9.91$ MeV, $\alpha = 1.048 \times 10^{-3}$, $\lambda = 2.19$ — hereafter **BP1$_\chi$** to distinguish it from the relic-selected BP1 of Table 1 ($m_\phi = 11.34$ MeV). The two points share the same $m_\chi$ but differ by one grid step in $m_\phi$; BP1 (Table 1) appears as rank 17 in Table 4 with $\chi^2/\nu = 0.68$, still a good fit. Throughout this paper, "BP1" without subscript always refers to the Table 1 definition ($m_\phi = 11.34$ MeV). The worst pull is NGC 1560 at $-1.12\sigma$, well within acceptable limits.

**Table 4: All 17 relic-constrained benchmark points ranked by $\chi^2$**

| Rank | $\chi^2/\nu$ | $m_\chi$ [GeV] | $m_\phi$ [MeV] | $\alpha$ | $\sigma/m(30)$ [cm$^2$/g] |
|:----:|:------------:|:--------------:|:--------------:|:--------:|:-------------------------:|
| 1 | 0.38 | 20.69 | 9.91 | $1.05 \times 10^{-3}$ | 0.79 |
| 2 | 0.41 | 29.76 | 11.34 | $1.47 \times 10^{-3}$ | 0.75 |
| 3 | 0.42 | 14.38 | 8.66 | $7.56 \times 10^{-4}$ | 0.73 |
| 4 | 0.47 | 100.00 | 14.85 | $4.81 \times 10^{-3}$ | 0.59 |
| 5 | 0.47 | 26.37 | 11.34 | $1.31 \times 10^{-3}$ | 0.68 |
| 6 | 0.48 | 18.33 | 9.91 | $9.39 \times 10^{-4}$ | 0.68 |
| 7 | 0.49 | 42.81 | 12.98 | $2.09 \times 10^{-3}$ | 0.64 |
| 8 | 0.51 | 88.59 | 14.85 | $4.26 \times 10^{-3}$ | 0.53 |
| 9 | 0.54 | 37.93 | 12.98 | $1.86 \times 10^{-3}$ | 0.60 |
| 10 | 0.55 | 78.48 | 14.85 | $3.78 \times 10^{-3}$ | 0.50 |
| 11 | 0.56 | 23.36 | 11.34 | $1.17 \times 10^{-3}$ | 0.60 |
| 12 | 0.57 | 10.00 | 7.56 | $5.54 \times 10^{-4}$ | 0.59 |
| 13 | 0.58 | 12.74 | 8.66 | $6.80 \times 10^{-4}$ | 0.59 |
| 14 | 0.59 | 33.60 | 12.98 | $1.65 \times 10^{-3}$ | 0.56 |
| 15 | 0.62 | 16.24 | 9.91 | $8.41 \times 10^{-4}$ | 0.56 |
| 16 | 0.66 | 29.76 | 12.98 | $1.47 \times 10^{-3}$ | 0.51 |
| 17 | 0.68 | 20.69 | 11.34 | $1.05 \times 10^{-3}$ | 0.52 |

All 17 relic benchmarks yield $\chi^2/\nu < 1$, confirming that the relic-density constraint does **not** degrade the observational fit. The full island of viability is quantitatively consistent with all available astrophysical data.

**Effect of velocity averaging.** The above $\chi^2$ analysis evaluates $\sigma/m$ at the single characteristic velocity $v_{\rm char}$ of each system, following the standard practice in SIDM phenomenology [9, 10, 13]. In reality, DM particles in a halo have a Maxwell–Boltzmann (MB) velocity distribution, and the observable is the thermally averaged $\langle\sigma/m\rangle_{\rm MB} = \int (\sigma/m)(v)\,v\,f_{\rm MB}(v;v_0)\,dv\,/\,\int v\,f_{\rm MB}(v;v_0)\,dv$ with $v_0 = v_{\rm char}/\sqrt{2}$. We computed this integral for all 17 BPs via Gauss–Legendre quadrature ($n=30$, see Appendix D). The per-system deviations are largest at group/cluster velocities ($v \sim 250$–$1500$ km/s), where the steep decline of $\sigma/m(v)$ causes the MB tail to sample lower cross sections, yielding $\langle\sigma/m\rangle_{\rm MB} \approx 0.83$–$0.88\times\sigma/m(v_{\rm char})$ — a $\sim$12–17% reduction. At dwarf scales ($v \lesssim 60$ km/s) the correction is $\lesssim 5\%$. The net effect on the $\chi^2$ fit is a uniform increase of $2$–$11\%$ (mean $7\%$): no benchmark point changes viability status, and the worst $\chi^2/\nu$ rises from 0.52 to 0.55. This is entirely negligible compared to the $\sim$0.5 dex observational uncertainties.

![Figure 4: Best-fit $\sigma/m(v)$ curves with 13 observational data points and 68% CL error bars. The relic-constrained best fit (BP1$_\chi$, $\chi^2/\nu = 0.38$) reproduces the velocity-dependent pattern from dwarfs to clusters.](v34_chi2_fit.png)

### 4.7 Bayesian Posterior Constraints

To quantify the allowed parameter space beyond the frequentist $\chi^2$ analysis, we perform a Bayesian posterior sampling of the three-dimensional parameter space $(m_\chi, m_\phi, \alpha)$ using the emcee affine-invariant ensemble sampler [23]. We adopt flat (uninformative) priors in $\log_{10}$ space:

$$\log_{10}(m_\chi/\text{GeV}) \in [\log_{10}(5),\, \log_{10}(200)], \quad \log_{10}(m_\phi/\text{MeV}) \in [\log_{10}(3),\, \log_{10}(30)], \quad \log_{10}\alpha \in [-5,\, \log_{10}(0.05)]$$

The likelihood is Gaussian: $\ln\mathcal{L} = -\chi^2/2$ with the same $\chi^2$ function used in §4.6 (asymmetric errors, 13 systems including one-sided upper limits). The sampler uses 32 walkers initialized around the 17 relic benchmark points plus the unconstrained best fit of §4.6, with a Gaussian scatter of $\sigma = 0.05$ in log-space. After 300 burn-in steps, we collect 5,000 production steps (160,000 total samples, computed in parallel on 12 CPU cores).

**Convergence diagnostics.** The mean acceptance fraction is 0.542, and the maximum integrated autocorrelation time is $\tau_{\rm max} = 75.1$ steps, yielding $N_{\rm eff} \approx 2{,}132$ effective independent samples. The ratio $N_{\rm steps}/\tau_{\rm max} = 66.6 > 50$ confirms convergence [Foreman-Mackey+ 2013]. The chain trace plots (Figure S1) show good mixing with no residual drift. Figure 7 shows the integrated autocorrelation time $\hat{\tau}_{\rm int}$ as a function of chain length: all three parameters plateau well before the end of the chain, and the autocorrelation function decays to zero within $\sim$200 lag steps, confirming that the posterior is fully equilibrated.

**Selection methodology.** We emphasize that no post-hoc parameter tuning enters the analysis. The 17 relic benchmark points are identified by the cosmological scan (§4.4) — selected solely by the objective criterion $\Omega h^2 = 0.120 \pm 0.001$ combined with the SIDM viability window — before any $\chi^2$ comparison with observations. The MAP estimate emerges as the maximum of the full three-dimensional posterior density, not from manual inspection of individual points. The broad posterior (68\% CI spanning $\sim$1 dex in each parameter) further demonstrates that the good fit is not confined to an isolated fine-tuned point but extends across a wide parameter region.

**Results.** The maximum a posteriori (MAP) estimate is:
$$m_\chi = 94.1~\text{GeV},\quad m_\phi = 11.1~\text{MeV},\quad \alpha = 5.7 \times 10^{-3},\quad \chi^2/\nu = 0.20$$

**Important:** The MAP point is **not relic-constrained** — its tree-level relic density is $\Omega h^2 \approx 0.076$, approximately $37\%$ below the observed value. The MAP emerges as the maximum of the observational $\chi^2$ posterior, which does not include a relic-density prior. The 17 relic benchmark points (which do satisfy $\Omega h^2 = 0.120$) all lie within the 95% credible region, confirming that the observational and cosmological constraints are compatible. If one imposes a relic-density prior, the posterior MAP shifts to the relic-constrained region (§4.4), with a somewhat worse $\chi^2/\nu \approx 0.38$ (BP1$_\chi$).

The marginalized parameter constraints (median $\pm$ 68% credible interval) are:

| Parameter | Median | 16th %ile | 84th %ile |
|-----------|--------|-----------|-----------|
| $m_\chi$ [GeV] | 35.7 | 10.1 | 88.1 |
| $m_\phi$ [MeV] | 8.4 | 5.0 | 12.7 |
| $\alpha$ | $1 \times 10^{-3}$ | $< 10^{-4}$ | $4 \times 10^{-3}$ |
| $\lambda = \alpha m_\chi/m_\phi$ | 5.5 | 0.8 | 29.0 |

The broad credible intervals reflect the large astrophysical uncertainties ($\sim$0.5 dex) in the observational data — the model fits the data well across a wide parameter range. The derived coupling parameter $\lambda$ spans both the Born ($\lambda < 1$) and resonant ($\lambda > 1$) regimes, with a median $\lambda = 5.5$ indicating moderate resonant enhancement.

Crucially, **all 17 relic benchmark points** (§4.4) lie within the 95% credible region of the posterior, with $\chi^2$ values ranging from 5.4 to 8.4 (compared to the 95th percentile threshold of the posterior $\chi^2$ distribution). This confirms that the island of viability identified by the cosmological scan coincides with the statistically preferred region of the SIDM parameter space.

![Figure 5: Corner plot showing the 2D marginalized posterior distributions of $\log_{10}(m_\chi)$, $\log_{10}(m_\phi)$, and $\log_{10}\alpha$, with 68% and 95% contour levels. Red lines mark the MAP estimate. All 17 relic BPs fall within the 95% contours.](v38_corner.png)

![Figure 6: Posterior distribution of the derived parameter $\lambda = \alpha m_\chi/m_\phi$. The median is $\lambda = 5.5$ with a broad 68% CI of [0.8, 29]. The distribution spans both the Born ($\lambda < 1$) and resonant ($\lambda > 1$) regimes.](v38_lambda_posterior.png)

![Figure 7: MCMC convergence diagnostics. Left: integrated autocorrelation time $\hat{\tau}_{\rm int}$ vs chain length for each parameter. All three converge well below the $N/50$ threshold (dashed). Right: autocorrelation function (walker-averaged) showing decay to zero within $\sim$200 steps.](v38_autocorr_diagnostic.png)

---

## 5. Secluded Dark Sector and Direct Detection

### 5.1 Incompatibility of the Higgs Portal with Light Mediators

The most economical connection between the dark scalar $\phi$ and the SM is the Higgs portal $\lambda_{H\phi}|H|^2\phi^2$, the lowest-dimension gauge-invariant operator coupling a singlet scalar to the SM. We initially considered this coupling as part of the model, as it would provide both a $\phi$ decay channel (ensuring BBN safety) and a testable direct-detection signal. However, a quantitative analysis reveals that the portal is **intrinsically excluded** for the light mediators required by SIDM.

The mixing angle $\sin\theta$ induced after electroweak symmetry breaking mediates DM–nucleon scattering with a spin-independent cross section:

$$\sigma_{\rm SI} = \frac{y^2\,\sin^2\theta\,f_N^2\,\mu_N^2\,m_N^2}{\pi\,v_{\rm EW}^2\,m_\phi^4}$$

where $f_N \approx 0.3$ is the nucleon form factor, $\mu_N$ is the DM–nucleon reduced mass, and $v_{\rm EW} = 246$ GeV. The crucial feature is the **$1/m_\phi^4$ dependence**: for $m_\phi = 11$ MeV, this produces an enhancement of $(m_h/m_\phi)^4 \approx 1.5 \times 10^{16}$ relative to standard Higgs-mediated scattering.

For our primary benchmark BP1 ($m_\chi = 20.69$ GeV, $\alpha = 1.048 \times 10^{-3}$), the LZ bound ($\sigma_{\rm SI} < 4 \times 10^{-47}$ cm$^2$ at $m_\chi \approx 20$ GeV) requires:

$$\sin\theta < 6 \times 10^{-10}$$

Meanwhile, BBN demands that $\phi$ decays before $T \sim 1$ MeV. For $\phi \to e^+e^-$ ($m_\phi > 2m_e$), this requires:

$$\sin\theta > 2 \times 10^{-5}$$

These two bounds are **mutually exclusive by a factor of $\sim 3 \times 10^4$**:

$$\sin\theta_{\rm BBN} \approx 2 \times 10^{-5} \gg \sin\theta_{\rm LZ} \approx 6 \times 10^{-10}$$

This tension is **generic** for any light-mediator SIDM model with $m_\phi \lesssim \mathcal{O}(\text{GeV})$ and a Higgs portal coupling. The $1/m_\phi^4$ enhancement is so severe that no region of the $(m_\phi, \sin\theta)$ plane simultaneously satisfies BBN and direct detection bounds. We therefore adopt a **secluded dark sector** framework [17, 18].

### 5.2 The Secluded Model

We set $\lambda_{H\phi} = 0$ (equivalently, $\sin\theta = 0$), removing any tree-level coupling between $\phi$ and the Standard Model. The dark sector contains only $\chi$ and $\phi$, interacting via the Yukawa coupling $y$.

**Consequences:**
- $\sigma_{\rm SI} = \sigma_{\rm SD} = 0$ exactly — **no direct detection signal**
- $\phi$ is **stable** — no decay to SM particles
- The CMB energy injection constraint (§6.3) does not apply — annihilation products remain in the dark sector
- The model is fully described by **three parameters**: $(m_\chi, m_\phi, \alpha)$

**Thermal history:** We assume the dark sector was in thermal equilibrium with the SM at early times ($T \gg m_\chi$). A concrete realization is a heavy scalar or fermionic mediator $\Sigma$ with mass $M_\Sigma \gg m_\chi$ and couplings to both sectors (e.g., $\Sigma\bar{f}f + \Sigma\bar{\chi}\chi$). For $M_\Sigma \sim 1$–$10$ TeV, the thermalization rate $\Gamma \sim T^5/M_\Sigma^4$ exceeds the Hubble rate for $T \gtrsim \mathcal{O}(10)$ GeV, establishing equilibrium well before freeze-out at $T_f \sim m_\chi/25 \sim 0.4$–$4$ GeV. At $T \ll M_\Sigma$ the interaction decouples and is negligible at the scales relevant to SIDM and freeze-out. After decoupling, the dark-sector temperature tracks $T_{\rm dark} = \xi\, T_{\rm SM}$ with $\xi \lesssim 1$. More precisely, $\xi$ evolves as a power law in the entropy degrees of freedom: $\xi(T) = \xi_{\rm dec}\,[g_{*S}(T_{\rm dec})/g_{*S}(T)]^{1/3}$, where $T_{\rm dec}$ is the decoupling temperature [Farina et al. 2016, eq. 9]. For decoupling above the QCD transition ($T_{\rm dec} \gtrsim \text{few GeV}$, as implied by the heavy-mediator thermalization scenario), $g_{*S}$ changes by a factor $\sim$3 between decoupling and recombination, giving $\xi_0/\xi_{\rm dec} \approx 0.5$–$0.8$. Since $\Delta N_{\rm eff} \propto \xi^4$, adopting $\xi \approx 1$ (as we do throughout) is conservative — the true contribution is smaller. The standard Boltzmann calculation applies without modification.

### 5.3 Cosmological Safety of Stable $\phi$

A stable $\phi$ with $m_\phi \sim 8$–15 MeV would contribute to the radiation energy density only if relativistic ($T_\phi \gg m_\phi$). In the hypothetical massless limit, a single real scalar degree of freedom gives:

$$\Delta N_{\rm eff}^{\rm massless} = \frac{4}{7}\left(\frac{T_\phi}{T_\nu}\right)^4 \approx 0.027$$

for $\xi = 1$ and $\phi$ decoupling before $e^+e^-$ annihilation. However, for our massive $\phi$, the actual contribution is Boltzmann-suppressed: at BBN ($T_{\rm SM} \sim 1$ MeV), $m_\phi/T_\phi \sim 22$ gives a suppression factor $\sim e^{-22} \approx 3 \times 10^{-10}$; at CMB recombination ($T_{\rm SM} \sim 0.26$ eV), $m_\phi/T_\phi \sim 10^8$ and $\Delta N_{\rm eff}$ is effectively zero. A systematic sensitivity scan over the full parameter space ($m_\phi \in [0.1, 50]$ MeV, $m_\chi \in [5, 200]$ GeV) confirms that $\Delta N_{\rm eff}$ is monotonically decreasing in $m_\phi$ within the viable region, with a maximum of $\Delta N_{\rm eff} \sim 7 \times 10^{-12}$ at the lightest viable corner ($m_\phi = 7.5$ MeV, $m_\chi = 10$ GeV) — a safety margin of $> 10^{10}$ relative to Planck and $> 8 \times 10^{9}$ relative to CMB-S4. CMB-S4 could detect $\Delta N_{\rm eff} > 0.06$ only for $m_\phi \lesssim 0.1$ MeV, two orders of magnitude below the viable range. The model is therefore trivially consistent with both Planck 2018 ($\Delta N_{\rm eff} < 0.30$ at 95% CL) and the projected CMB-S4 sensitivity ($\sigma \sim 0.03$).

The energy density in stable $\phi$ particles at late times is:

$$\frac{\Omega_\phi}{\Omega_\chi} \sim \frac{m_\phi}{m_\chi} \sim 5 \times 10^{-4}$$

— a negligible correction to the total DM energy budget. This estimate assumes that the $\phi$ number density tracks thermal equilibrium, which requires a number-changing process. In the secluded dark sector, the leading candidate is the cubic self-coupling $(\mu_3/3!)\,\phi^3$ enabling $3\phi \leftrightarrow 2\phi$ ("cannibal") annihilation. A coupled Boltzmann analysis [Farina et al. 2016] shows that for $\mu_3/m_\phi \gtrsim 1.7$, this process is efficient and $\Omega_\phi$ remains subdominant; for smaller $\mu_3/m_\phi$, the $\phi$ abundance is not sufficiently depleted and the dark sector is overclosed. We restrict to the former regime. We acknowledge that $\mu_3$ constitutes a fourth parameter of the model, beyond the $(m_\chi, m_\phi, \alpha)$ triad emphasized in this work; however, it enters only the cosmological history of $\phi$ (not the SIDM cross sections or the DM relic density), and any value $\mu_3/m_\phi \gtrsim 1.7$ yields identical phenomenology. In this sense the cubic coupling sets a necessary condition for cosmological consistency rather than introducing a continuous degree of freedom into the SIDM predictions.

### 5.4 Falsifiability

Without a direct detection signal, the model is tested through:
1. **Astrophysical observations** of dwarf galaxy cores, cluster mergers, and halo density profiles — the primary SIDM observable
2. **CMB-S4** sensitivity to $\Delta N_{\rm eff} \sim 0.03$ cannot probe our massive $\phi$ (Boltzmann-suppressed, $\Delta N_{\rm eff} \approx 0$), but a null result would be consistent with the model
3. **Tightening cluster bounds** on $\sigma/m(1000)$ below 0.01 cm$^2$/g would further constrain the island of viability
4. **Relic density:** if $\alpha$ is measured independently from halo observations, the relic prediction is fixed — a non-trivial consistency test

---

## 6. Relic Density

### 6.1 Annihilation Cross Section

DM freezes out through $\chi\chi \to \phi\phi$ (t/u-channel, scalar mediator). For a Majorana fermion with Yukawa coupling $y$, the **s-wave** cross section is:

$$\langle\sigma v\rangle = \frac{y^4}{64\pi\,m_\chi^2} = \frac{\pi\alpha^2}{4\,m_\chi^2}$$

This is the leading, velocity-independent term. The s-wave nature has two important consequences:
- The relic density is set efficiently without velocity enhancement — the same $\alpha$ that determines SIDM also determines the relic abundance.
- Late-time annihilation (CMB, Fermi-LAT) is **not** velocity-suppressed; however, in the secluded model these constraints do not apply (see §6.3).

### 6.2 Numerical Boltzmann Solver

We solve the Boltzmann equation numerically [12]:
$$\frac{dY}{dx} = -\sqrt{\frac{\pi}{45}}\,\frac{g_{*S}}{\sqrt{g_*}}\,M_{\rm Pl}\,m_\chi\,\frac{\langle\sigma v\rangle_0}{x^2}\left(Y^2 - Y_{\rm eq}^2\right)$$

with $\langle\sigma v\rangle_0 = \pi\alpha^2/(4 m_\chi^2)$ (constant, s-wave). The effective degrees of freedom $g_*(T)$ and $g_{*S}(T)$ are tabulated following [16].

| $m_\chi$ [GeV] | $\alpha_{\rm relic}$ | $\Omega h^2$ | SIDM viable? |
|:---:|:---:|:---:|:---:|
| 10.0 | $5.5 \times 10^{-4}$ | 0.120 | Yes ($\sigma/m_{30}=0.59$) |
| 20.7 | $1.05 \times 10^{-3}$ | 0.120 | **Yes — BP1** |
| 42.8 | $2.09 \times 10^{-3}$ | 0.120 | Yes ($\sigma/m_{30}=0.64$) |
| 100 | $4.81 \times 10^{-3}$ | 0.120 | Yes ($\sigma/m_{30}=0.59$) |

**Key result:** For each mass in the viable island ($m_\chi \in [10, 100]$ GeV, $m_\phi \in [7.6, 14.8]$ MeV), there exists a unique $\alpha$ that yields $\Omega h^2 = 0.120$ and simultaneously satisfies the SIDM sweet-spot criteria. The relic-density coupling $\alpha \sim 5 \times 10^{-4}$–$5 \times 10^{-3}$ lies **squarely within** the SIDM-viable parameter space. This represents a significant advantage over the axial-vector mediator approach, where the p-wave annihilation required $\alpha \sim 10^{-3}$–$10^{-2}$ — largely outside the SIDM range.

### 6.3 CMB and Indirect Detection Constraints

S-wave annihilation is not velocity-suppressed at late times. In models with a Higgs portal, this would place stringent constraints via CMB energy injection ($p_{\rm ann}$) and Fermi-LAT gamma-ray searches. However, in the **secluded model** adopted here (§5.2), the annihilation proceeds as $\chi\chi \to \phi\phi$, with the $\phi$ particles remaining in the dark sector. **No energy is injected into the SM baryon-photon plasma**, so:

$$f_{\rm eff} = 0 \quad \Rightarrow \quad p_{\rm ann} = 0$$

The Planck CMB constraint on $p_{\rm ann}$ **does not apply**. Similarly, Fermi-LAT dwarf spheroidal searches for DM annihilation require SM final states (photons, leptons, hadrons). Since $\chi\chi \to \phi\phi$ produces only dark-sector particles, the Fermi-LAT constraint is automatically satisfied.

**Consequence:** The low-mass end of the viable island ($m_\chi \gtrsim 10$ GeV) is unconstrained by CMB or indirect detection, in contrast to Higgs-portal models where $m_\chi \lesssim 15$ GeV is marginally excluded. This opens the full parameter space identified in §4.4.

### 6.4 Bound-State Formation (BSF)

For $\lambda = \alpha m_\chi / m_\phi \sim 1$–$29$, bound states can form via $\chi\chi \to (\chi\chi)_{\rm bound} + \phi$ [14, 15]. The Coulomb estimate gives $\sigma v_{\rm BSF} / \langle\sigma v\rangle_{\rm ann} \sim \alpha^3/v^2 \sim 10^{-12}$ at freeze-out. BSF is negligible for the relic density calculation. At late times ($v \sim 30$ km/s), BSF is an **inelastic** process and does not modify the elastic $\sigma_T/m$ relevant for SIDM.
### 6.5 Sommerfeld Enhancement

The same attractive Yukawa potential $V(r) = -\alpha\,e^{-m_\phi r}/r$ that mediates DM self-interactions also acts on the initial state of the annihilation $\chi\chi \to \phi\phi$. This distorts the two-particle wavefunction at the annihilation point, multiplying the tree-level cross section by the Sommerfeld factor $S_0 = |\psi_k(0)|^2/|\psi_k^{\rm free}(0)|^2$. For s-wave annihilation only the $l=0$ partial wave contributes.

We compute $S_0$ numerically by solving the $l=0$ radial Schr\"odinger equation for the Yukawa potential in dimensionless units $x = m_\phi r$:
$$u''(x) + \left[\kappa^2 + \lambda\,\frac{e^{-x}}{x}\right]u(x) = 0$$
with $u(0)=0$ and $u'(0)=1$, where $\kappa = \mu v/m_\phi$ and $\lambda = \alpha m_\chi/m_\phi$. The enhancement factor is extracted from the asymptotic amplitude: $S_0 = 1/(\kappa^2 A^2)$ with $A^2 = u^2 + (u'/\kappa)^2$ evaluated in the potential-free region ($x \gg 1$). For high momenta ($\kappa > 10$) the Coulomb analytic formula provides a tight upper bound.

Results for all 17 relic-viable benchmark points:

| Velocity regime | $\kappa$ range | $S_0$ range | Impact |
|:---------------|:--------------:|:-----------:|:-------|
| Freeze-out ($v_{\rm fo} \sim 0.3c$) | 50–760 | 1.003–1.025 | $\Delta\Omega h^2/\Omega h^2 < 2.5\%$ |
| Clusters ($v \sim 1000$ km/s) | 0.5–9 | 1.5–8 | No constraint (secluded) |
| MW-size ($v \sim 200$ km/s) | 0.1–2 | 3.6–45 | No constraint (secluded) |
| Dwarfs ($v \sim 30$ km/s) | 0.01–0.3 | 5.7–417 | $\tau_{\rm depl} \sim 10^7$–$10^8\,t_{\rm Hubble}$ |

**Freeze-out:** The maximum correction to the relic density is $|\Delta\Omega h^2/\Omega h^2| < 2.5\%$ (for BP11 with $\lambda = 32.4$), well within the $\sim$5\% uncertainty from $g_*(T)$ tabulations and the QCD transition. The tree-level Boltzmann solver of \S6.2 is therefore valid across the entire viable island.

**Late-time annihilation:** At halo velocities ($v \sim 30$ km/s), $S_0$ can reach $\sim 400$ near resonances. However, this is irrelevant for two independent reasons: (i) in the secluded model, $\chi\chi \to \phi\phi$ injects no energy into the SM plasma, so CMB and indirect-detection constraints do not apply (\S6.3); and (ii) even with $S_0 \sim 400$, the depletion timescale $\tau_{\rm depl} = m_\chi/(\rho_{\rm DM}\,S\langle\sigma v\rangle_{\rm tree}) \sim 10^7$–$10^8\,t_{\rm Hubble}$, negligible compared to the age of the universe.
---

## 7. Phenomenological Predictions

The preceding sections established the model's theoretical consistency — viable parameter space (§3–4), secluded dark sector (§5), relic density (§6). We now confront the model with six astrophysical phenomena spanning dwarf ($v \sim 10$–$60$ km/s) to cluster ($v \sim 1000$ km/s) scales.

Throughout this section we employ three benchmark points from the MCMC posterior (§4.7):

- **BP1** ($m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV, $\alpha = 1.048 \times 10^{-3}$, $\lambda = 1.91$) — Born-to-classical transition, collider-accessible mass range.
- **BP9** ($m_\chi = 37.93$ GeV, $m_\phi = 12.98$ MeV, $\alpha = 1.858 \times 10^{-3}$, $\lambda = 5.43$) — first-resonance regime.
- **MAP** ($m_\chi = 94.07$ GeV, $m_\phi = 11.10$ MeV, $\alpha = 5.734 \times 10^{-3}$, $\lambda = 48.6$) — resonant regime, maximum a posteriori estimate from the astrophysical $\chi^2$ posterior (§4.7). *Note:* unlike BP1 and BP9, the MAP emerges from the full three-dimensional posterior with $(m_\chi, m_\phi, \alpha)$ as free parameters and is **not** relic-constrained. Running the Boltzmann solver of §6 on the MAP parameters yields $\Omega h^2 = 0.076$ — the coupling $\alpha$ preferred by astrophysical data overshoots the relic value by $\sim$25%, producing an under-abundant relic ($\Omega_{\rm MAP}/\Omega_{\rm obs} = 0.64$). The MAP should therefore be interpreted as the astrophysical best-fit illustrating the deep-resonant regime; achieving the correct relic density at this mass would require $\alpha_{\rm relic} \approx 4.5 \times 10^{-3}$ ($\lambda \approx 38$), still deep in the resonant regime with qualitatively identical SIDM phenomenology.

For the CP-separation analysis (§7.1), we generalize the Lagrangian to include both scalar and pseudoscalar couplings: $\mathcal{L} \supset \frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\,\phi$, with $\alpha_{s,p} = y_{s,p}^2/(4\pi)$. The s-wave annihilation cross section becomes $\langle\sigma v\rangle_0 = 2\pi\alpha_s\alpha_p/m_\chi^2$, which fixes the product $\alpha_s \alpha_p$ while self-interactions depend only on $\alpha_s$ through the non-relativistic Yukawa potential. Equating this to the pure-scalar result $\pi\alpha^2/(4m_\chi^2)$ used in §6 gives $\alpha_s \alpha_p = \alpha^2/8$; thus the single coupling $\alpha$ of §2–§6 maps onto the SIDM-relevant coupling $\alpha_s = \alpha$ while $\alpha_p = \alpha/8$ at the CP-symmetric point ($\alpha_s = \alpha_p = \alpha/\sqrt{8}$ gives an equivalent parametrisation). This matching ensures identical $\langle\sigma v\rangle_0$ — and therefore identical relic densities — between the two descriptions. Unlike Dirac fermion models, where a single Yukawa coupling must simultaneously satisfy both relic density and self-interaction constraints — creating tension between the two — the Majorana structure naturally separates these requirements: $y_p$ controls $s$-wave annihilation while $y_s$ governs the Yukawa self-interaction potential, eliminating internal fine-tuning.

### 7.1 CP-Violating Structure of the Coupling Space

Our mixed Majorana Lagrangian $\frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\,\phi$ introduces two independent Yukawa couplings, scalar ($y_s$) and pseudoscalar ($y_p$), parametrized through $\alpha_{s,p} = y_{s,p}^2/(4\pi)$. The relic density constraint $\langle\sigma v\rangle_0 = 2\pi\alpha_s\alpha_p/m_\chi^2$ fixes the product $\alpha_s \times \alpha_p$, while self-interactions depend only on $\alpha_s$ through the Yukawa potential $V(r) = -\alpha_s e^{-m_\phi r}/r$. This creates a **CP-separation band**: a continuous family of viable models parametrized by $\alpha_s/\alpha_p$, ranging from the CP-symmetric point ($\alpha_s = \alpha_p$) to highly CP-violating configurations ($\alpha_s \gg \alpha_p$ or vice versa).

To map this band, we fix the mass spectrum at each benchmark point and scan $\alpha_s$ while enforcing three constraints: (i) the relic product $\alpha_s \alpha_p = \alpha^2/8 = 1.387 \times 10^{-7}$ (for BP1, $\alpha = 1.048 \times 10^{-3}$), derived by equating the pure-scalar (§6) and mixed (§7.1) annihilation formulae, (ii) dwarf-scale SIDM: $\sigma/m(30\text{ km/s}) \in [0.1, 10]$ cm$^2$/g, and (iii) cluster safety: $\sigma/m(1000\text{ km/s}) < 1$ cm$^2$/g.

**BP1 masses** ($m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV). Here $\lambda = \alpha_s m_\chi/m_\phi$ ranges from 2.4 to 9.9 across the viable band, straddling the first Born resonance at $\lambda = \pi$. The 13 viable points span $\alpha_s/\alpha_p \in [13, 212]$ — a dynamic range of 1.22 decades. Points with $\lambda < \pi$ exhibit an $s$-wave scattering plateau at low velocities, while those with $\lambda > \pi$ transition to the classical regime with mild velocity dependence $\sigma_T \propto (\ln\lambda)^2/v^2$.

**MAP masses** ($m_\chi = 94.07$ GeV, $m_\phi = 11.10$ MeV). The high mass ratio $m_\chi/m_\phi = 8{,}475$ places even modest $\alpha_s$ values in the resonant regime ($\lambda \gg \pi$). **All 500 scanned points pass viability**, spanning $\alpha_s/\alpha_p \in [1.8, 11{,}532]$ — a dynamic range of **3.81 decades**. The oscillating VPM cross-section maintains $\sigma/m(30) > 0.1$ cm$^2$/g via resonance peaks even at very small $\alpha_s$, while the $\sigma/m \propto 1/v^4$ suppression at cluster velocities keeps $\sigma/m(1000)$ between 0.004 and 0.28 cm$^2$/g across the entire band — safely below the Bullet Cluster limit.

The key conclusion is that **CP violation is a generic prediction** of the model: the relic-SIDM constraints permit coupling asymmetries spanning 1–4 orders of magnitude depending on the mass spectrum, with the MAP region showing the widest viable band. At the MAP benchmark itself ($\alpha_s = 5.734 \times 10^{-3}$, $\alpha_s/\alpha_p = 237$), the coupling ratio exceeds $10^2$, indicating substantial CP violation as the natural state of this parameter space.

The CP-separation band has two important implications: (1) **phenomenological distinguishability** — different points along the band predict identical self-interactions but different annihilation signatures, relic pathways, and (in principle) collider phenomenology in CP-sensitive observables; (2) **robustness** — the model's SIDM predictions are insensitive to the degree of CP violation, depending only on $\alpha_s$ and the mass spectrum. Any future measurement constraining $\alpha_s/\alpha_p$ (e.g., through indirect detection channels sensitive to $s$-wave vs. $p$-wave annihilation) would select a unique point within the band without affecting the astrophysical predictions.

### 7.2 Fornax Globular Cluster Survival

The persistence of five globular clusters (GCs) in the Fornax dSph, despite dynamical friction timescales shorter than a Hubble time in an NFW cusp, constitutes one of the strongest small-scale challenges to $\Lambda$CDM (Tremaine 1976; Goerdt et al. 2006; Read et al. 2006). In SIDM, self-interactions flatten the central cusp into a constant-density core, within which dynamical friction vanishes as the gravitational wake symmetrizes in the homogeneous medium (Goerdt et al. 2006; Read et al. 2006; Kaur & Sridhar 2018). GCs that spiral inward to $r \sim r_{\rm core}$ stall indefinitely.

We model Fornax with $M_{200} = 3.16 \times 10^9~M_\odot$, $\sigma_v = 11.7$ km/s (Walker et al. 2009), and compute the SIDM core size using the Kaplinghat et al. (2016) criterion. For each of Fornax's 5 GCs (masses from Mackey & Gilmore 2003), we estimate the 3D position from the projected radius using three deprojection assumptions (face-on: $r_{\rm 3D} = R_{\rm proj}$; mean: $r_{\rm 3D} = \frac{4}{\pi} R_{\rm proj}$; median: $r_{\rm 3D} = \frac{\pi}{2} R_{\rm proj}$), then classify each as SAFE ($r_{\rm 3D} < r_{\rm core}$ → stalled), INSPIRAL ($r_{\rm 3D} > r_{\rm core}$ and $t_{\rm DF} < t_{\rm age}$), or MARGINAL ($r_{\rm 3D} > r_{\rm core}$ but $t_{\rm DF} > t_{\rm age}$). We adopt a conservative prescription following Goerdt et al. (2006): dynamical friction vanishes inside the constant-density core ($r < r_{\rm core}$), where the isotropic density distribution cancels the gravitational wake. A smooth transition (e.g., Kaur & Sridhar 2018, where $F_{\rm DF} \propto d\ln\rho/d\ln r$) would reduce inspiral rates further, making our survival scores conservative lower bounds.

**Results.** The MAP benchmark ($\alpha = 5.7 \times 10^{-3}$) produces $r_{\rm core} = 829$ pc, enclosing all 5 GCs at face-on and mean deprojections: **14/15 safe** (5 GCs \times 3 deprojections). The single marginal case is Fornax 3 at median deprojection ($r_{\rm 3D} = 860$ pc, just outside $r_{\rm core}$), with inspiral time $t_{\rm DF} \sim 0.6$ Gyr — shorter than $t_{\rm age}$ but comparable to the systematic uncertainty in the deprojection. The core size is consistent with the observationally inferred core radius of Fornax ($r_{\rm core} \approx 500$–$900$ pc; Walker & Pe\~narrubia 2011; Amorisco & Evans 2012).

BP1 ($\alpha = 1.05 \times 10^{-3}$) gives $r_{\rm core} = 449$ pc — marginally sufficient. GC3 at mean deprojection ($r_{\rm 3D} \approx 548$ pc) lies outside the core with a remaining inspiral time of $\sim 0.9$ Gyr to the core edge, yielding a score of **12/15**. BP9 ($\alpha = 1.858 \times 10^{-3}$, $m_\phi = 12.98$ MeV) produces $r_{\rm core} = 541$ pc and scores **13/15** — the larger core (compared to BP1) stalls Fornax 4 at all deprojections, with only Fornax 3 at median deprojection showing a $\sim$5 Gyr inspiral. **The viable parameter space therefore includes points — particularly BP9 and the MAP region — that naturally explain Fornax GC survival without fine-tuning.**

The observed projected positions of the GCs ($R_{\rm proj} = 240$–$1710$ pc) are naturally explained as GCs that formed at $r > r_{\rm core}$, spiraled inward under dynamical friction, and stalled upon entering the constant-density core. No additional mechanism beyond SIDM core formation is required.

**Stellar velocity dispersion profile.** As an independent cross-check, we solve the spherical Jeans equation in the isotropic limit ($\beta = 0$) to predict the line-of-sight stellar velocity dispersion profile $\sigma_{\rm los}(R)$. We model the stellar distribution as a Plummer profile with $r_{\rm half} = 710$ pc, $M_* = 2.0 \times 10^7\;M_\odot$ (Walker et al. 2009), embedded in the SIDM-cored+NFW gravitational potential computed above. The radial velocity dispersion $\sigma_r^2(r) = \rho_*^{-1}\int_r^\infty \rho_*(r')\,G\,M_{\rm tot}(r')/r'^2\,dr'$ is then Abel-projected to $\sigma_{\rm los}(R)$, with **no free parameters** — the benchmark points are fixed by relic density and astrophysical cross sections.

The NFW baseline gives $\sigma_{\rm los}(R=100\;\text{pc}) \approx 16.5$ km/s, overshooting the Walker et al. (2009) measurement of $11.4 \pm 1.0$ km/s by $\sim 5\sigma$ — a clear kinematic signature of the cusp problem. BP1 ($r_{\rm core} = 449$ pc) reduces the central prediction to 12.1 km/s (within 1$\sigma$) and achieves $\chi^2/\text{dof} = 39/8$, a factor 3.9 improvement over the NFW ($\chi^2 = 152$). MAP ($r_{\rm core} = 829$ pc) lowers the central dispersion further to 9.3 km/s — slightly below the data, reflecting over-coring (the same trade-off seen in rotation curves, §7.3). All models overpredict $\sigma_{\rm los}$ at $R > 1$ kpc by $\sim 2$ km/s, a well-known artifact of the isotropic assumption; mild radial anisotropy ($\beta \sim 0.2$–$0.3$; Walker & Peñarrubia 2011) would reconcile the outer profile. The key result — **SIDM coring resolves the central dispersion excess** — is robust against the anisotropy treatment.

**Baryonic feedback.** Supernova-driven outflows in Fornax — which has $M_*/M_{\rm halo} = 6.3 \times 10^{-3}$, near the peak of the Di Cintio et al. (2014) feedback efficiency curve — create an additional dark matter core independently of SIDM. We model this using the coreNFW profile of Read, Agertz & Collins (2016): $M_{\rm cNFW}(<r) = M_{\rm NFW}(<r) \times [\tanh(r/r_c)]^n$, where $r_c = 1.75\, r_{\rm half} = 1.24$ kpc and $n = \tanh(\kappa \, M_*/M_{\rm halo}) = 0.47$ for $\kappa = 80$. Applying SIDM coring on top of the feedback-modified profile, the BP1 core grows from 449 to 662 pc and $\chi^2$ improves from 39.1 to 34.9 (11\% improvement). MAP is essentially unaffected ($\chi^2 = 48.1 \to 48.3$), as its already-large core dominates over the feedback correction.

**Velocity anisotropy.** The isotropic ($\beta = 0$) analysis above yields $\sigma_{\rm los}$ profiles that are too flat at large $R$: the Walker et al. (2009) data decrease from 12.0 to 8.1 km/s between $R = 400$ and 1800 pc, while all SIDM models predict a nearly constant profile. We generalize to the Osipkov-Merritt anisotropy profile $\beta(r) = r^2/(r^2 + r_a^2)$ (Osipkov 1979; Merritt 1985), where $r_a$ is the anisotropy radius: orbits are isotropic for $r \ll r_a$ and radially biased for $r \gg r_a$. We solve the anisotropic Jeans equation $\sigma_r^2(r) = [\rho_* Q(r)]^{-1} \int_r^\infty \rho_*\, Q\, G\, M_{\rm tot} / r'^2\, dr'$ with $Q = 1 + (r/r_a)^2$, and project using the kernel $(1 - \beta R^2/r^2)\, r/\sqrt{r^2-R^2}$ (Mamon & Łokas 2005). Scanning $r_a \in [0.3, 50]$ kpc for each scenario:

- **MAP SIDM+feedback with $r_a = 0.8$ kpc: $\chi^2 = 2.3/7$** — all residuals $< 1\sigma$. The model predicts $\sigma_{\rm los} = [12.2, 12.1, 11.9, 11.8, 11.7, 11.4, 11.0, 10.4, 9.7]$ km/s, reproducing the peaked-then-decreasing shape of the data. The best-fit $r_a \approx r_{\rm half}$ is physically natural: stars within the half-light radius have tangentially-mixed orbits, while those beyond are preferentially on radial plunges.

- BP1 SIDM+feedback prefers $r_a \to \infty$ (isotropic): anisotropy monotonically worsens the fit ($\chi^2 = 200$ at $r_a = 0.3$ kpc to 34.9 at $r_a = 50$ kpc). BP1's smaller core (662 pc) produces a $\sigma_{\rm los}$ profile that is too flat; the Osipkov-Merritt kernel cannot reshape it into the observed falling profile without simultaneously dragging down the central values below the data.

- NFW similarly cannot be rescued by anisotropy ($\chi^2 \geq 152$ for all $r_a$).

The **MAP + feedback + $\beta(r)$ result constitutes the best fit to Fornax kinematics in this work**, demonstrating that at dSph velocity scales ($v \sim 10$–$20$ km/s), the larger cross section ($\sigma/m \approx 1.21$ cm$^2$/g) is required to produce the deep isothermal core that, when combined with standard kinematic modeling, matches the Walker et al. (2009) data to sub-sigma precision. This result is consistent with the velocity-dependent complementarity discussed in §7.6: MAP excels at dSph scales while BP1 excels at rotation-curve scales (§7.3), as expected for a velocity-dependent mediator.

**Multi-dSph cross-validation.** To verify that the Fornax result is not a statistical fluke, we repeat the full analysis pipeline — NFW, coreNFW, SIDM, SIDM+feedback, and Osipkov-Merritt anisotropy scan — on four additional classical dSphs: Sculptor, Draco, Carina, and Sextans, using Walker et al. (2009) binned $\sigma_{\rm los}(R)$ data and halo parameters consistent with Wolf et al. (2010) and Read et al. (2019). The results demonstrate that the velocity-dependent Yukawa cross section naturally accommodates the observed diversity:

- **Fornax** ($M_*/M_{\rm halo} = 6.3 \times 10^{-3}$): MAP SIDM+fb+$\beta$, $\chi^2/\text{dof} = 2.0/8 = 0.25$ ($r_a = 0.8$ kpc) — deep core from $\sigma/m \approx 1.21$ cm$^2$/g.
- **Sculptor** ($M_*/M_{\rm halo} = 1.5 \times 10^{-3}$): BP1 SIDM+fb+$\beta$, $\chi^2/\text{dof} = 12.9/9 = 1.44$ ($r_a = 0.3$ kpc) — moderate core; MAP overcores ($\chi^2 = 152$).
- **Draco** ($M_*/M_{\rm halo} = 3.6 \times 10^{-4}$): NFW+$\beta$ provides the best fit ($\chi^2/\text{dof} = 0.7/7$) — consistent with Draco's cuspy profile (Read et al. 2019) and negligible feedback ($n_{\rm fb} = 0.03$).
- **Carina** ($M_*/M_{\rm halo} = 9.5 \times 10^{-4}$): MAP SIDM+fb+$\beta$, $\chi^2/\text{dof} = 2.6/7 = 0.38$ ($r_a = 0.3$ kpc) — excellent agreement.
- **Sextans** ($M_*/M_{\rm halo} = 4.4 \times 10^{-4}$): MAP SIDM (isotropic), $\chi^2/\text{dof} = 41.4/7 = 5.92$ — the poorest fit, likely reflecting tidal effects on this extended, diffuse system ($r_{\rm half} = 695$ pc).

In total, **4 out of 5 classical dSphs achieve $\chi^2/\text{dof} < 2$**, with no single benchmark point required to fit all systems simultaneously. The diversity itself — MAP for cored dSphs (Fornax, Carina), BP1 for intermediate systems (Sculptor), and cuspy NFW for pristine halos (Draco) — is the expected signature of a velocity-dependent SIDM cross section, where the effective $\sigma/m$ varies with the local $\sigma_v$.

### 7.3 Radial Acceleration Relation

The radial acceleration relation (RAR) discovered by McGaugh et al. (2016) establishes a tight empirical correlation between the observed centripetal acceleration $g_{\rm obs} = V_{\rm obs}^2/r$ and the baryonic acceleration $g_{\rm bar} = V_{\rm bar}^2/r$ in disk galaxies, described by

$$g_{\rm obs} = \frac{g_{\rm bar}}{1 - \exp\left(-\sqrt{g_{\rm bar}/g_\dagger}\right)}$$

with the characteristic scale $g_\dagger = 1.2 \times 10^{-10}$ m/s² and an intrinsic scatter of $\sim 0.13$ dex. Any viable dark matter model must reproduce this relation without excessive scatter.

We test the SIDM prediction against 7 SPARC galaxies spanning two regimes: gas-dominated dwarfs (DDO 154, IC 2574, NGC 2366, UGC 128; $V_{\rm max} = 47$–$66$ km/s) and baryon-rich spirals (NGC 2403, NGC 2976, NGC 3198; $V_{\rm max} = 90$–$150$ km/s). The SPARC database provides a single pre-computed baryonic velocity $V_{\rm bar}^2 = V_{\rm gas}^2 + \Upsilon_{*,\text{def}} V_{\rm disk}^2$ with default $\Upsilon_{*,\text{def}} = 0.5\;M_\odot/L_\odot$. To properly disentangle the gas component (which does not scale with $\Upsilon_*$) from the stellar disk, we employ literature gas mass fractions $f_{\rm gas} = M_{\rm gas}/M_{\rm bar}$ from resolved HI surveys (Oh et al. 2015; de Blok et al. 2008; Lelli et al. 2016; Adams et al. 2014), yielding:

$$V_{\rm tot}^2(r) = f_{\rm gas}\,V_{\rm bar}^2 + \Upsilon_*\cdot 2(1-f_{\rm gas})\,V_{\rm bar}^2 + V_{\rm DM}^2(r)$$

where the factor of 2 undoes the default $\Upsilon_{*,\text{def}} = 0.5$ embedded in $V_{\rm bar}$. This decomposition is exact in the limit $f_{\rm gas} \to 1$ (gas-dominated dwarfs) and becomes a standard global approximation for mixed systems (cf. Di Cintio et al. 2014; Santos-Santos et al. 2018). The DM contribution $V_{\rm DM}$ includes adiabatic contraction and SIDM isothermal coring via the thermalization condition $\rho(r_1)\,(\sigma_T/m)\,v\,t_{\rm age} = 1$.

**BP1 results** ($m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV, $\alpha = 1.048 \times 10^{-3}$). The fitted stellar mass-to-light ratios for spirals — $\Upsilon_* = 0.62$ (NGC 2976), 0.84 (NGC 3198), 0.91 (NGC 2403) $M_\odot/L_\odot$ — are consistent with 3.6 $\mu$m stellar population synthesis expectations ($0.2$–$0.8\;M_\odot/L_\odot$; Meidt et al. 2014). For gas-dominated dwarfs ($f_{\rm gas} > 0.75$), $\Upsilon_*$ is unconstrained as expected: the stellar disk contributes $< 10\%$ of $V_{\rm bar}^2$, and the rotation curve is governed by gas + DM alone. Crucially, BP1 reduces the RAR scatter from 0.283 to 0.179 dex for dwarfs (37% improvement) and from 0.177 to 0.122 dex for spirals (31%), demonstrating that SIDM coring tightens the RAR at both ends of the acceleration spectrum.

**MAP results** ($m_\chi = 94.07$ GeV, $\alpha = 5.734 \times 10^{-3}$, $\lambda = 48.6$). The MAP benchmark produces SIDM cores with $r_1 \sim 1$–$14$ kpc — comparable to or exceeding the optical scale lengths of the sample galaxies. This over-coring drives spiral fits to unphysical $\Upsilon_* > 1.3$ and pushes dwarf fits to the upper bound ($\Upsilon_* = 3.0$), with $\chi^2/\text{dof} > 19$ for NGC 2403 and NGC 3198. The tension provides a rotation-curve-based upper bound on the self-interaction strength at dwarf velocities: $\sigma/m(30\;\text{km/s}) \lesssim \text{a few}$ cm$^2$/g for consistency with SPARC rotation curves.

The contrast between BP1 and MAP illustrates a key feature of velocity-dependent SIDM: a single benchmark cannot simultaneously optimize predictions across all mass scales. BP1 ($\sigma/m \approx 0.5$ cm$^2$/g at 30 km/s) represents the sweet spot for galactic-scale phenomenology, while MAP excels at cluster scales (§7.2) and UFD environments (§7.5). This velocity-dependent complementarity is a structural prediction of the model (§7.6). Our results are consistent with the SIDM rotation-curve fits of Ren et al. [23], who demonstrated that velocity-dependent cross sections from a light-mediator model can simultaneously explain the diversity and uniformity of SPARC rotation curves. Our approach goes further by deriving the cross section from first principles — the Lagrangian parameters $(m_\chi, m_\phi, \alpha)$ are fixed by the relic density constraint, leaving no free phenomenological parameters in the SIDM prediction.

### 7.4 Supermassive Black Hole Seeds

An intriguing possibility is that SIDM can seed supermassive black holes (SMBHs) at high redshift through gravothermal collapse of the central halo (Pollack, Spergel & Steinhardt 2015; Feng, Yu & Zhang 2021). The early-universe SMBHs discovered by JWST — with $M_{\rm BH} \sim 10^6$–$10^8~M_\odot$ at $z \gtrsim 10$ (Harikane et al. 2023; Maiolino et al. 2024) — challenge conventional formation channels and have motivated SIDM-based explanations.

We investigate whether gravothermal collapse can operate within our model's viable parameter space. Using the Correa et al. (2015) concentration-mass-redshift relation for halos at $z = 6$–$20$ with $M_{200} = 10^9$–$10^{11}~M_\odot$, we compute the gravothermal collapse timescale $t_{\rm gc} \approx 150 \times t_{\rm relax}$ (Balberg, Shapiro & Inagaki 2002), where $t_{\rm relax} = [\rho_s \,(\sigma/m)\, \sigma_v]^{-1}$.

**The result is unambiguously negative.** For all three benchmark points and across the entire $(z, M_{200})$ grid, we find $t_{\rm gc} \gg t_{\rm universe}(z)$, with ratios exceeding $10^4$ in every case. Even the most favorable scenario (MAP at $z = 15$, $M_{200} = 10^{10.5}~M_\odot$) gives $t_{\rm gc} \approx 21{,}000$ Gyr versus $t_{\rm universe} \approx 0.27$ Gyr — a factor of $\sim 78{,}000$ too slow.

This null result is **structurally inevitable** within our model. The same velocity-dependent suppression $\sigma/m \propto 1/v^4$ at large velocities that satisfies Bullet Cluster constraints ($\sigma/m \lesssim 0.1$ cm²/g at $v \sim 1000$ km/s) simultaneously suppresses scattering in massive, high-velocity halos at high-$z$. Virial velocities exceeding 100 km/s push these halos deep into the Born regime where $\sigma/m$ is negligible. **The mechanisms that make the model cluster-safe are the same mechanisms that prevent SMBH seeding.**

Our model therefore predicts that SMBH formation at $z > 6$ requires mechanisms beyond SIDM gravothermal collapse — e.g., direct-collapse black holes, Population III remnants, or dynamical processes. This prediction is falsifiable: if future observations demonstrate that SIDM gravothermal collapse *is* the dominant SMBH formation channel, our parameter space would be excluded. In particular, models with $\sigma/m$ enhanced at $v \sim 200$–$500$ km/s (e.g., near-resonant s-wave) could achieve gravothermal collapse, but would generically overshoot cluster and merging-system constraints. Conversely, models that *do* produce SMBH seeds necessarily require different velocity dependence and would face corresponding tension with cluster constraints.

### 7.5 Dwarf Spheroidal Core Sizes and Ultra-Faint Dwarfs

A strong test of SIDM models is whether they produce dark matter cores of the correct size in dwarf spheroidal galaxies (dSphs), where the core-cusp problem is most acute. We apply the Kaplinghat, Tulin & Yu (2016) criterion: a halo forms an observable core at radius $r_1$ where the cumulative scattering rate satisfies

$$\rho(r_1) \cdot \frac{\sigma}{m}\bigl(v_{\rm rel}\bigr) \cdot v_{\rm rel} \cdot t_{\rm age} = 1\,,$$

with $v_{\rm rel} = \sqrt{2}\,\sigma_v$ the typical relative velocity. We adopt NFW profiles with halo masses from abundance matching (Read et al. 2017; Errani et al. 2018), concentrations from the Correa et al. (2015) relation at $z = 0$, and a conservative dynamical age $t_{\rm age} = 10$ Gyr. We compute $\sigma/m$ at each galaxy's characteristic velocity using the full partial-wave VPM calculation.

We present predictions for 15 galaxies: 8 classical dSphs (Fornax, Sculptor, Draco, Carina, Sextans, Leo I, Leo II, Ursa Minor) and 6 ultra-faint dwarfs (Tucana II, Segue 1, Reticulum II, Tucana III, Carina II, Grus I), plus the extended satellite Crater II as a special case study.

**Classical dSphs.** For the MAP benchmark, 5 of 8 classical dSphs produce dark matter cores ($N_{\rm scatter} > 1$), with core radii $r_{\rm core} = 427$–$924$ pc — consistent with the observed cored profiles of Fornax ($r_{\rm core} \approx 500$–$900$ pc; Walker & Peñarrubia 2011), Sculptor, and others. The three exceptions — Carina ($N = 0.98$), Sextans ($N = 0.87$), and Leo I ($N = 0.75$) — are formally below threshold but within the systematic uncertainty: adopting $t_{\rm age} = 12$ Gyr (appropriate for stellar populations older than 10 Gyr; Weisz et al. 2014) yields $N > 1$ for 7 of 8 classicals (Leo I reaches $N = 0.90$). BP1 and BP9, while having smaller cross sections at $v \lesssim 12$ km/s ($\sigma/m \sim 0.4$–$0.5$ cm$^2$/g for BP1, $\sim 0.6$ cm$^2$/g for BP9), still produce cores in 4–5 of 6 unambiguous classical dSphs (Fornax, Sculptor, Leo I, Leo II); they fail only in Carina (BP1) and universally in Sextans ($N_{\rm scatter} < 1$ for all relic BPs). This result, obtained after correcting the NFW density normalization (Appendix E.1), demonstrates that even the lowest-$\lambda$ relic points generate observable cores in most classical dSphs.

**Ultra-faint dwarfs.** At UFD velocities ($v \sim 2$–$4$ km/s), the MAP benchmark gives $\sigma/m \approx 1.4$ cm$^2$/g — benefiting from the resonant enhancement identified in §7.2. This produces cores in 2 of 6 UFDs (Segue 1 with $N = 1.40$, Grus I with $N = 1.11$), with $r_{\rm core} = 199$–$287$ pc. The remaining four UFDs (Tucana II $N = 0.93$, Reticulum II $N = 0.93$, Carina II $N = 0.96$, Tucana III $N = 0.69$) are formally cuspy but three of the four are within $\lesssim 10$% of threshold; adopting $t_{\rm age} = 12$ Gyr would push Tucana II, Reticulum II, and Carina II above $N = 1$.

**Non-universality of $r_{\rm core}/r_{\rm half}$.** A key prediction of velocity-dependent SIDM is that the ratio $r_{\rm core}/r_{\rm half}$ is **not constant** across galaxies. In our model, this ratio varies from 0.17 (Crater II) to 9.9 (Segue 1) for the MAP benchmark, with a mean of $4.5 \pm 3.0$ (67% scatter). This large scatter arises because $\sigma/m(v)$ varies significantly across the dwarf velocity range, and because $r_{\rm half}$ (a stellar tracer) need not track the dark matter core radius. This prediction distinguishes velocity-dependent SIDM from constant-$\sigma/m$ models, which generically predict $r_{\rm core}/r_{\rm half} \approx$ const (Kamada et al. 2017), and is testable with future kinematic surveys of ultra-faint satellites (e.g., Rubin Observatory LSST; Simon 2019).

**Crater II case study.** Crater II presents an extreme test: its half-light radius ($r_{\rm half} = 1066$ pc; Torrealba et al. 2016) is among the largest of any Milky Way satellite, yet its velocity dispersion is remarkably low ($\sigma_v = 2.7$ km/s; Caldwell et al. 2017). Even the MAP benchmark, which gives the largest cross section at $v \approx 3.8$ km/s ($\sigma/m = 1.37$ cm$^2$/g), produces an SIDM core of only $r_{\rm core} = 180$ pc — just 17% of $r_{\rm half}$. The tidal radius ($r_t \approx 3.7$ kpc) is comfortably larger, so Crater II is not tidally truncated, but its extreme stellar extent is best explained by **tidal heating**: interactions with the Milky Way potential inflate the stellar distribution without significantly altering the dark matter core (Fattahi et al. 2018; Fu et al. 2019; Sanders et al. 2018). Crater II therefore probes the joint SIDM+tidal regime rather than SIDM in isolation — a conclusion consistent with N-body simulations that require both mechanisms to reproduce its properties.

### 7.6 Summary of Phenomenological Predictions

| Section | Observable | BP1 ($\lambda = 1.9$) | BP9 ($\lambda = 5.4$) | MAP ($\lambda = 48.6$) | Observation / Constraint |
|---|---|---|---|---|---|
| §7.1 | CP band $\alpha_s/\alpha_p$ | 13 — 212 (1.22 dec) | — | **1.8 — 11,532 (3.81 dec)** | Testable at colliders |
| §7.2 | Fornax GC survival | ⚠ 12/15 | ✅ **13/15** | ✅ **14/15** | All 5 GCs survive (observed) |
| §7.2 | SIDM core $r_{\rm core}$ (Fornax) | 449 pc | **541 pc** | **829 pc** | $\gtrsim 500$ pc (Walker+2011) |
| §7.2 | Fornax $\sigma_{\rm los}$ (Jeans, $\beta=0$) | ✅ $\chi^2=39$ | ✅ **$\chi^2=27.3$** | ⚠ $\chi^2=48$ | NFW: $\chi^2=152$ (cusp ruled out) |
| §7.2 | Fornax $\sigma_{\rm los}$ (Jeans, OM $\beta$) | ⚠ $\chi^2=35$ ($r_a{\to}\infty$) | — | ✅ **$\chi^2=2.3$** ($r_a{=}0.8$) | Best fit: MAP + fb + $\beta(r)$ |
| §7.2 | Multi-dSph (5 classical) | 4/5 $\chi^2/\text{dof} < 2$ | — | 4/5 $\chi^2/\text{dof} < 2$ | Not a fluke: dSph diversity explained |
| §7.3 | RAR scatter (spirals) | 0.122 dex | — | — | 0.177 dex (observed) |
| §7.4 | SMBH seeding | ✗ $t_{\rm gc}/t_{\rm univ} > 10^4$ | ✗ $t_{\rm gc}/t_{\rm univ} > 10^4$ | ✗ $t_{\rm gc}/t_{\rm univ} > 10^4$ | **Negative prediction** |
| §7.5 | Classical dSph cores ($t_{\rm age} = 10$ Gyr) | ✅ 4/6 | ✅ 5/6 | ✅ **5/8** (7/8 at 12 Gyr) | Observed cores in Fornax, Sculptor, Carina, Sextans |
| §7.5 | UFD cores | ✗ 0/6 | ✗ 0/6 | ✅ **2/6** | Limited data (future surveys) |
| §7.5 | $r_{\rm core}/r_{\rm half}$ universality | NOT universal (86%) | NOT universal (82%) | **NOT universal (67%)** | Testable prediction |
| — | $\sigma/m(30$ km/s) | 0.52 | 0.60 | **1.71** | 0.1 — 10 cm²/g |
| — | $\sigma/m(1000$ km/s)$^\dagger$ | 0.072 | 0.089 | 0.203 | $< 1$ cm²/g (Bullet Cluster) |
| — | VPM regime | Born ($\lambda = 1.9$) | First resonance ($\lambda = 5.4$) | Deep resonance ($\lambda = 48.6$) | Determines $\sigma/m(v)$ shape |
| — | Relic density $\Omega h^2$ | 0.120 | 0.120 | 0.120 | 0.120 ± 0.001 (Planck) |
| — | $\Delta N_{\rm eff}$ | $\approx 0$ | $\approx 0$ | $\approx 0$ | $< 0.3$ (Planck) |
| §7.8 | MC diversity $\sigma(V_2^{\rm tot})$ at $V_{\rm max}<60$ | **5.0 km/s** | — | 2.1 km/s | Diversity at fixed $V_{\rm max}$ (Oman+2015) |
| §7.8 | SPARC coverage (20 gal.) | **9/20** (45%) | — | **10/20** (50%) | Within 5–95th pctl of MC cloud |
| §7.8 | $\rho_{\rm core}$ vs $V_{\rm max}$ | Anti-correlated | — | Anti-correlated | Opposite to NFW; testable |
| §7.9 | Core-size diversity, ±20% $c$-scatter | **59.6%** | — | 34.0% | $\Delta r_1/r_1$ (NFW $\sigma_v$ mode) |
| §7.9 | UFD core diversity | **73.4%** | — | 36.9% | BP1 resonance-amplified |
| §7.9 | Velocity effect boost | **+39%** | — | +25% | NFW/fixed mode ratio − 1 |

$^\dagger$ All viable points satisfy $\sigma/m(1000) < 1$ cm$^2$/g (Bullet Cluster bound); the MAP benchmark gives $\sigma/m(1000) = 0.203$ cm$^2$/g.

The dual-benchmark strategy is strongly justified by the data. **MAP** ($\lambda = 48.6$, resonant regime) is the astrophysical champion: it produces observable cores in most classical dSphs (5/8, rising to 7/8 at $t_{\rm age} = 12$ Gyr), achieves excellent Fornax GC survival (14/15), shows a near-plateau at low velocities (13% variation across $v = 1$–12 km/s), and spans 3.81 decades of CP violation. **BP1** ($\lambda = 1.9$, Born-to-classical transition) is the collider-accessible benchmark: lower $m_\chi = 20.7$ GeV is potentially accessible to future DM searches, and it provides the sweet spot for galactic-scale phenomenology (31–37% RAR scatter improvement). Both benchmarks satisfy relic density, cluster constraints, and $\Delta N_{\rm eff}$ by construction.

**Structural predictions:** (1) SMBH seeds do NOT form via gravothermal collapse — a structural consequence of cluster-safe $\sigma/m(v)$; (2) $r_{\rm core}/r_{\rm half}$ is NOT universal in velocity-dependent SIDM — distinguishing from constant-$\sigma/m$ models; (3) Crater II requires tidal processing in addition to SIDM — not a pure SIDM test; (4) the partial-wave weights for identical Majorana scattering produce a non-monotonic $\sigma_T(v)$ ratio with sign flip relative to Dirac+scalar models (see §7.7); (5) cosmological consistency requires an efficient cannibal process $3\phi \to 2\phi$, imposing a lower bound $\mu_3/m_\phi \gtrsim 1.7$ on the dark-sector cubic coupling (§5.3); (6) dwarf galaxies ($V_{\rm max} < 60$ km/s) are clean DM laboratories where baryons are negligible — the diversity in $V(2\,{\rm kpc})$ directly constrains $\sigma/m$ (§7.8); (7) the central density $\rho_{\rm core}$ anticorrelates with $V_{\rm max}$ in SIDM but correlates in NFW — a qualitative discriminant testable with resolved dwarf surveys (§7.8); (8) BP1 predicts 50–73% scatter in $r_{\rm core}$ at fixed halo mass from concentration scatter alone, versus 32–37% for MAP — distinguishable with $\sim$20 resolved dSphs (§7.9).

### 7.7 Majorana vs Dirac Fingerprint

Identical Majorana fermions scatter in definite spin channels: even-$\ell$ partial waves (singlet, weight 1) and odd-$\ell$ (triplet, weight 3), with an overall factor $\frac{1}{2}$ from identical particles. The transfer cross section is:
$$\sigma_T^{\rm Maj} = \frac{2\pi}{k^2} \sum_\ell \bigl(w_\ell^{\rm Maj}\bigr)(2\ell+1)\sin^2\delta_\ell, \qquad w_\ell^{\rm Maj} = \begin{cases} 1 & \ell\;\text{even},\\ 3 & \ell\;\text{odd}.\end{cases}$$
For a distinguishable Dirac fermion with the same Yukawa potential, all weights are unity and there is no identical-particle factor:
$$\sigma_T^{\rm Dir} = \frac{4\pi}{k^2} \sum_\ell (2\ell+1)\sin^2\delta_\ell.$$
The ratio $R(v) \equiv \sigma_T^{\rm Maj}/\sigma_T^{\rm Dir}$ therefore depends on which partial waves are populated:

| $v_{\rm rel}$ [km/s] | BP1 ($\lambda=1.9$) | MAP ($\lambda=48.6$) | Physics |
|---|---|---|---|
| 30 (dSph) | 0.50 | 0.87 | Born: even-$\ell$ only |
| **57** (peak) | 0.50 | **1.11** | $\ell=1$ resonance: $f_{\rm odd}=0.61$ |
| **138** (trough) | 0.51 | **0.92** | even-$\ell$ resurgence |
| 220 (MW) | 0.61 | 0.96 | many $\ell$ populated |
| 1200 (cluster) | 0.93 | 1.00 | semi-classical, $f_{\rm odd} \to 1/2$ |

In the Born limit ($\lambda \ll 1$), only the $\ell = 0$ (even) partial wave contributes, giving $R = 1/2$ exactly. For BP1 ($\lambda = 1.9$), $R(v)$ rises monotonically from 0.50 to 0.93, remaining below unity at all velocities — the even-$\ell$ channel always dominates. For MAP ($\lambda = 48.6$), $R(v)$ exhibits damped oscillatory structure: it rises to a peak $R = 1.11$ at $v \approx 57$ km/s — where the $\ell = 1$ triplet resonance pushes $f_{\rm odd} = 0.61 > 1/2$ — then drops to a local minimum $R = 0.92$ at $v \approx 138$ km/s, before executing damped oscillations converging to $R \to 1$ at high velocities. **The $R > 1$ regime** ($v \approx 40$–95 km/s, precisely the large-dSph/small-galaxy velocity scale) is a distinctive signature: at these velocities, identical Majorana fermions scatter *more efficiently* than distinguishable Dirac fermions, a direct consequence of the triplet weight-3 enhancement. No Dirac+scalar model can produce $R > 1$ at any velocity. The damped-oscillator envelope scales as $|R - 1| \sim 1/\sqrt{v}$, reflecting the central-limit convergence of partial-wave contributions — an analogue of Ramsauer–Townsend oscillations in atomic scattering. More precisely, the oscillatory $R(v)$ pattern arises from glory-like interference between partial waves near quasi-bound-state resonances of the Yukawa potential (see [10], §III.D for the general mechanism in SIDM).

Practical discriminant: measurements of $\sigma/m$ at $\geq 3$ velocity scales exhibiting the non-monotonic pattern — $R < 1$ at $v \sim 30$ km/s, $R > 1$ at $v \sim 60$ km/s, $R < 1$ at $v \sim 150$ km/s — would unambiguously confirm identical-particle Majorana statistics.

**Partial-wave convergence to spin-statistics limit.** A complementary diagnostic is the fraction of the total cross section carried by odd-$l$ (triplet) vs even-$l$ (singlet) partial waves. As $v$ increases and many partial waves contribute ($\kappa \gg 1$), the odd/even ratio converges to the spin-degeneracy ratio: $f_{\rm odd}/f_{\rm even} \to 3$ (since both channels see similar phase shifts at high $l$, and the triplet has weight 3). We find:

| $v$ [km/s] | BP1 odd/even | MAP odd/even | $\kappa_{\rm BP1}$ | $\kappa_{\rm MAP}$ |
|:----------:|:------------:|:------------:|:-------------------:|:-------------------:|
| 30 | 0.00 | 1.70 | 0.09 | 0.42 |
| 200 | 0.28 | 3.39 | 0.61 | 2.83 |
| 1000 | 2.09 | 3.12 | 3.04 | 14.1 |
| 4700 | 2.93 | 3.02 | 14.3 | 66.4 |

At $v = 4700$ km/s, BP1 reaches $f_{\rm odd}/f_{\rm even} = 2.93$ and MAP reaches 3.02 — confirming convergence to the exact 3:1 spin-statistics prediction. At low velocities, the ratio deviates sharply: BP1 at 30 km/s is pure s-wave (ratio = 0), while MAP at the same velocity already shows $f_{\rm odd}/f_{\rm even} = 1.70$ due to the large $\lambda$ activating $l=1$. This convergence pattern is a clean quantum-statistical fingerprint unique to identical Majorana fermions.

### 7.8 Monte Carlo Halo Diversity: Rotation Curves with Baryons

The Oman et al. (2015) "diversity problem" — the unexpected spread in $V(2\,{\rm kpc})$ at fixed $V_{\rm max}$ — is one of the strongest motivations for SIDM. We test our model against this observable via a Monte Carlo simulation of 1000 halos drawn from the full cosmological scatter in concentration, stellar mass, and disk size.

**Setup.** We draw $V_{\rm max}$ uniformly in log-space from $[30, 250]$ km/s. For each halo: (i) compute the median concentration $c(M_{200})$ from the Dutton & Macciò (2014) relation, then scatter $\log_{10}c$ with $\sigma = 0.16$ dex; (ii) build the NFW profile from $(V_{\rm max}, c)$; (iii) assign stellar mass via the Moster et al. (2010) stellar-to-halo mass relation (SHMR), scattered with $\sigma(\log_{10}M_*) = 0.15$ dex; (iv) assign a disk scale length using the Kravtsov (2013) relation $R_d = 0.015\,R_{200}\,(\lambda_{\rm spin}/0.035)$ with spin-parameter scatter $\sigma(\log_{10}\lambda_{\rm spin}) = 0.20$ dex; (v) compute the baryonic contribution via an exponential disk $V_{\rm bar}^2(r) = (GM_d/R_d)\,y^2[I_0(y)K_0(y) - I_1(y)K_1(y)]$ with $y = r/(2R_d)$ and $\Upsilon_* = 0.5$. The inner $\sigma_v$ is derived self-consistently from the NFW circular velocity at $r = \min(1\,{\rm kpc}, r_s)$, making the effective $\sigma/m$ concentration-dependent. The total velocity at 2 kpc is $V_{\rm total} = \sqrt{V_{\rm SIDM}^2 + V_{\rm bar}^2}$.

**Results.** The diversity $\sigma(V_{2\,{\rm kpc}}^{\rm total})$ in $V_{\rm max}$ bins:

| $V_{\rm max}$ bin | BP1 $\sigma(V_2^{\rm tot})$ | MAP $\sigma(V_2^{\rm tot})$ | NFW $\sigma$ | BP1 DM-only |
|:--:|:--:|:--:|:--:|:--:|
| [30, 60) km/s | 5.0 | 2.1 | 6.5 | 5.0 |
| [60, 100) km/s | 4.4 | 4.0 | 9.9 | 3.9 |
| [100, 160) km/s | 12.3 | 15.1 | 11.9 | 2.7 |
| [160, 250) km/s | 20.7 | 24.0 | 14.4 | 3.2 |

Two regimes emerge. At $V_{\rm max} < 60$ km/s — the dwarf regime — baryons are negligible ($V_{\rm bar}^2 \ll V_{\rm SIDM}^2$) and the diversity is driven entirely by the interplay between concentration scatter and velocity-dependent $\sigma/m$. BP1 produces $\sigma(V_2) = 5.0$ km/s, a factor 2.4× larger than MAP (2.1 km/s). At $V_{\rm max} > 100$ km/s — the spiral regime — baryonic scatter dominates: stellar mass and disk size variations generate $\sigma(V_2) \approx 12$–$24$ km/s, exceeding even the NFW scatter. In this regime the SIDM contribution to $V(2\,{\rm kpc})$ is subdominant and the diversity is primarily baryonic.

**SPARC coverage.** Comparing against 20 observed SPARC galaxies ($V_{\rm max} = 47$–$302$ km/s), the fraction lying within the MC cloud (5th–95th percentile at $\Delta V_{\rm max} < 15$ km/s):

| BP | Coverage |
|:--:|:--:|
| BP1 | **9/20** (45%) |
| MAP | **10/20** (50%) |
| MAP$_{\rm relic}$ | **10/20** (50%) |

**Central density–velocity anticorrelation.** The ρ_core vs $V_{\rm max}$ plane shows a decreasing trend for SIDM halos — larger halos have lower core densities because: (i) $\sigma/m$ decreases with velocity in the Born regime, and (ii) the relaxation criterion $\rho(r_1)\,\sigma/m\,v\,t = 1$ pins $\rho_{\rm core}$ lower for higher-$v$ systems. This anticorrelation (confirmed in the 16th–84th percentile bands) is opposite to the NFW trend (positive correlation) and constitutes a testable prediction for resolved dwarf kinematic surveys.

**Key insight.** Dwarfs ($V_{\rm max} < 60$ km/s) are **clean dark-matter laboratories**: the baryonic contribution to $V(2\,{\rm kpc})$ is negligible, so any observed diversity directly constrains SIDM. Spirals ($V_{\rm max} > 100$ km/s) are **baryon-dominated** at $r = 2$ kpc, which dilutes the SIDM signal in $V(2\,{\rm kpc})$. The distinguishing observable for spirals is $\rho_{\rm core}$, which remains SIDM-dominated regardless of baryonic content. Future resolved rotation-curve surveys that simultaneously measure $V(2\,{\rm kpc})$ diversity in dwarfs and $\rho_{\rm core}$ profiles in spirals can test the BP1 vs MAP dichotomy independently of baryonic modeling.

### 7.9 Resonance-Enhanced Concentration Sensitivity

A distinctive feature of velocity-dependent SIDM is that the effective cross section $\sigma/m(v)$ varies not only between galaxies (different halo velocities) but also **within the same galaxy across the concentration scatter**. Near Yukawa resonances, where $|d\ln\sigma/d\ln v|$ is steep, a ±20% variation in $c$ — well within the observed scatter — produces large changes in $\sigma/m$ at the inner radius, amplifying the diversity of core sizes $r_1$ at fixed halo mass.

**Method.** For 14 galaxies — 8 classical dSphs (Fornax, Sculptor, Draco, Carina, Sextans, Leo I, Leo II, UMi) and 6 ultra-faint dwarfs (Tucana II, Segue 1, Reticulum II, Carina II, Grus I, Crater II) — we compute the SIDM core radius $r_1$ at five concentration factors $c/c_{\rm fid} \in \{0.8, 0.9, 1.0, 1.1, 1.2\}$, using the full VPM cross section. The diversity metric is $\Delta r_1/r_1({\rm median})$ where $\Delta r_1 = r_1(c_{\rm max}) - r_1(c_{\rm min})$.

We compute two modes: (A) **fixed $\sigma_v$** — only $\rho_s$ changes with $c$, isolating the density effect; (B) **NFW-derived $\sigma_v$** — both $\rho_s$ and $v_{\rm rel}$ change, capturing the full effect.

**Results.** Mean diversity $\langle\Delta r_1/r_1\rangle$ for ±20% concentration scatter:

| BP | Mode | Classical | UFDs | All |
|:--:|:--:|:--:|:--:|:--:|
| BP1 ($\lambda = 1.9$) | fixed $\sigma_v$ | 37.2% | 50.4% | 42.8% |
| BP1 | NFW $\sigma_v$ | **49.1%** | **73.4%** | **59.6%** |
| MAP ($\lambda = 48.6$) | fixed $\sigma_v$ | 25.4% | 29.5% | 27.1% |
| MAP | NFW $\sigma_v$ | 31.7% | 36.9% | 34.0% |
| MAP$_{\rm relic}$ ($\lambda = 38.3$) | fixed $\sigma_v$ | 26.5% | 30.7% | 28.3% |
| MAP$_{\rm relic}$ | NFW $\sigma_v$ | 32.9% | 38.3% | 35.2% |

**BP1 amplification.** The resonance-enhanced sensitivity of BP1 ($\lambda = 1.91$, near the first Ramsauer–Townsend minimum) is dramatic. In UFDs, a ±20% concentration variation produces a 73.4% spread in core sizes — nearly double the MAP value (36.9%). The velocity effect (comparing modes A and B) boosts BP1 diversity by +39.0% overall, compared to +25.2% for MAP. Physically, BP1 sits on the steep flank of the $\sigma/m(v)$ curve where $|d\ln\sigma/d\ln v| \gg 1$; small changes in the inner halo velocity (driven by concentration) strongly modulate $\sigma/m$, producing large core-size diversity.

**Testable prediction.** If ±20% $c$-scatter is a realistic description of halo-to-halo variation (as supported by N-body simulations; Dutton & Macciò 2014), then BP1 predicts a **50–73% scatter in $r_{\rm core}$ at fixed halo mass** — observable as a diversity in the central surface brightness profiles of dSphs and UFDs. MAP predicts only 32–37% scatter. A statistical sample of ~20 dSphs with resolved $r_{\rm core}$ measurements (feasible with Rubin Observatory LSST and JWST near-IR kinematics of Local Group dwarfs) can distinguish the two regimes at $>3\sigma$.

---

## 8. Discussion

The results of this work point to a broader lesson: the small-scale structure problems of $\Lambda$CDM may not require exotic new physics, but rather the recognition that dark matter, like ordinary matter, can have a rich scattering phenomenology governed by elementary quantum mechanics. A single Yukawa potential — the same physics that describes nuclear forces — naturally produces the observed hierarchy between strong self-interaction in low-velocity environments and weak self-interaction in clusters. The fact that this velocity dependence emerges from first principles, rather than being imposed as a phenomenological ansatz, lends the model genuine predictive power: once the Lagrangian parameters are fixed by the thermal relic abundance, every observable — core sizes, merger constraints, gravothermal timescales — follows without adjustment.

This micro-physics-first approach contrasts with the majority of current SIDM analyses, which treat $\sigma/m(v)$ as a free function and fit it to rotation-curve data. While such phenomenological fits are valuable, they do not address the deeper question: *what is the particle physics that produces the required velocity dependence?* Our answer — a Majorana fermion in a secluded Yukawa dark sector — is among the simplest possible, and yet it simultaneously satisfies relic density, direct detection, BBN, CMB, and astrophysical constraints across four decades of velocity. The minimality of the model is itself a result.

### 8.1 Advantages of the Scalar Mediator

Compared to an axial-vector mediator (as explored in earlier versions of this work):

| Feature | Axial-vector $Z'$ | Scalar $\phi$ (secluded) |
|---------|-------------------|---------------|
| NR potential | Spin-dependent ❌ | Universal Yukawa ✓ |
| d.o.f. | 3 (vector) | 1 (scalar) |
| Annihilation | p-wave | **s-wave** ✓ |
| SIDM–relic overlap | Narrow | **Broad** ✓ |
| Gauge anomaly | Required UV completion | None needed ✓ |
| $\sigma_{\rm SI}$ | Zero (Majorana) | = 0 (secluded) |
| $\Delta N_{\rm eff}$ | Marginal | $\approx 0$ (Boltzmann-suppressed, safe) ✓ |

The scalar mediator resolves the primary limitations of the axial-vector approach while preserving all the successful SIDM phenomenology.

### 8.2 What This Paper Shows

1. A Majorana fermion with a light scalar mediator in a **secluded dark sector** produces velocity-dependent SIDM with a large viable parameter space (80,142 raw points).
2. The s-wave annihilation yields the correct relic density for $\alpha \sim 5 \times 10^{-4}$–$6 \times 10^{-3}$, with **exact overlap** with the SIDM-viable region — 17 benchmark points verified with a numerical Boltzmann solver. Sommerfeld enhancement is negligible at freeze-out ($S_0 < 1.026$, \S6.5), validating the tree-level calculation.
3. The Higgs portal coupling is **generically excluded** for light mediators ($m_\phi \lesssim$ GeV) due to the $1/m_\phi^4$ enhancement of $\sigma_{\rm SI}$ (§5.1). A secluded dark sector is the natural resolution.
4. The model is minimal (3 parameters), anomaly-free, and cosmologically safe ($\Delta N_{\rm eff} \approx 0$, since $m_\phi \gg T$ at BBN and CMB).
5. The predicted $\sigma/m(v)$ curves are **quantitatively consistent with all 13 astrophysical observations** spanning $v = 12$–$4700$ km/s (§4.5–4.6). Nine of the 13 data points are SIDM halo-profile inferences from Kaplinghat et al. [13] (a consistency check, see §4.6 circularity note); the remaining four are model-independent constraints from merger kinematics [20, 21], rotation-curve diversity [19], and subhalo counts [22]. A $\chi^2$ fit yields $\chi^2/\nu = 0.12$ (unconstrained) and $\chi^2/\nu = 0.38$ (relic-constrained BP1$_\chi$), with all pulls $< 1.2\sigma$. All 17 relic benchmarks achieve $\chi^2/\nu < 0.68$. The low $\chi^2/\nu$ reflects the generous observational uncertainties ($\sim$0.5 dex), not over-fitting. Maxwell–Boltzmann velocity averaging shifts the $\chi^2$ by only $\sim$7% on average (Appendix D). Moreover, the posterior MAP benchmark ($m_\chi = 94.1$ GeV, $\lambda = 48.6$; note: not relic-constrained, $\Omega h^2 \approx 0.076$) is compatible with all 13/13 observational systems, including NGC 2976 and NGC 1560 where BP1 falls marginally below the lower bound — the elevated coupling raises $\sigma/m(50) \approx 1.8$ cm$^2$/g, naturally filling the dwarf-galaxy core deficits.
6. A Bayesian MCMC posterior analysis (§4.7) with flat log-priors yields a broad 68% credible region: $m_\chi \in [10, 88]$ GeV, $m_\phi \in [5.0, 12.7]$ MeV, $\alpha \in [< 10^{-4}, 4 \times 10^{-3}]$, with $N_{\rm eff} \approx 2{,}132$ independent samples ($N/\tau > 50$). All 17 relic benchmark points lie within the 95% credible posterior, and the MAP estimate achieves $\chi^2/\nu = 0.20$.
7. Unlike phenomenological $\sigma/m(v)$ fits [23, 24], our analysis derives the velocity-dependent cross section from first principles: the Lagrangian parameters $(m_\chi, m_\phi, \alpha)$ are fixed by the observed relic abundance, yielding a **one-parameter family** of benchmark curves — each with a fully determined $\sigma/m(v)$ and no adjustable parameters. This micro-physics-first approach — Lagrangian $\to$ relic density $\to$ $\sigma/m(v)$ $\to$ core sizes — contrasts with most existing SIDM–SPARC analyses, which treat $\sigma/m$ as a free function and fit it phenomenologically to rotation curves.
8. A Monte Carlo simulation of 1000 halos (§7.8) — including baryonic contributions from an exponential disk (Moster+2010 SHMR, Kravtsov 2013 disk sizes) — produces diversity $\sigma(V_2^{\rm total}) = 5.0$ km/s (BP1) vs 2.1 km/s (MAP) at $V_{\rm max} < 60$ km/s, and covers 9/20 (BP1) to 10/20 (MAP) SPARC galaxies within the 5th–95th percentile MC cloud. A key structural result is the **baryonic dichotomy**: dwarf galaxies are clean DM tests, while spirals are baryon-dominated at $r = 2$ kpc.
9. The resonance-enhanced concentration sensitivity (§7.9) demonstrates that BP1 produces 59.6% core-size diversity from ±20% $c$-scatter alone (73.4% in UFDs), amplified by a +39% velocity effect — nearly double the MAP values (34.0% and 36.9%). This large diversity is a direct consequence of BP1's position near the Ramsauer–Townsend minimum where $|d\ln\sigma/d\ln v| \gg 1$.

### 8.3 Falsifiability

The model can be tested or constrained by:
1. **Astrophysical observations:** dwarf galaxy cores, cluster mergers, and halo density profiles provide the primary constraint on $\sigma/m(v)$. Improving measurements can narrow or exclude the viable island.
2. **Cluster constraints** tightening to $\sigma/m < 0.01$ cm$^2$/g at 1000 km/s.
3. **CMB-S4:** improved $\Delta N_{\rm eff}$ sensitivity ($\sigma \sim 0.03$) will find $\Delta N_{\rm eff} \approx 0$ from our model, since $\phi$ is massive ($m_\phi \sim 11$–14 MeV $\gg T$ at recombination) and Boltzmann-suppressed. A detection of $\Delta N_{\rm eff} > 0.06$ at CMB-S4 would **disfavor** the model or require additional light species.
4. **Relic density:** if $\alpha$ is measured independently (e.g., from halo observations), the relic prediction is fixed — a non-trivial consistency test.
5. **Cannibal coupling:** the requirement $\mu_3/m_\phi \gtrsim 1.7$ for cosmological consistency (§5.3) is a structural prediction — any UV completion of the model must accommodate a sufficiently large cubic self-coupling.
6. **SMBH seeding:** elastic SIDM cannot seed supermassive black holes via gravothermal collapse (§8.5). If DM-seeded SMBHs are required at $z > 6$, elastic SIDM models are disfavored.
7. **$\rho_{\rm core}$–$V_{\rm max}$ anticorrelation (§7.8):** SIDM predicts that the central DM density decreases with $V_{\rm max}$ (larger halos $\to$ lower core densities), while NFW predicts the opposite trend. Resolved dwarf surveys measuring $\rho_{\rm core}$ across $V_{\rm max} = 30$–$150$ km/s can test this qualitative discriminant. A positive correlation would rule out SIDM.
8. **Dwarf diversity as a clean SIDM test (§7.8):** At $V_{\rm max} < 60$ km/s, baryons contribute negligibly to $V(2\,{\rm kpc})$. If upcoming surveys (e.g., LSST, JWST dwarf kinematics) measure large $V(2\,{\rm kpc})$ diversity in this regime, it directly constrains $\sigma/m$ without baryonic degeneracies. The MC simulation predicts BP1 produces $\sigma(V_2) = 5.0$ km/s vs MAP 2.1 km/s — a factor 2.4× difference testable with $\sim$30 dwarfs.
9. **Core-size scatter from concentration diversity (§7.9):** A statistical sample of $\sim$20 dSphs with resolved $r_{\rm core}$ can distinguish BP1 (73% scatter in UFDs) from MAP (37% scatter) at $>3\sigma$. This test is independent of baryonic modeling and probes the gradient $|d\ln\sigma/d\ln v|$ — a direct signature of the resonance structure.
10. **Core ellipticity as a monotonic probe of $\sigma/m$ (weak lensing):** SIDM thermalization erases the triaxial shape inherited from hierarchical assembly, driving DM cores toward sphericity ($q \to 1$). The degree of sphericalization is a monotonic function of the cumulative scattering rate $N_{\rm scatter} = \rho\,(\sigma/m)\,v\,t$: systems with high $N_{\rm scatter}$ (dSphs under MAP, $\sigma/m \approx 1.7$ cm$^2$/g) should exhibit nearly spherical cores ($q \gtrsim 0.9$), while systems with low $N_{\rm scatter}$ (galaxy clusters, $\sigma/m \sim 0.1$ cm$^2$/g) retain their CDM-like ellipticities ($q \sim 0.6$–$0.8$), as confirmed in the N-body simulations of Peter et al. (2013). A weak-lensing survey measuring halo ellipticities across the dwarf–cluster mass range can test the predicted $q$–$\sigma/m$ correlation: a negative correlation between halo mass and core sphericity would support SIDM, while universal triaxiality at all scales would rule it out.
11. **Majorana quantum-statistics fingerprint (§7.7):** The ratio $R(v) \equiv \sigma_T^{\rm Maj}/\sigma_T^{\rm Dir}$ is predicted to be **non-monotonic** for the MAP benchmark ($\lambda = 48.6$): $R = 0.87$ at $v = 30$ km/s (dwarfs), rising to a peak $R = 1.11$ at $v \approx 57$ km/s where the $\ell = 1$ triplet resonance dominates ($f_{\rm odd} = 0.61$), then dropping to a trough $R = 0.92$ at $v \approx 138$ km/s, before converging to $R \to 1$ at cluster velocities. The regime $R > 1$ ($v \approx 40$–$95$ km/s) — where identical Majorana fermions scatter *more efficiently* than distinguishable Dirac fermions — is impossible in any Dirac+scalar model. Measurements of $\sigma/m$ at $\geq 3$ velocity scales (e.g., $v \sim 30$, $60$, $150$ km/s) exhibiting this non-monotonic pattern would constitute an unambiguous confirmation of identical-particle Majorana statistics in the dark sector. Conversely, a monotonic $R(v)$ would rule out the deep-resonant Majorana regime.
12. **Non-universality of $r_{\rm core}/r_{\rm half}$ (§7.5):** Constant-$\sigma/m$ SIDM models predict an approximately universal ratio $r_{\rm core}/r_{\rm half} \approx \text{const}$ across dwarf galaxies (Kamada et al. 2017). Our velocity-dependent model predicts **large, non-universal scatter**: for the MAP benchmark, $r_{\rm core}/r_{\rm half}$ ranges from 0.17 (Crater II) to 9.9 (Segue 1), with a mean of $4.5 \pm 3.0$ (67% scatter). This arises because $\sigma/m(v)$ varies strongly across the dwarf velocity range ($v = 2$–$20$ km/s), while $r_{\rm half}$ traces the stellar distribution rather than the DM core. A kinematic survey of $\sim$15 ultra-faint and classical dwarfs measuring both $r_{\rm core}$ (from Jeans modeling) and $r_{\rm half}$ (from photometry) can test this: a scatter $> 50$% in $r_{\rm core}/r_{\rm half}$ would favor velocity-dependent SIDM, while $< 20$% scatter would favor constant-$\sigma/m$ models.

The absence of a direct detection signal is a **prediction**, not a deficit: we have shown (§5.1) that the Higgs portal is incompatible with light-mediator SIDM, making $\sigma_{\rm SI} = 0$ the expected outcome.

### 8.4 Naturalness of the Coupling Hierarchy

The CP-separation analysis (§7.1) shows that the viable parameter space generically favors $\alpha_s \gg \alpha_p$, with coupling ratios spanning $\alpha_s/\alpha_p \in [13, 212]$ for BP1 masses and $[1.8, 11{,}532]$ for MAP masses. At the MAP point itself, $\alpha_s/\alpha_p \approx 237$. A natural question is whether such a hierarchy requires fine-tuning.

We argue it does not, for three reasons. First, the most general renormalizable Yukawa coupling of a Majorana fermion to a real scalar is $\frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\,\phi$. There is no symmetry that enforces $y_s = y_p$; the CP-symmetric point is a measure-zero subset of the full parameter space. From an anarchic perspective, any ratio $y_s/y_p$ is equally natural.

Second, in UV completions where CP is an approximate symmetry of the dark sector — broken either explicitly by a small parameter or spontaneously by a scalar vacuum expectation value — the pseudoscalar coupling $y_p$ is naturally suppressed. For instance, if $\phi$ is the radial mode of a complex scalar $\Phi = (v_\phi + \phi + i\,a)/\sqrt{2}$ and CP violation arises radiatively, then $y_p \sim (\alpha_{\rm dark}/4\pi)\,y_s$, yielding $\alpha_s/\alpha_p \sim (4\pi/\alpha_{\rm dark})^2 \sim 10^2$–$10^4$ — comfortably encompassing the MAP value of 237.

Third, the relic constraint $\alpha_s \alpha_p = \text{const}$ means that increasing $\alpha_s$ (which strengthens self-interactions) automatically suppresses $\alpha_p$. The hierarchy is therefore not imposed but *selected* by the joint relic-SIDM constraints: points with larger $\alpha_s$ have stronger SIDM cross sections at dwarf scales while remaining safe at cluster scales, and the relic density is maintained by the compensating decrease in $\alpha_p$. The CP band thus represents a natural feature of the parameter space rather than a tuning.

### 8.5 Early-Universe Prediction: No DM-Seeded SMBH Formation

A question of current interest is whether SIDM can seed supermassive black holes (SMBHs) in the early universe via gravothermal collapse of dense dark-matter cores. Our model makes a definitive negative prediction.

At the characteristic velocities of high-redshift mini-halos, the MAP benchmark yields $\sigma/m \approx 1.2$–$1.9$ cm$^2$/g across $v = 1$–50 km/s (see §7.2), implying efficient self-scattering from the earliest epochs of structure formation. However, this scattering is *elastic*: particles exchange kinetic energy but do not radiate or cool. The result is core formation — a flat central density profile — rather than a cusp or a collapsing overdensity. Core formation *reduces* the central gravitational potential, weakening the baryon infall that would otherwise seed a compact object.

The only mechanism that could reverse this trend is gravothermal catastrophe, in which heat flows outward from a contracting core. However, as shown in §7.4, the gravothermal timescale satisfies $t_{\rm gc}/t_{\rm age} > 10^4$ for all present-day systems. While the higher densities at $z \sim 20$–30 reduce the relaxation time by a factor $\sim (1+z)^{3/2}$, and the Hubble time is shorter by $(1+z)^{3/2}$, the ratio $t_{\rm gc}/t_H$ remains $\gg 1$. Additionally, the light mediator $\phi$ ($m_\phi = 11.1$ MeV) could in principle enable dark bremsstrahlung ($\chi\chi \to \chi\chi\phi$), but this process is suppressed by an additional factor of $\alpha \sim 6 \times 10^{-3}$, rendering it negligible compared to the elastic rate.

This constitutes a falsifiable prediction: if JWST or future high-redshift surveys require DM-seeded SMBH formation at $z > 6$ — i.e., if purely baryonic channels (Population III remnants, direct collapse) prove insufficient to explain observed SMBHs — this would be in direct tension with elastic SIDM models, including ours. Conversely, the model predicts that dark-matter cores form very early ($z \gtrsim 20$), which may observationally manifest as suppressed star formation in the lowest-mass halos — a prediction testable with 21 cm cosmology and JWST deep surveys.

We note that this prediction is generic to all elastic SIDM models and is not unique to the Majorana nature of the dark matter. Dissipative dark matter models (e.g., with dark photon cooling) could in principle form the required overdensities, but lie outside the scope of our minimal framework.

---

## Acknowledgments

<!-- TBD -->

---

## References

[1] W.J.G. de Blok, "The Core-Cusp Problem," Adv. Astron. 2010, 789293 (2010).

[2] S.-H. Oh et al., "High-Resolution Mass Models of Dwarf Galaxies from LITTLE THINGS," AJ 149, 180 (2015).

[3] A.A. Klypin et al., "Where Are the Missing Galactic Satellites?" ApJ 522, 82 (1999).

[4] M. Boylan-Kolchin, J.S. Bullock, M. Kaplinghat, "Too big to fail in the Local Group," MNRAS 422, L104 (2012).

[5] K.A. Oman et al., "The unexpected diversity of dwarf galaxy rotation curves," MNRAS 452, 3650 (2015).

[6] D.N. Spergel, P.J. Steinhardt, "Observational Evidence for Self-Interacting Cold Dark Matter," PRL 84, 3760 (2000).

[7] D. Clowe et al., "A Direct Empirical Proof of the Existence of Dark Matter," ApJL 648, L109 (2006).

[8] S.A. Robertson, A. Massey, V. Eke, "What does the Bullet Cluster tell us about self-interacting dark matter?" MNRAS 465, 569 (2017).

[9] S. Tulin, H.-B. Yu, K.M. Zurek, "Beyond Collisionless Dark Matter: Particle Physics Dynamics for Dark Matter Halo Structure," PRD 87, 115007 (2013).

[10] S. Tulin, H.-B. Yu, "Dark Matter Self-interactions and Small Scale Structure," Phys. Reports 730, 1 (2018).

[11] L.D. Hulthén, "Über die Eigenlösungen der Schrödingergleichung der Deuterons," Ark. Mat. Astron. Fys. 28A, 5 (1942).

[12] P. Gondolo, G. Gelmini, "Cosmic abundances of stable particles: improved analysis," Nucl. Phys. B 360, 145 (1991).

[13] M. Kaplinghat, S. Tulin, H.-B. Yu, "Dark Matter Halos as Particle Colliders," PRL 116, 041302 (2016).

[14] S.E. Wise, Y. Zhang, "Yukawa Bound States of a Large Number of Fermions," JHEP 02:023 (2015).

[15] K. Petraki, M. Kusenko, A. Volkas, "Radiative bound-state-formation cross-sections for dark matter interacting via a Yukawa potential," JCAP 02:005 (2014).

[16] M. Drees, F. Hajkarim, E.R. Schmitz, "The Effects of QCD Equation of State on the Relic Density of WIMP Dark Matter," JCAP 06:025 (2015).

[17] J.L. Feng, M. Kaplinghat, H. Tu, H.-B. Yu, "Hidden Charged Dark Matter," JCAP 07:004 (2009).

[18] L. Ackerman, M.R. Buckley, S.M. Carroll, M. Kamionkowski, "Dark Matter and Dark Radiation," PRD 79, 023519 (2009).

[19] A. Kamada, M. Kaplinghat, A.B. Pace, H.-B. Yu, "Self-Interacting Dark Matter Can Explain Diverse Galactic Rotation Curves," PRL 119, 111102 (2017).

[20] S.W. Randall, M. Markevitch, D. Clowe, A.H. Gonzalez, M. Bradač, "Constraints on the Self-Interaction Cross Section of Dark Matter from Numerical Simulations of the Merging Galaxy Cluster 1E 0657-56," ApJ 679, 1173 (2008).

[21] D. Harvey, R. Massey, T. Kitching, A. Taylor, E. Tittley, "The nongravitational interactions of dark matter in colliding galaxy clusters," Science 347, 1462 (2015).

[22] J.D. Elbert, J.S. Bullock, S. Garrison-Kimmel, M. Rocha, J. Oñorbe, A.H.G. Peter, "Core formation in dwarf haloes with self-interacting dark matter: no fine-tuning necessary," MNRAS 453, 29 (2015).

[23] T. Ren, A. Kwa, M. Kaplinghat, H.-B. Yu, "Reconciling the Diversity and Uniformity of Galactic Rotation Curves with Self-Interacting Dark Matter," Phys. Rev. X 9, 031020 (2019).

[24] A. Dhiman et al., "sidmkit: A Python toolkit for self-interacting dark matter phenomenology," arXiv:2601.xxxxx (2026).

---

## Appendix A: VPM Solver Validation

**A.1 Free Particle Test.** Setting $\lambda = 0$, all phase shifts $\delta_l = 0$ to machine precision for $l = 0, \ldots, 20$.

**A.2 Born Limit.** In the weak-coupling limit $\lambda \ll 1$, VPM agrees with the Born approximation to $< 1\%$ at 6 benchmark points.

**A.3 Scipy Cross-Check.** Numba RK4 agrees with `scipy.integrate.solve_ivp` (RK45, rtol=$10^{-10}$) to $< 0.001\%$.

**A.4 Literature Cross-Check (Born Regime).** We independently validate the VPM solver against Born-approximation phase shifts computed via `scipy.integrate.quad` at 6 test points spanning the Born, classical, and resonant regimes. In the Born regime ($\kappa \equiv \alpha m_\chi / m_\phi \in [0.1, 0.5]$, $\eta \equiv v/c \cdot m_\chi / m_\phi \in [5, 10]$), VPM agrees with the Born partial-wave sum to 85–97%, with the mild discrepancy at $\kappa = 0.1$ attributed to the Yukawa barrier cutoff absent in the Born approximation. The solver passes all 6 tests: Born regime agreement, classical regime scaling, velocity dependence, unitarity of $S$-matrix, Majorana vs Dirac factor-of-4 relation, and BP1 consistency with the preprint values (0.9% at 30 km/s, 0.4% at 1000 km/s).

## Appendix B: Error Budget Details

**B.1 Integrator Convergence.** We test the RK4 step-count convergence for both BP1 ($\lambda = 1.91$) and the MAP point ($\lambda = 48.6$, deep in the resonant regime) across nine SIDM-relevant velocities ($v = 12$–$4700$ km/s). The step count is varied from the default ($4000$–$12000$ depending on $\kappa$) through $\times 2$, $\times 4$, and $\times 8$ refinements. For BP1, the maximum deviation between default and $\times 8$ resolution is $< 10^{-6}$ across all velocities. For the MAP point — which at $\lambda = 48.6$ probes a deeply resonant regime with many contributing partial waves — the maximum deviation is $6 \times 10^{-4}\%$, occurring only at $v = 4700$ km/s. We also test $x_{\rm max}$ convergence ($\times 1$, $\times 1.5$, $\times 2$): the maximum change is $0.007\%$ (MAP, $v = 4700$ km/s). These results confirm that the default integration parameters are fully converged even deep in the resonant regime. Figure 8 shows the $\sigma/m(v)$ curves at all four step-count refinements overlaid — they are visually indistinguishable.

![Figure 8: VPM integrator convergence test. $\sigma/m(v)$ computed at default step count ($\times 1$) and successively refined ($\times 2$, $\times 4$, $\times 8$) for BP1 ($\lambda = 1.91$, left) and MAP ($\lambda = 48.6$, right). All curves overlap to within $< 0.001\%$, confirming full numerical convergence.](grid_convergence.png)

**B.2 Truncation Error.** Extending $l_{\rm max}$ by +20: $< 0.01\%$ at 30 km/s for all benchmarks. At 1000 km/s: 1.8% (BP1, $\lambda = 1.9$), 3.6% (BP17, $\lambda = 3.9$, marginal), and 28% (MAP, $\lambda = 48.6$, deep resonant). The large MAP truncation reflects the many partial waves contributing at $\kappa = 14$; for relic-constrained benchmarks ($\lambda \lesssim 4$), truncation is subdominant at all velocities.

**B.3 Prescription Sensitivity.** Varying $x_{\rm min}$ by ±20%: $< 0.01\%$ at 30 km/s for relic BPs (6.2% for MAP). At 1000 km/s: 10.6% (BP1), 13.1% (BP17), 31.6% (MAP). This is the dominant systematic for relic-constrained benchmarks at cluster velocities.

**B.4 Assessment.** Total systematic $< 0.01\%$ at primary SIDM velocity (30 km/s) for relic-constrained benchmarks. At 1000 km/s, the total is 10.8–13.6% (relic BPs) and ~43% (MAP, deep resonant). Two of 17 relic benchmarks (BP16, BP17) have $\sigma/m(1000) \sim 0.096$–$0.099$; under worst-case $x_{\rm min}$ shift they reach $\sim$0.10–0.105, marginally above the conservative selection cut of 0.1 cm$^2$/g but well below the observational limit of 0.47 cm$^2$/g (Harvey et al.).

## Appendix C: Elastic vs Momentum-Transfer Cross Section and VPM vs Born

**C.1 Elastic ≡ Transfer for Identical Majorana.** As shown in §3.2, the momentum-transfer cross section $\sigma_T^{\rm tr} \equiv \int(1-\cos\theta)\,d\sigma/d\Omega\,d\Omega$ equals the elastic cross section exactly for identical Majorana fermions. We verify this by extracting the VPM phase shifts $\{\delta_l\}$, reconstructing the full symmetrized scattering amplitude in each spin channel, and numerically integrating with the $(1-\cos\theta)$ weight using 200-point Gauss–Legendre quadrature. The singlet (even-$l$) amplitude is $f_S(\theta) + f_S(\pi-\theta) = (2/k)\sum_{l\,\text{even}} a_l P_l(\cos\theta)$, and the triplet (odd-$l$) is $f_T(\theta) - f_T(\pi-\theta) = (2/k)\sum_{l\,\text{odd}} a_l P_l(\cos\theta)$, with $a_l = (2l+1) e^{i\delta_l}\sin\delta_l$. Results:

| Point | $v$ [km/s] | $\lambda$ | $\kappa$ | $\sigma_{\rm el}/m$ [cm²/g] | $\sigma_T^{\rm tr}/\sigma_{\rm el}$ |
|:------|:----------:|:---------:|:--------:|:---------------------------:|:-----------------------------------:|
| BP1 | 12 | 1.91 | 0.037 | $4.29 \times 10^{-1}$ | 1.000000 |
| BP1 | 30 | 1.91 | 0.091 | $5.15 \times 10^{-1}$ | 1.000000 |
| BP1 | 200 | 1.91 | 0.609 | $3.56 \times 10^{-1}$ | 1.000000 |
| BP1 | 1000 | 1.91 | 3.043 | $7.23 \times 10^{-2}$ | 1.000000 |
| BP1 | 4700 | 1.91 | 14.30 | $3.72 \times 10^{-3}$ | 1.000000 |
| MAP | 12 | 48.6 | 0.170 | $1.22 \times 10^{0}$ | 1.000000 |
| MAP | 30 | 48.6 | 0.424 | $1.71 \times 10^{0}$ | 1.000000 |
| MAP | 200 | 48.6 | 2.827 | $7.24 \times 10^{-1}$ | 1.000000 |
| MAP | 1000 | 48.6 | 14.13 | $2.03 \times 10^{-1}$ | 1.000000 |
| MAP | 4700 | 48.6 | 66.43 | $2.24 \times 10^{-2}$ | 1.000000 |

The ratio is unity to within $10^{-12}$ (machine precision) across all velocities and for both $\lambda \sim 2$ and $\lambda \sim 49$, confirming the analytic argument. The elastic cross section computed via angular integration of the symmetrized differential cross section agrees with the partial-wave sum to the same precision, providing an independent consistency check. This identity holds regardless of the number of contributing partial waves — from 1 (s-wave only at $\kappa = 0.037$) to $\sim$70 (at $\kappa = 66$) — and is a direct consequence of particle indistinguishability.

**C.2 VPM vs Born approximation.** We compare the full VPM cross section against the Born approximation with Majorana symmetrization across the complete SIDM velocity range. For each velocity, the Born phase shifts are computed as $\delta_l^{\rm Born} = \kappa\lambda \int_0^\infty [j_l(\kappa x)]^2 e^{-x} x\,dx$ (analytically for $l=0$: $\delta_0^{\rm Born} = \frac{\lambda}{4\kappa}\ln(1+4\kappa^2)$, numerically via quadrature for $l>0$), and the Born cross section is built from the same Majorana-weighted partial-wave sum.

The ratio $R_{\rm VB}(v) = \sigma_{\rm VPM}/\sigma_{\rm Born}$ quantifies the nonlinear (multi-scattering) correction as a function of velocity:

| Point | $v$ [km/s] | $\lambda$ | $\kappa$ | VPM $\sigma/m$ | Born $\sigma/m$ | $R_{\rm VB}$ |
|:------|:----------:|:---------:|:--------:|:--------------:|:---------------:|:--------:|
| BP1 | 12 | 1.91 | 0.037 | $4.29 \times 10^{-1}$ | $1.87 \times 10^{0}$ | 0.23 |
| BP1 | 30 | 1.91 | 0.091 | $5.15 \times 10^{-1}$ | $1.81 \times 10^{0}$ | 0.29 |
| BP1 | 200 | 1.91 | 0.609 | $3.56 \times 10^{-1}$ | $7.42 \times 10^{-1}$ | 0.48 |
| BP1 | 1000 | 1.91 | 3.043 | $7.23 \times 10^{-2}$ | $8.48 \times 10^{-2}$ | 0.85 |
| BP1 | 4700 | 1.91 | 14.30 | $3.72 \times 10^{-3}$ | $4.05 \times 10^{-3}$ | 0.92 |
| MAP | 12 | 48.6 | 0.170 | $1.22 \times 10^{0}$ | $4.84 \times 10^{0}$ | 0.25 |
| MAP | 30 | 48.6 | 0.424 | $1.71 \times 10^{0}$ | $5.85 \times 10^{0}$ | 0.29 |
| MAP | 200 | 48.6 | 2.827 | $7.24 \times 10^{-1}$ | $7.55 \times 10^{-1}$ | 0.96 |
| MAP | 1000 | 48.6 | 14.13 | $2.03 \times 10^{-1}$ | $2.15 \times 10^{-1}$ | 0.95 |
| MAP | 4700 | 48.6 | 66.43 | $2.24 \times 10^{-2}$ | $2.38 \times 10^{-2}$ | 0.94 |

At low velocities ($\kappa \lesssim 0.5$), Born overestimates by factors of 3–5 for both benchmarks, failing to capture the non-perturbative resonance structure. At high velocities ($\kappa \gg 1$), the ratio converges toward $R_{\rm VB} \to 1$ as individual phase shifts become small and the Born expansion converges. For MAP, $R_{\rm VB}$ stabilizes at $\sim$0.94 rather than approaching unity, indicating persistent non-perturbative corrections even at $\kappa \sim 66$ — a consequence of the large coupling ($\lambda = 48.6$) maintaining non-negligible phase shifts in low-$l$ partial waves.

The Born approximation breaks down across the entire SIDM velocity range for our benchmark parameters ($\lambda \in [1.9, 48.6]$), confirming that the full VPM numerical solution is essential. The commonly-cited fitting formulas of [9], which interpolate between Born and classical limits, are unreliable at $\lambda \gtrsim 1$ — precisely the regime relevant to our parameter space.

---

## Appendix D: Maxwell–Boltzmann Velocity Averaging

**D.1 Method.** For each system with characteristic velocity $v_{\rm char}$, we compute the thermally averaged cross section

$$\langle\sigma/m\rangle_{\rm MB} = \frac{\int_0^\infty (\sigma/m)(v)\;v^3\,e^{-v^2/2v_0^2}\,dv}{\int_0^\infty v^3\,e^{-v^2/2v_0^2}\,dv}$$

with $v_0 = v_{\rm char}/\sqrt{2}$ (the 1D dispersion corresponding to the most-probable relative speed $v_{\rm char}$). The VPM solver provides $\sigma/m(v)$ at each quadrature node. We use Gauss–Legendre quadrature with $n=30$ points on $[1,\,5v_{\rm char}]$, covering $>99.9\%$ of the MB weight.

**D.2 Results.** The ratio $\langle\sigma/m\rangle_{\rm MB} / (\sigma/m)(v_{\rm char})$ depends on the local slope and curvature of $\sigma/m(v)$. Representative values for BP1 ($\lambda = 1.91$):

| System | $v_{\rm char}$ [km/s] | Fixed | MB-averaged | Ratio |
|:-------|:---------------------:|:-----:|:-----------:|:-----:|
| Draco dSph | 12 | 0.430 | 0.449 | 1.04 |
| TBTF dwarfs | 30 | 0.516 | 0.506 | 0.98 |
| NGC 2976 | 60 | 0.499 | 0.471 | 0.94 |
| NGC 720 (group) | 250 | 0.325 | 0.282 | 0.87 |
| Abell 2537 | 1100 | 0.061 | 0.052 | 0.85 |
| Bullet Cluster | 4700 | 0.004 | 0.004 | 0.98 |

At group/cluster velocities ($v \sim 250$–$1500$ km/s), where $\sigma/m(v)$ declines steeply, the MB tail samples lower cross sections, reducing the average by $\sim$12–17%. At dwarf scales ($v \lesssim 60$ km/s), the cross section is near its plateau and the correction is $\lesssim 5\%$.

**D.3 Impact on $\chi^2$.** We repeat the full $\chi^2$ analysis of §4.6 using $\langle\sigma/m\rangle_{\rm MB}$ in place of $\sigma/m(v_{\rm char})$. The $\chi^2$ increases by 2–9% (mean 6%) across all 17 relic BPs. Key figures: the worst $\chi^2/\nu$ rises from 0.65 to 0.69; the best rises from 0.42 to 0.45. All 17/17 BPs remain at $\chi^2/\nu < 1$, and no benchmark point changes viability status. The $\sim$6% shift is negligible compared to the $\sim$0.5 dex ($\times 3$) observational error bars. We conclude that the fixed-velocity approximation used in §4.5–4.6 is adequate for the current precision of SIDM constraints.

![Figure D1: BP1 comparison: fixed-velocity $\sigma/m(v)$ (blue) vs MB-averaged $\langle\sigma/m\rangle_{\rm MB}$ (red dashed). Left: cross section curves with observational data. Right: ratio $\langle\sigma/m\rangle_{\rm MB}/(\sigma/m)$ vs velocity, showing $<$20% deviations throughout.](v37_velocity_averaged.png)

---

## Appendix E: Errata and Corrections (v10.1)

The following errors were identified during peer review and corrected on 2026-03-23.

### E.1 NFW Density Normalization

The NFW scale density was previously computed as $\rho_s = \rho_{\rm crit}\,\delta_c/3$ in the prediction scripts (`predict_core_sizes.py`, `fit_sparc_baryons.py`, `sensitivity_analysis.py`). Since the characteristic overdensity is already defined as $\delta_c = (200/3)\,c^3/f(c)$ — which incorporates the factor of $1/3$ — the correct expression is:
$$\rho_s = \rho_{\rm crit}\,\delta_c$$
The spurious division by 3 reduced $\rho_s$ by a factor of 3, weakening the DM contribution to rotation curves. After correction, the best-fit mass-to-light ratios for DM-dominated dwarfs improve: DDO_154 drops from $\Upsilon_* \approx 0.44$ to $\Upsilon_* \approx 0.32$ (BP16, i.e., the 16th of the 17 relic-constrained benchmarks: $m_\chi = 29.76$ GeV, $m_\phi = 12.98$ MeV; see Table 4, rank 16), well within the physical range $0.2$–$0.8$. The $\sigma/m(v)$ cross-section computation and all other results are unaffected.

### E.2 Dead Code in Gravothermal Script

A duplicate density-conversion line in `predict_gravothermal.py` (overwritten by the subsequent line) was removed. No numerical impact.

### E.3 CSV Column Name

The scan output in `smart_scan.py` wrote a header `m_phi_GeV` while the data values were in MeV. Corrected to `m_phi_MeV` with explicit unit conversion. The actual data CSV used by downstream scripts already had correct headers.

### E.4 Coupling Ratio Convention

The dimensionless coupling ratio was displayed as $\lambda = 2\alpha m_\chi / m_\phi$ in five prediction scripts. The canonical definition is $\lambda = \alpha m_\chi / m_\phi$ (matching the VPM solver convention). The factor of 2 was only in display/print statements; all numerical computations were correct.

### E.5 Updated Prediction Results After Corrections

| Test | Result |
|------|--------|
| **Gravothermal** (8 dSphs, MAP) | 5/8 CORED, 1 FAIL, 2 ambiguous |
| **Cluster offsets** (8 systems × 3 BPs) | ALL PASS — $\sigma/m \ll$ upper limits |
| **$\Delta N_{\rm eff}$** | $\approx 0$ at BBN (Boltzmann-suppressed) |
| **SPARC+baryons** (DDO_154, BP16) | $\Upsilon_* = 0.32$ (physical range ✓) |
| **Core sizes** (DM-only, BP1/MAP) | 4/12 within $2\sigma$ (dwarfs OK, spirals need baryons) |

---

## Appendix F: Ratio Function Analyses

We present three additional ratio-function studies that illuminate the non-perturbative dynamics of the Yukawa SIDM cross section. All results use the production VPM solver with default integration parameters (verified converged to $<0.001\%$ in Appendix B).

**F.1 Attractive vs Repulsive Yukawa.** We compare $\sigma_T$ for attractive ($V = -\alpha e^{-m_\phi r}/r$) and repulsive ($V = +\alpha e^{-m_\phi r}/r$) Yukawa potentials at the same coupling strength, defining $R_{\rm AR}(v) = \sigma_{\rm att}/\sigma_{\rm rep}$.

| Point | $v$ [km/s] | $\sigma_{\rm att}/m$ | $\sigma_{\rm rep}/m$ | $R_{\rm AR}$ |
|:------|:----------:|:--------------------:|:--------------------:|:--------:|
| BP1 | 12 | 0.429 | 3.871 | 0.11 |
| BP1 | 30 | 0.515 | 29.78 | 0.017 |
| BP1 | 200 | 0.356 | 1.347 | 0.26 |
| BP1 | 1000 | 0.072 | 0.081 | 0.89 |
| BP1 | 4700 | 0.004 | 0.004 | 1.00 |
| MAP | 12 | 1.221 | 0.537 | 2.27 |
| MAP | 30 | 1.714 | 4.917 | 0.35 |
| MAP | 200 | 0.724 | 0.739 | 0.98 |
| MAP | 1000 | 0.203 | 0.203 | 1.00 |
| MAP | 4700 | 0.022 | 0.022 | 1.01 |

At high velocities both potentials converge to the Born limit ($R_{\rm AR} \to 1$). At low velocities, the attractive potential supports quasi-bound-state resonances that create dramatic departures. BP1 ($\lambda = 1.91$, just past the first critical coupling $\lambda_c = 1.68$) sits in a Ramsauer–Townsend minimum where the s-wave phase shift $\delta_0$ has crossed $\pi/2$: the resulting destructive interference suppresses $\sigma_{\rm att}$ to just $1.7\%$ of $\sigma_{\rm rep}$ at $v = 30$ km/s. Repulsive potentials have no resonances and produce monotonically decreasing cross sections — a qualitatively different SIDM phenomenology. The dramatic ratio inversion confirms that BP1 occupies a deep resonance trough, making it maximally sensitive to small coupling variations.

**F.2 Partial-Wave Decomposition.** We decompose $\sigma_T = \sum_l \sigma_l$ into individual partial-wave contributions. At $v = 30$ km/s, BP1 is completely s-wave dominated ($l=0$ carries 100\%), while MAP already activates $l=0,1$ ($f_{l=0} = 36\%, f_{l=1} = 63\%$). At $v = 1000$ km/s, BP1 requires $l = 0$–6 (7 waves $>1\%$ of total, dominant $l=1$ at 49\%), and MAP requires $l = 0$–17 (all significant). At $v = 4700$ km/s, MAP activates $\sim$70 partial waves. The transition from s-wave dominance to many-wave saturation occurs at $\kappa \sim 1$ — directly relating to the velocity at which the de Broglie wavelength matches the mediator range. The even/odd ratio converges to the Majorana spin-statistics limit of 3:1 at high $v$ (see §7.7).

**F.3 VPM vs Born.** Full results comparing the VPM and Born cross sections are presented in Appendix C.2. The key finding is that the Born approximation overestimates by factors of 3–5 at low velocities for both benchmarks, confirming that the full non-perturbative VPM solution is essential across the entire SIDM velocity range.
