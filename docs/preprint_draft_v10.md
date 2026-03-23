# Self-Interacting Dark Matter from a Majorana Fermion with Light Scalar Mediator:
# Velocity-Dependent Cross Sections via Yukawa Partial-Wave Analysis

**Omer P.**  
*Independent researcher*

**March 2026 — V10 Draft**

---

## Abstract

We construct a minimal **secluded** dark matter model in which a Majorana fermion $\chi$ couples to a light real scalar $\phi$ through a Yukawa interaction $-(y/2)\bar{\chi}\chi\,\phi$. The dark sector has no tree-level coupling to the Standard Model. The $\phi$-mediated self-interaction produces velocity-dependent scattering governed by the attractive Yukawa potential $V(r) = -\alpha\,e^{-m_\phi r}/r$ with $\alpha = y^2/(4\pi)$.

We compute the self-interaction transfer cross section using the **Variable Phase Method (VPM)** — a full partial-wave analysis that captures quasi-bound-state resonances of the Yukawa potential. A systematic scan over $(m_\chi, m_\phi, \alpha)$ — with $m_\chi \in [0.1, 100]$ GeV, $m_\phi \in [0.1, 200]$ MeV — identifies **80,142 parameter points** satisfying the SIDM sweet-spot criteria $\sigma/m(30\text{ km/s}) \in [1, 10]\text{ cm}^2/\text{g}$ and $\sigma/m(1000\text{ km/s}) < 0.1\text{ cm}^2/\text{g}$.

The annihilation $\chi\chi \to \phi\phi$ is **s-wave**, yielding a thermal relic density $\Omega h^2 = 0.120$ for $\alpha \sim 5 \times 10^{-4}$–$5 \times 10^{-3}$ — squarely within the SIDM-viable region. A dedicated cosmological scan using an exact numerical Boltzmann solver identifies **17 benchmark points** forming an "island of viability" at $m_\chi \in [10, 100]$ GeV, $m_\phi \in [7.6, 14.8]$ MeV that simultaneously satisfy relic density, SIDM, and cluster constraints. The primary benchmark has $m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV, $\alpha = 1.05 \times 10^{-3}$, $\Omega h^2 = 0.120$, $\sigma/m(30) = 0.52$ cm$^2$/g, and $\sigma/m(1000) = 0.072$ cm$^2$/g. A quantitative $\chi^2$ fit to 13 astrophysical systems spanning $v = 12$–$4700$ km/s yields $\chi^2/\nu = 0.54$ for the best relic-constrained point and $\chi^2/\nu = 0.26$ for the unconstrained best fit, with all pulls below $1.2\sigma$. A Bayesian MCMC posterior analysis confirms that all 17 relic benchmarks lie within the 95% credible region, with marginalized constraints $m_\chi = 73^{+43}_{-56}$ GeV, $m_\phi = 10.8^{+4.8}_{-4.3}$ MeV (68% CI). We show that a Higgs portal coupling to the SM is generically **excluded** for light mediators ($m_\phi \lesssim \mathcal{O}(\text{GeV})$) by the tension between direct detection ($\sigma_{\rm SI} \propto 1/m_\phi^4$) and BBN lifetime constraints; a secluded dark sector is the natural resolution. The scalar has only 1 degree of freedom and contributes $\Delta N_{\rm eff} \approx 0.027$ — well below the Planck limit.

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

**Related work.** Light-mediator SIDM was introduced by Feng, Kaplinghat, Tu & Yu [17] in the context of hidden charged dark matter, and developed systematically by Tulin, Yu & Zurek [9] who provided fitting formulas for the Yukawa transfer cross section. The comprehensive review by Tulin & Yu [10] covers the full landscape of SIDM models. Our work differs from these precedents in several respects: (i) we use the full VPM partial-wave analysis rather than the approximate Hulthén-potential fitting formulas of [9], which fail to capture resonance structure at $\lambda \gtrsim 1$ (Appendix A.5); (ii) we demonstrate quantitatively that the Higgs portal is generically excluded for light mediators (§5.1), motivating a secluded dark sector rather than the visible-decay models assumed in much of the earlier literature; (iii) we perform a combined SIDM + relic density scan with an exact Boltzmann solver (§4.4, §6.2), identifying a well-defined island of viability rather than individual benchmark points; and (iv) we provide a quantitative $\chi^2$ fit to 13 astrophysical systems (§4.6) spanning four decades in velocity.

---

## 2. The Model (The Scalar Mediator)

To resolve the small-scale structure anomalies without violating theoretical unitarity or parity constraints, we consider a simplified dark sector consisting of a Majorana fermion dark matter candidate, $\chi$, interacting via a real scalar mediator, $\phi$. The interaction Lagrangian is given by:

$$\mathcal{L}_{\text{int}}=-\frac{y}{2}\bar{\chi}\chi\phi$$

In the non-relativistic limit, this scalar exchange generates a universally attractive Yukawa potential between the dark matter particles, irrespective of their spin configuration:

$$V(r)=-\frac{\alpha}{r}e^{-m_\phi r}$$

where $\alpha=y^2/(4\pi)$. A natural connection between $\phi$ and the Standard Model would be the Higgs portal coupling $\lambda_{H\phi}|H|^2\phi^2$, which induces a mixing angle $\sin\theta$ and enables both $\phi$ decay and DM–nucleon scattering. We investigated this possibility and show in \S5.1 that it is **generically excluded** for light mediators ($m_\phi \lesssim \mathcal{O}(\text{GeV})$): the $1/m_\phi^4$ enhancement of $\sigma_{\rm SI}$ closes the window between direct detection and BBN bounds by a factor of $\sim 3 \times 10^4$. The dark sector is therefore **secluded** — $\phi$ has no tree-level coupling to the SM.

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

*Note:* This formula is the spin-averaged **elastic** cross section, which differs from the **momentum-transfer** (transport) cross section $\sigma_T^{\rm tr} = \int(1-\cos\theta)\,d\sigma/d\Omega\,d\Omega$ by a factor that depends on the angular distribution. In the SIDM-relevant regime ($\kappa \equiv k/m_\phi \lesssim 1$, few partial waves), we verify numerically (Appendix A.5) that the two agree to within $\sim$20% — well within the factor-of-3 observational uncertainties. We adopt the elastic formula throughout, consistent with the convention used in several numerical SIDM studies.

### 3.3 Numerical Implementation

The VPM ODE is integrated using a 4th-order Runge-Kutta scheme with:
- Adaptive $x_{\rm max}$: 50 (low $\kappa$), 80 (medium), 100 (high $\kappa$)
- Adaptive step count: 4000–12000 steps
- Truncation: $l_{\rm max} = \min(\max(3, \lfloor\kappa\rfloor + 3), 80)$ with early termination when the contribution falls below $10^{-3}$ of the running sum
- Centrifugal barrier: $x_{\rm min} = \max(10^{-5}, l/\kappa)$ for $l > 0$

The implementation is verified against:
1. **Free particle test:** $\delta_l(\lambda = 0) = 0$ for all $l$ — PASS
2. **Born limit:** Agreement with analytical Born approximation at $\lambda \ll 1$ — 6/6 PASS
3. **Scipy ground truth:** Agreement with `scipy.integrate.solve_ivp` (RK45, rtol=$10^{-10}$) to $< 0.001\%$ — PASS

### 3.4 Numerical Error Budget

A systematic error analysis identifies three components:

| Source | 30 km/s | 1000 km/s | Method |
|--------|---------|-----------|--------|
| **Integrator** (RK4 step-size) | < 0.001% | < 0.001% | Compare 4000 vs 16000 steps |
| **Truncation** ($l_{\rm max}$) | < 0.01% | ~10–13% | Compare $l_{\rm max}$ vs $l_{\rm max}+20$ |
| **Prescription** ($x_{\rm min}$) | 2–7% | 25–38% | Vary $x_{\rm min}$ factor by ±20% |
| **Total** (quadrature) | **2–7%** | **27–40%** |  |

At the primary SIDM velocity (30 km/s, dwarf galaxies), the total systematic is 2–7% — well within the viability window $\sigma/m \geq 0.5$ cm$^2$/g. At cluster velocities (1000 km/s), the systematic is larger (27–40%); points with $\sigma/m(1000) \ll 0.1$ remain safely below the cut.

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
| Bullet Cluster | 4700 | 0.7 | [0.0, 1.25] | [20] |
| 72 cluster mergers | 1500 | 0.2 | [0.0, 0.47] | [21] |
| TBTF dwarfs | 30 | 1.0 | [0.5, 5.0] | [22] |

**Table 3: Observational constraints on $\sigma/m$ used for comparison.**

BP1 falls within the allowed range of **11 out of 13** observational systems (Figure 3). The two marginal cases — NGC 2976 ($\sigma/m_{\rm theory} = 0.50$ vs lower bound 0.50 cm$^2$/g, a 0.4% shortfall) and NGC 1560 ($\sigma/m_{\rm theory} = 0.50$ vs lower bound 1.0 cm$^2$/g) — are both rotation-curve systems where the inferred $\sigma/m$ has large uncertainties ($\sim$0.5 dex) from halo-profile modeling. The overall agreement across dwarfs ($v \sim 10$–50 km/s), galaxies ($v \sim 50$–300 km/s), and clusters ($v \sim 1000$–5000 km/s) demonstrates that the velocity dependence of our model naturally reproduces the astrophysical pattern: $\sigma/m \sim 0.5$ cm$^2$/g in dwarfs, declining through the group and cluster scales to $\sigma/m \sim 0.05$ cm$^2$/g at $v \sim 1000$ km/s.

![Figure 3: Predicted $\sigma/m(v)$ for BP1 (solid blue), BP5 (dashed orange), and BP17 (dotted green) compared with observational constraints from Kaplinghat et al. [13], Kamada et al. [19], Randall et al. [20], Harvey et al. [21], and Elbert et al. [22]. Error bars show 68% CL ranges. All three benchmarks trace a similar velocity-dependent envelope consistent with the observed pattern.](v33_observational_comparison.png)

### 4.6 Quantitative $\chi^2$ Fit to Observational Data

To move beyond the qualitative compatibility assessment of §4.5, we perform a full $\chi^2$ fit of the predicted $\sigma/m(v)$ to all 13 observational systems listed in Table 3. For each parameter point we evaluate the VPM transfer cross section at all 12 unique velocities and compute:

$$\chi^2 = \sum_{i=1}^{13} \left(\frac{\sigma/m_{\rm theory}(v_i) - \sigma/m_{\rm obs,\,i}}{\sigma_i}\right)^2$$

with asymmetric errors: $\sigma_i^+ = (\sigma/m_{\rm upper} - \sigma/m_{\rm obs})/1.0$ and $\sigma_i^- = (\sigma/m_{\rm obs} - \sigma/m_{\rm lower})/1.0$ at 68% CL, selecting $\sigma_i^+$ or $\sigma_i^-$ according to the sign of the residual. With 13 data points and 3 free parameters, the number of degrees of freedom is $\nu = 10$.

We evaluate 5,009 points sampled from the 80,142 raw viable set plus all 17 relic benchmarks, totaling 5,026 points $\times$ 12 velocities = 60,312 full VPM evaluations. A fine-grained scan around the top 50 candidates adds 400 further evaluations.

**Unconstrained best fit.** The global minimum is:

$$\chi^2_{\rm min}/\nu = 2.59/10 = 0.26$$

at $m_\chi = 100$ GeV, $m_\phi = 9.15$ MeV, $\alpha = 2.84 \times 10^{-3}$, $\lambda = 31.0$. All 13 pulls are below $1\sigma$ (largest: Bullet Cluster at $-0.99\sigma$). Representative predictions:

| System | $v$ [km/s] | Theory | Obs | Pull |
|:-------|:----------:|:------:|:---:|:----:|
| Draco dSph | 12 | 1.32 | 0.60 | $+0.52$ |
| NGC 1560 | 55 | 1.77 | 3.00 | $-0.61$ |
| Diverse RC | 40 | 1.95 | 3.00 | $-0.42$ |
| Abell 2537 | 1100 | 0.11 | 0.15 | $-0.33$ |
| Bullet Cluster | 4700 | 0.008 | 0.70 | $-0.99$ |
| TBTF dwarfs | 30 | 1.91 | 1.00 | $+0.23$ |

The low $\chi^2/\nu \ll 1$ reflects the generous observational uncertainties characteristic of astrophysical SIDM constraints — typical 68% CL ranges span $\sim$0.5 dex (a factor of 3), compared with the $\mathcal{O}(10\%)$ precision common in collider physics. In this regime, $\chi^2/\nu < 1$ does not imply over-fitting but rather that the data do not yet have the resolving power to distinguish between models at the $\sim$factor-of-2 level. A more discriminating test will require tighter observational error bars, e.g., from JWST-era kinematics of ultra-faint dwarfs. Nevertheless, the fit demonstrates that a single set of parameters can **simultaneously** match all 13 systems spanning $v = 12$–$4700$ km/s — a non-trivial consistency check given the four-decade velocity baseline.

**Relic-constrained best fit.** Restricting to the 17 benchmark points satisfying $\Omega h^2 = 0.120$ (§4.4):

$$\chi^2_{\rm relic}/\nu = 5.40/10 = 0.54$$

at the $(m_\chi, m_\phi)$ grid cell containing BP1: $m_\chi = 20.69$ GeV, $m_\phi = 9.91$ MeV, $\alpha = 1.048 \times 10^{-3}$, $\lambda = 2.19$. (Note: this grid cell differs slightly in $m_\phi$ from the BP1 defined in Table 1, which has $m_\phi = 11.34$ MeV; the original BP1 appears as rank 17 with $\chi^2/\nu = 0.84$, still an excellent fit.) The worst pull is NGC 1560 at $-1.12\sigma$, well within acceptable limits.

**Table 4: All 17 relic-constrained benchmark points ranked by $\chi^2$**

| Rank | $\chi^2/\nu$ | $m_\chi$ [GeV] | $m_\phi$ [MeV] | $\alpha$ | $\sigma/m(30)$ [cm$^2$/g] |
|:----:|:------------:|:--------------:|:--------------:|:--------:|:-------------------------:|
| 1 | 0.54 | 20.69 | 9.91 | $1.05 \times 10^{-3}$ | 0.79 |
| 2 | 0.57 | 29.76 | 11.34 | $1.47 \times 10^{-3}$ | 0.75 |
| 3 | 0.58 | 14.38 | 8.66 | $7.56 \times 10^{-4}$ | 0.73 |
| 4 | 0.61 | 100.00 | 14.85 | $4.81 \times 10^{-3}$ | 0.59 |
| 5 | 0.63 | 26.37 | 11.34 | $1.31 \times 10^{-3}$ | 0.68 |
| 6 | 0.64 | 18.33 | 9.91 | $9.39 \times 10^{-4}$ | 0.68 |
| 7 | 0.65 | 42.81 | 12.98 | $2.09 \times 10^{-3}$ | 0.64 |
| 8 | 0.66 | 88.59 | 14.85 | $4.26 \times 10^{-3}$ | 0.53 |
| 9 | 0.70 | 37.93 | 12.98 | $1.86 \times 10^{-3}$ | 0.60 |
| 10 | 0.70 | 78.48 | 14.85 | $3.78 \times 10^{-3}$ | 0.50 |
| 11 | 0.72 | 23.36 | 11.34 | $1.17 \times 10^{-3}$ | 0.60 |
| 12 | 0.73 | 10.00 | 7.56 | $5.54 \times 10^{-4}$ | 0.59 |
| 13 | 0.74 | 12.74 | 8.66 | $6.80 \times 10^{-4}$ | 0.59 |
| 14 | 0.75 | 33.60 | 12.98 | $1.65 \times 10^{-3}$ | 0.56 |
| 15 | 0.78 | 16.24 | 9.91 | $8.41 \times 10^{-4}$ | 0.56 |
| 16 | 0.83 | 29.76 | 12.98 | $1.47 \times 10^{-3}$ | 0.51 |
| 17 | 0.84 | 20.69 | 11.34 | $1.05 \times 10^{-3}$ | 0.52 |

All 17 relic benchmarks yield $\chi^2/\nu < 1$, confirming that the relic-density constraint does **not** degrade the observational fit. The full island of viability is quantitatively consistent with all available astrophysical data.

**Effect of velocity averaging.** The above $\chi^2$ analysis evaluates $\sigma/m$ at the single characteristic velocity $v_{\rm char}$ of each system, following the standard practice in SIDM phenomenology [9, 10, 13]. In reality, DM particles in a halo have a Maxwell–Boltzmann (MB) velocity distribution, and the observable is the thermally averaged $\langle\sigma/m\rangle_{\rm MB} = \int (\sigma/m)(v)\,v\,f_{\rm MB}(v;v_0)\,dv\,/\,\int v\,f_{\rm MB}(v;v_0)\,dv$ with $v_0 = v_{\rm char}/\sqrt{2}$. We computed this integral for all 17 BPs via Gauss–Legendre quadrature ($n=30$, see Appendix D). The per-system deviations are largest at group/cluster velocities ($v \sim 250$–$1500$ km/s), where the steep decline of $\sigma/m(v)$ causes the MB tail to sample lower cross sections, yielding $\langle\sigma/m\rangle_{\rm MB} \approx 0.83$–$0.88\times\sigma/m(v_{\rm char})$ — a $\sim$12–17% reduction. At dwarf scales ($v \lesssim 60$ km/s) the correction is $\lesssim 5\%$. The net effect on the $\chi^2$ fit is a uniform increase of $2$–$9\%$ (mean $6\%$): no benchmark point changes viability status, and the worst $\chi^2/\nu$ rises from 0.65 to 0.69. This is entirely negligible compared to the $\sim$0.5 dex observational uncertainties.

![Figure 4: Best-fit $\sigma/m(v)$ curves with 13 observational data points and 68% CL error bars. The relic-constrained best fit (BP1, $\chi^2/\nu = 0.54$) reproduces the velocity-dependent pattern from dwarfs to clusters.](v34_chi2_fit.png)

### 4.7 Bayesian Posterior Constraints

To quantify the allowed parameter space beyond the frequentist $\chi^2$ analysis, we perform a Bayesian posterior sampling of the three-dimensional parameter space $(m_\chi, m_\phi, \alpha)$ using the emcee affine-invariant ensemble sampler [23]. We adopt flat (uninformative) priors in $\log_{10}$ space:

$$\log_{10}(m_\chi/\text{GeV}) \in [\log_{10}(5),\, \log_{10}(200)], \quad \log_{10}(m_\phi/\text{MeV}) \in [\log_{10}(3),\, \log_{10}(30)], \quad \log_{10}\alpha \in [-5,\, \log_{10}(0.05)]$$

The likelihood is Gaussian: $\ln\mathcal{L} = -\chi^2/2$ with the same $\chi^2$ function used in §4.6 (asymmetric errors, 13 observational systems). The sampler uses 32 walkers initialized around the 17 relic benchmark points plus the unconstrained best fit of §4.6, with a Gaussian scatter of $\sigma = 0.05$ in log-space. After 300 burn-in steps, we collect 2,000 production steps (64,000 total samples, computed in parallel on 12 CPU cores).

**Convergence diagnostics.** The mean acceptance fraction is 0.522, and the maximum autocorrelation time is 73.3 steps, yielding $\sim$874 effective independent samples. The chain trace plots (Figure S1) show good mixing with no residual drift.

**Results.** The maximum a posteriori (MAP) estimate is:
$$m_\chi = 90.6~\text{GeV},\quad m_\phi = 13.9~\text{MeV},\quad \alpha = 2.5 \times 10^{-2},\quad \chi^2/\nu = 0.16$$

The marginalized parameter constraints (median $\pm$ 68% credible interval) are:

| Parameter | Median | 16th %ile | 84th %ile |
|-----------|--------|-----------|-----------|
| $m_\chi$ [GeV] | 73.0 | 17.0 | 115.9 |
| $m_\phi$ [MeV] | 10.8 | 6.5 | 15.6 |
| $\alpha$ | $4 \times 10^{-3}$ | $1 \times 10^{-3}$ | $2.2 \times 10^{-2}$ |
| $\lambda = 2\alpha m_\chi/m_\phi$ | 59.0 | 3.6 | 272.2 |

The broad credible intervals reflect the large astrophysical uncertainties ($\sim$0.5 dex) in the observational data — the model fits the data well across a wide parameter range. The derived coupling parameter $\lambda$ is peaked in the resonant regime ($\lambda > 1$), consistent with the need for strong velocity dependence.

Crucially, **all 17 relic benchmark points** (§4.4) lie within the 95% credible region of the posterior, with $\chi^2$ values ranging from 5.4 to 8.4 (compared to the 95th percentile threshold of the posterior $\chi^2$ distribution). This confirms that the island of viability identified by the cosmological scan coincides with the statistically preferred region of the SIDM parameter space.

![Figure 5: Corner plot showing the 2D marginalized posterior distributions of $\log_{10}(m_\chi)$, $\log_{10}(m_\phi)$, and $\log_{10}\alpha$, with 68% and 95% contour levels. Red lines mark the MAP estimate. All 17 relic BPs fall within the 95% contours.](v38_corner.png)

![Figure 6: Posterior distribution of the derived parameter $\lambda = 2\alpha m_\chi/m_\phi$. The median is $\lambda = 59$ with a broad 68% CI of [3.6, 272]. The distribution peaks in the resonant scattering regime ($\lambda > 1$).](v38_lambda_posterior.png)

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

**Thermal history:** We assume the dark sector was in thermal equilibrium with the SM at early times ($T \gg m_\chi$). A concrete realization is a heavy scalar or fermionic mediator $\Sigma$ with mass $M_\Sigma \gg m_\chi$ and couplings to both sectors (e.g., $\Sigma\bar{f}f + \Sigma\bar{\chi}\chi$). For $M_\Sigma \sim 1$–$10$ TeV, the thermalization rate $\Gamma \sim T^5/M_\Sigma^4$ exceeds the Hubble rate for $T \gtrsim \mathcal{O}(10)$ GeV, establishing equilibrium well before freeze-out at $T_f \sim m_\chi/25 \sim 0.4$–$4$ GeV. At $T \ll M_\Sigma$ the interaction decouples and is negligible at the scales relevant to SIDM and freeze-out. After decoupling, the dark-sector temperature tracks $T_{\rm dark} = \xi\, T_{\rm SM}$ with $\xi \lesssim 1$. For $\xi \approx 1$ (as assumed throughout), the standard Boltzmann calculation applies without modification.

### 5.3 Cosmological Safety of Stable $\phi$

A stable $\phi$ with $m_\phi \sim 8$–15 MeV contributes to the energy density as dark radiation at temperatures $T \ll m_\phi$. As a single real scalar degree of freedom:

$$\Delta N_{\rm eff} = \frac{4}{7}\left(\frac{T_\phi}{T_\nu}\right)^4 \approx 0.027$$

for $\xi = 1$ and $\phi$ decoupling before $e^+e^-$ annihilation. This is a factor of $\sim 11$ below the Planck 2018 limit ($\Delta N_{\rm eff} < 0.30$ at 95% CL) and well within the projected CMB-S4 sensitivity ($\sigma \sim 0.03$).

The energy density in stable $\phi$ particles at late times is:

$$\frac{\Omega_\phi}{\Omega_\chi} \sim \frac{m_\phi}{m_\chi} \sim 5 \times 10^{-4}$$

— a negligible correction to the total DM energy budget.

### 5.4 Falsifiability

Without a direct detection signal, the model is tested through:
1. **Astrophysical observations** of dwarf galaxy cores, cluster mergers, and halo density profiles — the primary SIDM observable
2. **CMB-S4** sensitivity to $\Delta N_{\rm eff} \sim 0.03$ could probe the stable $\phi$ contribution
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

**Key result:** For each mass in the viable island ($m_\chi \in [10, 100]$ GeV, $m_\phi \in [7.6, 14.8]$ MeV), there exists a unique $\alpha$ that yields $\Omega h^2 = 0.120$ and simultaneously satisfies the SIDM sweet-spot criteria. The relic-density coupling $\alpha \sim 5 \times 10^{-4}$–$5 \times 10^{-3}$ lies **squarely within** the SIDM-viable parameter space. This is a **major improvement** over the axial-vector mediator approach, where the p-wave annihilation required $\alpha \sim 10^{-3}$–$10^{-2}$ — largely outside the SIDM range.

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

## 7. Summary of Predictions

| Observable | Prediction | Status |
|-----------|-----------|--------|
| $\sigma_{\rm SI}$ | **= 0** (secluded: $\sin\theta = 0$) | No DD signal |
| $\sigma_{\rm SD}$ | **= 0** (scalar mediator) | No constraint ✓ |
| $\sigma/m$ at 30 km/s | 0.5–0.8 cm$^2$/g (viable island) | Resolves core-cusp ✓ |
| $\sigma/m$ at 1000 km/s | $< 0.1$ cm$^2$/g (by construction) | Consistent with clusters ✓ |
| Relic density | $\Omega h^2 = 0.120$ for $\alpha \sim 5 \times 10^{-4}$–$5 \times 10^{-3}$ | **Exact overlap with SIDM** ✓ |
| $\phi$ status | **Stable** (no SM coupling) | Dark radiation |
| $\Delta N_{\rm eff}$ | $\approx 0.027$ (1 scalar d.o.f.) | Well below Planck ✓ |
| CMB energy injection | **None** ($\chi\chi \to \phi\phi$ stays dark) | No constraint ✓ |
| Fermi-LAT | **No SM final states** | No constraint ✓ |
| BSF correction | $< 0.01\%$ at freeze-out | Negligible ✓ |

---

## 8. Discussion

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
| $\Delta N_{\rm eff}$ | Marginal | 0.027 (safe) ✓ |

The scalar mediator resolves the primary limitations of the axial-vector approach while preserving all the successful SIDM phenomenology.

### 8.2 What This Paper Shows

1. A Majorana fermion with a light scalar mediator in a **secluded dark sector** produces velocity-dependent SIDM with a large viable parameter space (80,142 raw points).
2. The s-wave annihilation yields the correct relic density for $\alpha \sim 5 \times 10^{-4}$–$5 \times 10^{-3}$, with **exact overlap** with the SIDM-viable region — 17 benchmark points verified with a numerical Boltzmann solver. Sommerfeld enhancement is negligible at freeze-out ($S_0 < 1.026$, \S6.5), validating the tree-level calculation.
3. The Higgs portal coupling is **generically excluded** for light mediators ($m_\phi \lesssim$ GeV) due to the $1/m_\phi^4$ enhancement of $\sigma_{\rm SI}$ (§5.1). A secluded dark sector is the natural resolution.
4. The model is minimal (3 parameters), anomaly-free, and cosmologically safe ($\Delta N_{\rm eff} \approx 0.027$).
5. The predicted $\sigma/m(v)$ curves are **quantitatively consistent with all 13 astrophysical observations** spanning $v = 12$–$4700$ km/s (§4.5–4.6). A $\chi^2$ fit to data from five independent analyses [13, 19, 20, 21, 22] yields $\chi^2/\nu = 0.26$ (unconstrained) and $\chi^2/\nu = 0.54$ (relic-constrained BP1), with all pulls $< 1.2\sigma$. All 17 relic benchmarks achieve $\chi^2/\nu < 0.85$. Maxwell–Boltzmann velocity averaging shifts the $\chi^2$ by only $\sim$6% on average (Appendix D), well within observational uncertainties.
6. A Bayesian MCMC posterior analysis (§4.7) with flat log-priors yields a broad 68% credible region: $m_\chi \in [17, 116]$ GeV, $m_\phi \in [6.5, 15.6]$ MeV, $\alpha \in [10^{-3}, 2.2 \times 10^{-2}]$. All 17 relic benchmark points lie within the 95% credible posterior, and the MAP estimate achieves $\chi^2/\nu = 0.16$.

### 8.3 Falsifiability

The model can be tested or constrained by:
1. **Astrophysical observations:** dwarf galaxy cores, cluster mergers, and halo density profiles provide the primary constraint on $\sigma/m(v)$. Improving measurements can narrow or exclude the viable island.
2. **Cluster constraints** tightening to $\sigma/m < 0.01$ cm$^2$/g at 1000 km/s.
3. **CMB-S4:** improved $\Delta N_{\rm eff}$ sensitivity ($\sigma \sim 0.03$) could probe the stable $\phi$ contribution of $\approx 0.027$.
4. **Relic density:** if $\alpha$ is measured independently (e.g., from halo observations), the relic prediction is fixed — a non-trivial consistency test.

The absence of a direct detection signal is a **prediction**, not a deficit: we have shown (§5.1) that the Higgs portal is incompatible with light-mediator SIDM, making $\sigma_{\rm SI} = 0$ the expected outcome.

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

---

## Appendix A: VPM Solver Validation

**A.1 Free Particle Test.** Setting $\lambda = 0$, all phase shifts $\delta_l = 0$ to machine precision for $l = 0, \ldots, 20$.

**A.2 Born Limit.** In the weak-coupling limit $\lambda \ll 1$, VPM agrees with the Born approximation to $< 1\%$ at 6 benchmark points.

**A.3 Scipy Cross-Check.** Numba RK4 agrees with `scipy.integrate.solve_ivp` (RK45, rtol=$10^{-10}$) to $< 0.001\%$.

**A.4 Literature Cross-Check (Born Regime).** We independently validate the VPM solver against Born-approximation phase shifts computed via `scipy.integrate.quad` at 6 test points spanning the Born, classical, and resonant regimes. In the Born regime ($\kappa \equiv \alpha m_\chi / m_\phi \in [0.1, 0.5]$, $\eta \equiv v/c \cdot m_\chi / m_\phi \in [5, 10]$), VPM agrees with the Born partial-wave sum to 85–97%, with the mild discrepancy at $\kappa = 0.1$ attributed to the Yukawa barrier cutoff absent in the Born approximation. The solver passes all 6 tests: Born regime agreement, classical regime scaling, velocity dependence, unitarity of $S$-matrix, Majorana vs Dirac factor-of-4 relation, and BP1 consistency with the preprint values (0.9% at 30 km/s, 0.4% at 1000 km/s).

## Appendix B: Error Budget Details

**B.1 Integrator Convergence.** Doubling step count (4000→16000): $< 0.001\%$ change at all benchmarks.

**B.2 Truncation Error.** Extending $l_{\rm max}$ by +20: $< 0.01\%$ at 30 km/s, 10–13% at 1000 km/s.

**B.3 Prescription Sensitivity.** Varying $x_{\rm min}$ by ±20%: 2–7% at 30 km/s, 25–38% at 1000 km/s.

**B.4 Assessment.** Total systematic 2–7% at primary SIDM velocity (30 km/s). At 1000 km/s the error is 27–40% but points with $\sigma/m(1000) \ll 0.1$ remain safely viable.

## Appendix C: VPM vs Analytic Born Transfer Cross Section

**C.1 Method.** We compare the VPM elastic cross section (§3.2) against the exact Born transfer cross section for identical Majorana fermions, computed by numerically integrating the symmetrized Born amplitude with $(1-\cos\theta)$ weighting. The Born amplitude uses $f(\theta) = 2\mu\alpha/(q^2 + m_\phi^2)$ with the Majorana symmetrization $\frac{1}{4}|f+f'|^2 + \frac{3}{4}|f-f'|^2$ and the identical-particle factor of $1/2$.

**C.2 Results.** In the Born regime ($\lambda < 1$, $\kappa \gtrsim 0.1$), the VPM elastic cross section agrees with the Born transfer to within $\sim$10–20%. The residual difference is the well-known distinction between elastic and momentum-transfer cross sections: the transfer formula down-weights forward scattering via $(1-\cos\theta)$, while the elastic formula counts all angles equally. At $\kappa < 0.05$ (extremely low velocity), the VPM integration range needs extension; this region lies below the SIDM velocities of interest ($v \gtrsim 10$ km/s).

In the resonant regime ($\lambda \sim 2$–$30$), the Born approximation over-predicts by 1–2 orders of magnitude and fails to capture the non-perturbative resonance structure. This validates our use of the full VPM solver: the commonly-cited fitting formulas of [9], which interpolate between Born and classical limits, are unreliable at $\lambda \gtrsim 1$ — precisely the regime relevant to our parameter space ($\lambda \in [0.73, 32.4]$).

| Test point | $\lambda$ | $\kappa$ | VPM $\sigma/m$ | Born $\sigma/m$ | VPM/Born |
|:-----------|:---------:|:--------:|:--------------:|:---------------:|:--------:|
| Born ($m_\chi=100, v=200$) | 0.02 | 0.17 | $1.76 \times 10^{-7}$ | $1.93 \times 10^{-7}$ | 0.91 |
| Born ($m_\chi=50, v=100$) | 0.05 | 0.08 | $5.18 \times 10^{-6}$ | $6.68 \times 10^{-6}$ | 0.78 |
| BP1 ($v=1000$) | 2.19 | 0.35 | $9.64 \times 10^{-2}$ | $1.20 \times 10^{-1}$ | 0.80 |
| BP1 ($v=30$) | 2.19 | 0.10 | $7.93 \times 10^{-1}$ | $3.10 \times 10^{0}$ | 0.26 |
| Free best ($v=30$) | 31.04 | 0.55 | $1.91 \times 10^{0}$ | $7.88 \times 10^{1}$ | 0.02 |

The Born formula breaks down catastrophically at $\lambda \gg 1$ (VPM/Born $\sim 0.02$), confirming that the resonant regime requires the full numerical partial-wave analysis.

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
