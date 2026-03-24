# Benchmark Selection Methodology — Transparency Document

## Overview

This document explains how the 17 benchmark points in this analysis are selected, why the process is deterministic and free of observational bias, and how it relates to standard practices in the SIDM literature.

---

## The Selection Pipeline

The benchmark points are produced by a **two-stage deterministic pipeline** with zero free parameters at the prediction stage:

### Stage 1: Relic Density Constraint → Fixes α

For each cell in a 20×30 grid over $(m_\chi, m_\phi)$ with $m_\chi \in [10, 100]$ GeV and $m_\phi \in [1, 50]$ MeV:

1. A numerical Boltzmann solver (4th-order Runge–Kutta) computes the relic density $\Omega h^2$ as a function of $\alpha$.
2. Bisection in $\alpha$ finds the **unique** coupling that yields $\Omega h^2 = 0.120$ (Planck 2018 central value).
3. This step is purely cosmological — **no astrophysical SIDM data enters**.

**Result:** Each $(m_\chi, m_\phi)$ cell maps to exactly one $\alpha_{\rm relic}$. There is no choice or optimization involved.

### Stage 2: SIDM Viability Filter

For each relic-constrained point $(m_\chi, m_\phi, \alpha_{\rm relic})$:

1. The VPM (Variable Phase Method) solver computes $\sigma_T/m$ at two reference velocities.
2. Points are kept if they satisfy:
   - $\sigma/m(30\text{ km/s}) \in [0.5, 10]$ cm²/g — dwarf-scale self-interaction
   - $\sigma/m(1000\text{ km/s}) < 0.1$ cm²/g — cluster-scale safety

These cuts come directly from the observational literature (see §References below), not from our analysis.

**Result:** 17 out of 600 grid cells survive. These are the benchmark points.

### Stage 3: Comparison with Observations (Post-hoc)

Only **after** the 17 points are fixed do we compute $\chi^2$ against 13 astrophysical systems. The $\chi^2$ is a **diagnostic**, not a selection criterion. No point is added, removed, or modified based on $\chi^2$.

---

## Why This Is Not "Cooking" Points

| Concern | Answer |
|---------|--------|
| "You picked the best-fitting points" | No. The 17 points are **all** points that satisfy relic + SIDM. We report all of them, not a subset. |
| "The cuts were tuned to get good fits" | No. The SIDM cuts (0.5–10 at dwarfs, <0.1 at clusters) are taken from Kaplinghat, Tulin & Yu (2016) and are standard in the field. |
| "$\chi^2/\nu = 0.38$ looks suspiciously good" | This is the **best** of 17 points (trial factor). The median $\chi^2/\nu$ across all 17 BPs should also be reported. Two of 13 systems are one-sided upper limits contributing $\chi^2 = 0$, and several systems have error bars spanning 0.5 dex. |
| "You have hidden free parameters" | At the prediction stage, $(m_\chi, m_\phi, \alpha)$ are fully determined by relic density + SIDM. The velocity-dependent $\sigma/m(v)$ curve is a **zero-parameter prediction**. |

---

## Why the Process Is Deterministic

Given:
- The Lagrangian $\mathcal{L} = \frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\,\phi + \frac{1}{2}m_\phi^2\phi^2$
- The grid resolution (20×30)
- The relic target $\Omega h^2 = 0.120$
- The SIDM cuts from the literature

→ The 17 benchmark points are **uniquely determined**. Running the pipeline again with the same inputs produces identical output. There is no random sampling, no MCMC initialization dependence, and no human intervention in the selection.

---

## Who Else Uses This Approach

This "cosmology-driven parameter selection" is the standard methodology in the SIDM literature:

| Reference | Method | Citation count |
|-----------|--------|---------------|
| Tulin, Yu & Zurek (2013) | Fix $(m_\chi, m_\phi)$ → compute $\sigma/m$ from Yukawa → compare to observations | ~900 |
| Kaplinghat, Tulin & Yu (2016) | SIDM-viable parameter space → $\chi^2$ fit to 7 astrophysical systems | ~800 |
| Kamada, Kaplinghat, Pace & Yu (2017) | SIDM model parameters → fit to diverse rotation curves | ~400 |
| Ren, Kwa, Kaplinghat & Yu (2019) | SIDM + baryonic physics → SPARC rotation curve fits | ~200 |

All of these works:
1. Define a particle physics model with $(m_\chi, m_\phi, \alpha)$
2. Compute $\sigma/m(v)$ from first principles (Born, classical, or partial-wave)
3. Compare predictions against astrophysical data

Our analysis adds the relic density as an **additional** constraint that further reduces the parameter space before the astrophysical comparison — making our procedure **more constrained**, not less.

---

## Statistical Transparency Measures

To address potential concerns about the $\chi^2$ values:

1. **All 17 BPs are reported** in Table 4 of the preprint, with individual $\chi^2/\nu$ — not just the best one.
2. **Trial factor:** when quoting the minimum $\chi^2/\nu$, it should be stated that this is the best of 17 independent points.
3. **One-sided limits:** 2 of 13 observations (Bullet Cluster, Harvey+15 cluster mergers) are upper bounds that contribute $\chi^2 = 0$ when the model undershoots the limit. This effectively reduces the number of constraining data points.
4. **Large error bars:** several observations have asymmetric uncertainties spanning factors of 3–5×, making it relatively easy for a smooth $\sigma/m(v)$ curve to pass through them.
5. **MCMC posterior:** a full Bayesian analysis (§4.7) confirms that all 17 relic BPs lie within the 95% credible region — the good fits are not statistical flukes.

---

## Key Distinction: Prediction vs. Fit

| Aspect | Phenomenological fit (e.g., Ren+ 2019) | Our approach |
|--------|----------------------------------------|-------------|
| Free parameters at prediction stage | $\sigma/m$ as free function | **Zero** — $\sigma/m(v)$ from Lagrangian |
| Relic density | Not imposed | $\Omega h^2 = 0.120$ imposed **before** comparing to SIDM data |
| $\chi^2$ role | Selection criterion | Post-hoc diagnostic |
| Parameter determination | Fit to rotation curves | Bisection on relic density |

The pipeline is: **Lagrangian → relic density → $\sigma/m(v)$ → comparison**. At no point does astrophysical SIDM data feed back into the parameter selection.

---

## References

### SIDM methodology
- S. Tulin, H.-B. Yu, K.M. Zurek, "Beyond Collisionless Dark Matter: Particle Physics Dynamics for Dark Matter Halo Structure," Phys. Rev. D 87, 115007 (2013). [arXiv:1302.3898](https://arxiv.org/abs/1302.3898)
- M. Kaplinghat, S. Tulin, H.-B. Yu, "Dark Matter Halos as Particle Colliders: Unified Solution to Small-Scale Structure Puzzles from Dwarfs to Clusters," Phys. Rev. Lett. 116, 041302 (2016). [arXiv:1508.03339](https://arxiv.org/abs/1508.03339)
- A. Kamada, M. Kaplinghat, A.B. Pace, H.-B. Yu, "Self-Interacting Dark Matter Can Explain Diverse Galactic Rotation Curves," Phys. Rev. Lett. 119, 111102 (2017). [arXiv:1611.02716](https://arxiv.org/abs/1611.02716)
- T. Ren, A. Kwa, M. Kaplinghat, H.-B. Yu, "Reconciling the Diversity and Uniformity of Galactic Rotation Curves with Self-Interacting Dark Matter," Phys. Rev. X 9, 031020 (2019). [arXiv:1808.05695](https://arxiv.org/abs/1808.05695)

### Observational constraints used in SIDM cuts
- S. Tulin, H.-B. Yu, "Dark Matter Self-interactions and Small Scale Structure," Phys. Rep. 730, 1 (2018). [arXiv:1705.02358](https://arxiv.org/abs/1705.02358) — Review establishing the σ/m ~ 1 cm²/g sweet spot.
- D. Harvey, R. Massey, T. Kitching, A. Taylor, E. Tittley, "The nongravitational interactions of dark matter in colliding galaxy clusters," Science 347, 1462 (2015). [arXiv:1503.07675](https://arxiv.org/abs/1503.07675) — Cluster merger upper bound σ/m < 0.47 cm²/g.
- S.W. Randall, M. Markevitch, D. Clowe, A.H. Gonzalez, M. Bradač, "Constraints on the Self-Interaction Cross Section of Dark Matter from Numerical Simulations of the Merging Galaxy Cluster 1E 0657-56," ApJ 679, 1173 (2008). [arXiv:0704.0261](https://arxiv.org/abs/0704.0261) — Bullet Cluster bound.

### Relic density
- Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters," A&A 641, A6 (2020). [arXiv:1807.06209](https://arxiv.org/abs/1807.06209) — $\Omega_c h^2 = 0.1200 \pm 0.0012$.
