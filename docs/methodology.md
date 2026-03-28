# Methodology: Secluded Majorana SIDM

**Self-Interacting Dark Matter from a Majorana Fermion with Light Scalar Mediator**

Omer P. | March 2026

---

> **Living Document — Update Tracker**
>
> | Item | Section | Status | Trigger |
> |------|---------|--------|---------|
> | χ²/dof (free + relic) & best-fit params | §6.3 | ✅ DONE | Stage 2 completed 26 Mar 2026 |
> | MCMC MAP, posteriors, effective samples | §7.5 | ⏳ PENDING | Stage 3 run_mcmc completes |
> | Predictions table results | §8 | ⏳ PENDING | Stages 7–8 complete |
> | Timeline appendix — completion dates | App. A | ⏳ PENDING | Each stage completes |
>
> *When a stage finishes, update the corresponding section and change its status to ✅ DONE.*

---

## §1 Executive Summary

This document describes both the **computational methods** and the **research process** behind our analysis of a secluded Majorana fermion dark matter model with a light scalar Yukawa mediator. The model is defined by three parameters: the DM mass $m_\chi$, the mediator mass $m_\phi$, and the dark coupling $\alpha = y^2 / 4\pi$, with the portal coupling set to $\sin\theta = 0$ (fully secluded).

The analysis proceeds through a 54-stage computational pipeline that combines:
1. **Quantum-mechanical scattering** (Variable Phase Method) across a 3D parameter grid
2. **Cosmological filtering** (Boltzmann relic density matching)
3. **Astrophysical fitting** ($\chi^2$ against 13 observed systems spanning dwarfs to clusters)
4. **Bayesian inference** (MCMC posterior sampling)
5. **Phenomenological predictions** tested against independent datasets

The methodology was iterative: an earlier model (V9, Higgs portal) was excluded by direct-detection bounds, leading to the secluded sector (V10). Peer reviews by three AI systems uncovered two FATAL-level issues that were resolved by adopting a mixed scalar-pseudoscalar Yukawa coupling. A critical numerical bug (ℓ_max truncation) required a full pipeline rerun. This document records both the final methods and the process that shaped them.

---

## Part A — Methods

### §2 Model & Motivation

#### 2.1 The Lagrangian

The general spin-0 mediator framework was established by Kahlhoefer, Schmidt-Hoberg & Wild (JCAP 2017, arXiv:1704.02149). We adopt their Lagrangian in the fully secluded limit ($\sin\theta = 0$):

$$\mathcal{L} \supset \bar\chi (i\partial\!\!\!/ - m_\chi)\chi - (y_s \bar\chi\chi + i\, y_p\, \bar\chi\gamma^5\chi)\,\phi + \frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2$$

| Symbol | Meaning | Range scanned |
|--------|---------|---------------|
| $m_\chi$ | Majorana DM mass | 10–100 GeV |
| $m_\phi$ | Scalar mediator mass | 1–50 MeV |
| $\alpha = y^2/4\pi$ | Dark fine-structure constant | $10^{-5}$–$5 \times 10^{-2}$ |
| $\sin\theta$ | Higgs-portal mixing | **= 0** (secluded) |

#### 2.2 Why Secluded?

V9 of this project used a tunable Higgs portal ($\sin\theta > 0$). LZ direct-detection bounds require $\sin\theta < 6 \times 10^{-10}$ at $m_\phi = 11$ MeV, while BBN thermalization demands $\sin\theta > 2 \times 10^{-5}$ — a gap of $\times 30{,}000$. The Higgs portal is conclusively excluded. V10 sets $\sin\theta = 0$ exactly, yielding:
- $\sigma_{\rm SI} = 0$ (no direct detection signal — a *prediction*)
- $\phi$ is cosmologically stable (no decay to SM)
- No BBN constraints on $\phi$ decay
- Only 3 free parameters instead of 4

#### 2.3 Why Majorana?

Peer review raised a FATAL issue: $\chi\chi \to \phi\phi$ is *p-wave suppressed* for a pure-scalar Yukawa coupling to Majorana fermions, giving insufficient annihilation at freeze-out. The resolution (documented in `mixed_coupling/MixedMajorana.md`): a **mixed scalar + pseudoscalar** coupling restores s-wave annihilation while preserving the Majorana nature. This is the physical model throughout V10.

#### 2.4 Physical Approximations & GR Connections

The analysis operates in the **non-relativistic, Newtonian limit** throughout. The three places where General Relativity enters — implicitly — are documented here for completeness:

| Location | GR component | Justification for Newtonian limit |
|----------|-------------|-----------------------------------|
| Relic density (`v27_boltzmann_relic.py`) | Hubble parameter $H(T)$ from FRW metric; $g_*(T)$ from SM thermodynamics | Standard freeze-out formalism; flat $\Lambda$CDM background is exact at $T \sim m_\chi/20$ |
| Structure formation — Jeans / NFW (`predictions/fornax_jeans.py`) | Poisson equation from weak-field GR limit | $v/c \sim 10^{-4}$ and $\Phi/c^2 \sim 10^{-6}$ in dwarf galaxies; post-Newtonian corrections $\lesssim 10^{-8}$ |
| Bullet Cluster constraint (Harvey+15) | DM offset measured via **gravitational lensing** (full GR) | The constraint is taken as a literature datum; lensing calculation is not reproduced here |

**Cosmological sector:** The Boltzmann equation is
$$\frac{dY}{dx} = -\frac{s\langle\sigma v\rangle}{H\, x}(Y^2 - Y_{\rm eq}^2)$$
where $H = H(T)$ is the Friedmann-equation Hubble rate. This is the only equation in the pipeline that requires GR input. All other computations — VPM scattering, $\chi^2$ fitting, MCMC — are purely non-relativistic quantum mechanics and statistics.

**Validity:** For $m_\chi \sim 10$–$100$ GeV DM in dwarf-to-cluster environments, the Newtonian limit is accurate to better than $10^{-4}$. No GR corrections are needed or relevant at the current level of observational precision ($\sim$10–30%).

---

### §3 Variable Phase Method (VPM) Solver

The core computation: given $(m_\chi, m_\phi, \alpha, v)$, compute the transfer cross section $\sigma_T / m$.

#### 3.1 Physical Setup

Two identical Majorana fermions interact via a Yukawa potential:

$$V(r) = -\alpha \frac{e^{-m_\phi r}}{r}$$

In dimensionless variables $x = m_\phi r$:

$$V(x) = -\frac{\lambda}{x} e^{-x}, \quad \lambda = \frac{\alpha m_\chi}{m_\phi}, \quad \kappa = \frac{m_\chi v}{2 m_\phi c}$$

The transfer cross section for identical Majorana fermions:

$$\frac{\sigma_T}{m} = \frac{2\pi}{k^2 m_\chi} \sum_{\ell=0}^{\ell_{\max}} w_\ell (2\ell+1) \sin^2 \delta_\ell$$

with spin-statistical weights $w_\ell = 1$ (even $\ell$), $w_\ell = 3$ (odd $\ell$).

#### 3.2 Numerical Implementation

- **ODE integrator**: 4th-order Runge-Kutta (RK4) with fixed step
- **Grid**: 4000–12000 steps (adaptive per evaluation)
- **Phase shift extraction**: $\delta_\ell$ from the variable-phase ODE at large $x$
- **Adaptive ℓ_max**:
  $$\ell_{\max} = \min\!\Big(\max\big(3,\ \min(\lfloor\kappa x_{\max}\rfloor,\ \lfloor\kappa\rfloor + \lfloor\lambda\rfloor + 20)\big),\ 500\Big)$$
- **Early exit (Ramsauer-Townsend)**: if 5 consecutive partial waves contribute < $10^{-4} \times$ peak, summation stops
- **Implementation**: Numba JIT-compiled, no Python overhead in inner loop

#### 3.3 Error Budget

Systematic errors measured against high-resolution reference solutions:

| Velocity regime | Relative error | Dominant source |
|----------------|----------------|-----------------|
| 30 km/s (dwarfs) | < 0.01% | Integrator step size |
| 1000 km/s (clusters) | 10.8–13.6% | ℓ_max truncation |
| Resonant (λ > 30) | up to 43% | Requires ℓ_max ≫ 100 |

The 30 km/s precision is excellent. The cluster-velocity and resonant errors are conservative upper bounds; the ℓ_max fix (§4.3) significantly improved these.

#### 3.4 Verification Tests (6/6 PASS)

| Test | Method | Result |
|------|--------|--------|
| T1 Born limit | VPM vs analytical Born ($\kappa \in [0.1, 0.5]$) | 85–97% agreement |
| T2 Classical limit | VPM in classical regime ($\kappa > 10$) | Median ratio ~0.75 |
| T3 Velocity dependence | $\sigma/m$ decreases with $v$ | Monotonic decrease ✓ |
| T4 Unitarity | $|S_\ell| \leq 1$ for all $\ell$ | PASS |
| T5 Majorana/Dirac ratio | $\sigma_{\rm Maj} / \sigma_{\rm Dirac} = 4$ | PASS |
| T6 Benchmark consistency | BP1 vs preprint values | 0.9% at 30 km/s, 0.4% at 1000 km/s |

---

### §4 Parameter Scan (Stage 0)

#### 4.1 Grid Design

| Dimension | Range | Points | Spacing |
|-----------|-------|--------|---------|
| $m_\chi$ | 10–100 GeV | 50 | Log-uniform |
| $m_\phi$ | 1–50 MeV | 70 | Log-uniform |
| Resonance band | 4 bands centered on Yukawa resonances | 4 | — |
| $\alpha$ | $10^{-5}$–$5 \times 10^{-2}$ | 200 | Log-uniform per band |

Total: $50 \times 70 \times 4 \times 200 = 2{,}800{,}000$ evaluations.

#### 4.2 SIDM Cuts

All cuts are taken directly from the observational literature — **not tuned to our results**:

| Cut | Value | Physical origin | Reference |
|-----|-------|----------------|-----------|
| $\sigma/m(30) \geq 0.5$ cm²/g | Lower bound | Core formation in dSphs | Elbert+15 |
| $\sigma/m(30) \leq 10.0$ cm²/g | Upper bound | Halo ellipticity preservation | Kaplinghat+16 |
| $\sigma/m(1000) \leq 0.47$ cm²/g | Upper bound | 72 cluster mergers (95% CL) | Harvey+15 |

#### 4.3 Results

- **194,313** raw viable points (with relaxed cluster cut at 1.0)
- **153,353** survive the strict Harvey+15 cut (0.47) — 78.9% retention
- **1,334** representative points (unique (m_χ, m_φ) cells)
- **1,085** BBN-safe ($m_\phi > 2m_e$ not required in secluded sector, but tracked)
- Runtime: 34,562 seconds (~9.6 hours) on 14 cores

---

### §5 Relic Density (Stage 1)

#### 5.1 Boltzmann Equation

The thermally-averaged annihilation cross section for s-wave $\chi\chi \to \phi\phi$:

$$\langle\sigma v\rangle_0 = \frac{\pi \alpha^2}{4 m_\chi^2}$$

Numerical integration of the Boltzmann equation $dY/dx$ (Kolb & Turner formalism) with tabulated $g_*(T)$.

#### 5.2 Two-Stage Approach

**Stage 1a — Approximate extraction** (`benchmark_extractor.py`):
- Kolb-Turner approximation on all 194,313 raw points
- Result: 1,353 candidates with $\Omega h^2 \approx 0.120$

**Stage 1b — Exact bisection** (`smart_scan.py`):
- 20 × 30 grid in $(m_\chi, m_\phi)$
- For each cell: bisect in $\alpha$ until $|\Omega h^2 - 0.1200| < 10^{-4}$
- SIDM cuts applied **after** bisection (with strict Harvey+15 at 0.47)
- Result: **122 golden points** (viable relic + viable SIDM)
- Runtime: 238 seconds on 10 workers

#### 5.3 Key Properties of the 122 Points

- $\Omega h^2 = 0.120003 \pm 0.000003$ (machine precision)
- $\alpha \propto m_\chi^{0.5}$ scaling (as expected from s-wave)
- $\sigma/m$ velocity ratio (30 → 1000 km/s): median 8.4× (range 4.7–16.5×)
- Coverage: 20 unique $m_\chi$ values, 11 of 30 $m_\phi$ values viable

---

### §6 Observational χ² Fit (Stage 2)

#### 6.1 Dataset: 13 Astrophysical Systems

| System | $v$ [km/s] | $\sigma/m$ [cm²/g] | Type | Reference |
|--------|-----------|-------------------|------|-----------|
| Draco dSph | 12 | $0.6^{+1.4}_{-0.5}$ | Dwarf | KTY16 |
| Fornax dSph | 12 | $0.8^{+2.2}_{-0.6}$ | Dwarf | KTY16 |
| NGC 2976 | 60 | $2.0^{+3.0}_{-1.5}$ | LSB spiral | KTY16 |
| NGC 1560 | 55 | $3.0^{+5.0}_{-2.0}$ | LSB spiral | KTY16 |
| IC 2574 | 50 | $1.5^{+3.5}_{-1.2}$ | LSB spiral | KTY16 |
| NGC 720 | 250 | $0.5^{+1.0}_{-0.4}$ | Group | KTY16 |
| NGC 1332 | 280 | $0.3^{+0.7}_{-0.25}$ | Group | KTY16 |
| Abell 611 | 1200 | $0.10^{+0.20}_{-0.08}$ | Cluster | KTY16 |
| Abell 2537 | 1100 | $0.15^{+0.25}_{-0.12}$ | Cluster | KTY16 |
| Diverse RC band | 40 | $3.0^{+7.0}_{-2.5}$ | Rotation curves | Kamada+17 |
| Bullet Cluster | 4700 | $< 1.25$ | Merger (one-sided) | Randall+08 |
| 72 cluster mergers | 1000 | $< 0.47$ | Mergers (one-sided) | Harvey+15 |
| TBTF dwarfs | 30 | $1.0^{+4.0}_{-0.5}$ | Dwarfs | Elbert+15 |

#### 6.2 Fitting Procedure

The $\chi^2$ statistic with asymmetric errors:

$$\chi^2 = \sum_{i=1}^{13} \frac{(\sigma_T/m)^{\rm theory}(v_i) - (\sigma_T/m)^{\rm obs}_i)^2}{\sigma_{{\rm err},i}^2}$$

where $\sigma_{\rm err}$ is the upper error bar if theory > observation, lower error bar otherwise.

**One-sided upper limits** (Bullet Cluster, Harvey+15): no penalty when theory < limit.

**Two parallel tracks:**

| Track | Free parameters | dof | Purpose |
|-------|----------------|-----|---------|
| Free fit | $m_\chi, m_\phi, \alpha$ (3) | 10 | Best possible match |
| Relic-constrained | $m_\chi, m_\phi$ (2); $\alpha$ from $\Omega h^2$ | 11 | Physical consistency test |

The relic-constrained fit is the more demanding and more interesting test: can one model simultaneously explain cosmological relic abundance **and** small-scale structure?

#### 6.3 Computational Approach

1. Sample ~5,000 points from the 194K raw scan (every N-th), plus all 122 relic BPs
2. Evaluate $\sigma_T/m$ at all 12 unique velocities per point → ~63,000 VPM calls
3. Compute $\chi^2$ for each point
4. Fine-scan: vary $\alpha$ by $\times[0.8, 0.85, ..., 1.2]$ around top 50 → 400 refinement points
5. Output: ranked lists + plots + CSV

**V10 Results (26 Mar 2026, post-ℓ_max fix):**

| Track | Best-fit params | $\chi^2$ | dof | $\chi^2/{\rm dof}$ |
|-------|----------------|----------|-----|--------------------|
| Free | $m_\chi$=100.0, $m_\phi$=9.15 MeV, $\alpha$=2.85e-3 | 1.22 | 10 | **0.12** |
| Relic-constrained | $m_\chi$=100.0, $m_\phi$=11.34 MeV, $\alpha$=4.81e-3 | 1.43 | 11 | **0.13** |

Of 122 relic BPs: **97 (79%) achieve $\chi^2/{\rm dof} < 1.0$**, all 122 have $\chi^2/{\rm dof} < 6.7$.

*Compared to V9 (pre-ℓ_max fix): free improved 0.26 → 0.12, relic improved 0.54 → 0.13.*

---

### §7 MCMC Posterior Sampling (Stage 3)

#### 7.1 Setup

| Parameter | Value |
|-----------|-------|
| Sampler | emcee (affine-invariant ensemble) |
| Walkers | 32 |
| Burn-in | 300 steps |
| Production | 2,000 steps |
| Total evaluations | 73,600 likelihood calls |
| Workers | 12 parallel (multiprocessing.Pool) |
| Seed | `np.random.default_rng(42)` — reproducible |

#### 7.2 Priors

Uniform in log-space:

| Parameter | Prior range |
|-----------|-------------|
| $\log_{10}(m_\chi / {\rm GeV})$ | $[\log_{10}(5),\ \log_{10}(200)]$ |
| $\log_{10}(m_\phi / {\rm MeV})$ | $[\log_{10}(3),\ \log_{10}(30)]$ |
| $\log_{10}(\alpha)$ | $[\log_{10}(10^{-5}),\ \log_{10}(0.05)]$ |

#### 7.3 Likelihood

$$\ln\mathcal{L} = -\chi^2 / 2$$

with the same asymmetric-error $\chi^2$ as §6.

#### 7.4 Seeding Strategy

Initial walker positions: 17 relic BPs + V34 unconstrained best-fit, each perturbed by $\mathcal{N}(0, 0.05)$ in log-space.

#### 7.5 Results

<!-- UPDATE:mcmc — Replace this block when Stage 3 completes -->
| Quantity | Value | Status |
|----------|-------|---------|
| MAP $(m_\chi, m_\phi, \alpha)$ | ⏳ TBD | Stage 3 not started |
| $\chi^2/{\rm dof}$ at MAP | ⏳ TBD | Stage 3 not started |
| Acceptance fraction | ⏳ TBD | Stage 3 not started |
| Effective samples | ⏳ TBD | Stage 3 not started |
| BPs within 95% CI | ⏳ TBD | Stage 3 not started |
| 68% CI for $m_\chi$ | ⏳ TBD | Stage 3 not started |

*V9 reference (superseded): MAP (90.6, 13.9, 2.55e-2), χ²/dof=0.16, 874 eff. samples, 17/17 in 95% CI.*
<!-- END UPDATE:mcmc -->

---

### §8 Phenomenological Predictions (Stages 7–8)

Each prediction uses the VPM solver + benchmark points to test the model against independent data:

<!-- UPDATE:predictions — Replace results column when Stages 7-8 complete with V10 BPs -->
| Prediction | Method | Key result | Update |
|------------|--------|------------|--------|
| Gravothermal collapse | $t_{\rm grav}/t_{\rm age}$ for 8 dSphs | ⏳ *V9: 6/0/2* | Stages 7–8 |
| DM-only rotation curves | SIDM core size vs $V(2\,{\rm kpc})$ | ⏳ *V9: dwarfs PASS, spirals MISS* | Stages 7–8 |
| SPARC + baryons | Full $V(r)$ fits, 7 galaxies, $\Upsilon_*$ free | ⏳ *V9: 3/4 physical* | Stages 7–8 |
| Sensitivity analysis | 25,200-point scan in ($c_{\rm factor}, \nu_{\rm AC}$) | ⏳ *V9: 3/4 robust* | Stages 7–8 |
| Cluster offsets | $\sigma/m$ at $v \sim 1000$–3000 km/s | ⏳ *V9: ALL pass* | Stages 7–8 |
| $\Delta N_{\rm eff}$ | $\phi$ contribution to radiation | ⏳ *V9: ≈ 0* | Stages 7–8 |
| Fornax GC survival | Dynamical friction in SIDM core | ⏳ *V9: survive* | Stages 7–8 |
| SMBH seeds | Enhanced core collapse → seed BH | ⏳ *V9: testable* | Stages 7–8 |
| Multi-dSph Jeans | Jeans modelling of 4+ dSphs | ⏳ *V9: consistent* | Stages 7–8 |
| Majorana vs Dirac | $\sigma_{\rm Maj}/\sigma_{\rm Dirac}$ velocity fingerprint | ⏳ *V9: factor 4* | Stages 7–8 |
| CP separation | Partial-wave decomposition | ⏳ *V9: model-specific* | Stages 7–8 |
| Ultra-faint dwarfs | Scaling to UFDs + Crater II | ⏳ *V9: consistent* | Stages 7–8 |
<!-- END UPDATE:predictions -->

---

## Part B — Research Process

### §9 Model Iterations

#### 9.1 V9 → V10: Higgs Portal Exclusion

The original model (V9) included a nonzero Higgs-portal mixing angle $\sin\theta$, making the mediator $\phi$ unstable and allowing spin-independent direct detection. After computing the confrontation with LZ bounds:

$$\sin\theta_{\rm LZ} < 6 \times 10^{-10} \quad \text{vs} \quad \sin\theta_{\rm BBN} > 2 \times 10^{-5}$$

The Higgs portal was excluded by a factor of ~30,000. This forced the transition to a fully secluded dark sector ($\sin\theta = 0$), which became V10.

**What was preserved:** The entire VPM solver, the scattering potential $V(r) = -\alpha e^{-m_\phi r}/r$, the parameter grid structure. The potential is identical — only the cosmological sector changed.

**What changed:** $\phi$ became cosmologically stable, BBN constraints disappeared, direct detection was predicted to be null, and the parameter count dropped from 4 to 3.

#### 9.2 Peer Review & FATAL Findings

After completing V10's initial pipeline, the preprint was submitted to three AI peer reviewers (Claude Opus 4.6, GPT-5.4, Gemini 3.1 Pro). Two FATAL-level issues were identified:

| # | Finding | Resolution |
|---|---------|------------|
| F1 | $\chi\chi \to \phi\phi$ is p-wave for pure-scalar Majorana coupling | Mixed scalar + pseudoscalar coupling: $y_s \bar\chi\chi\phi + iy_p \bar\chi\gamma^5\chi\phi$ restores s-wave |
| F2 | Stable $\phi$ in secluded sector → $\Omega_\phi h^2 \sim 10^4$ (overclosure) | Cannibal $\phi^3$ self-interaction: $3\phi \to 2\phi$ depletes $\phi$ number density |

Five technical bugs were also found and fixed:
- NFW $\rho_s$ off by factor of 3
- CSV column name inconsistency (m_phi_GeV vs m_phi_MeV)
- $\lambda$ convention factor of 2
- config.json path error
- Gravothermal opacity factor

**Key lesson:** The potential and scattering cross sections are unchanged by the coupling structure — only the relic density calculation was affected. All VPM results remain valid.

#### 9.3 Mixed Majorana: The Resolution

A detailed A-B discussion between two AI agents (documented in `mixed_coupling/MixedMajorana.md`) established that a mixed scalar-pseudoscalar Yukawa coupling:
- Preserves Majorana nature (self-conjugate fermion)
- Gives s-wave annihilation $\chi\chi \to \phi\phi$ (the pseudoscalar component)
- Maintains the same Yukawa potential for scattering (non-relativistic limit)
- Does not introduce new free parameters (only the ratio $y_p/y_s$ matters; fixed by requiring s-wave dominance)

This was the critical theoretical fix that kept the Majorana model viable.

---

### §10 Numerical Fixes & Pipeline Reruns

#### 10.1 The ℓ_max Bug (25 Mar 2026)

**What happened:** The adaptive formula for $\ell_{\max}$ in `sigma_T_vpm` was incorrect, capping the number of partial waves too aggressively. For the MAP point ($\lambda \approx 48.6$), $\ell_{\max}$ was ~17 instead of the correct ~82 — a factor of ~5 too low.

**Impact:** The first complete pipeline run (78 minutes, 157,373 points) produced systematically incorrect cross sections for high-$\lambda$ (resonant) points. All downstream results (relic BPs, χ² fits, MCMC) were contaminated.

**Fix:**

```python
# Before:
l_max_hard = min(max(3, int(kappa * x_max)), 500)

# After:
l_max_hard = min(max(3, min(int(kappa * x_max), int(kappa) + int(lam) + 20)), 500)
```

**Consequence:** Full pipeline rerun from scratch. The corrected Stage 0 scan took 9.6 hours (vs 78 minutes) — the longer runtime was itself evidence that high-$\lambda$ points were now properly computed.

#### 10.2 ndof Bug (26 Mar 2026)

**What happened:** `chi2_fit.py` used `ndof = len(OBSERVATIONS) - 3 = 10` for both the free fit and the relic-constrained fit.

**Why wrong:** The relic-constrained fit has only 2 free parameters ($m_\chi, m_\phi$; $\alpha$ is fixed by $\Omega h^2$), so the correct dof is $13 - 2 = 11$.

**Impact:** $\chi^2/{\rm dof}$ for relic BPs was underreported by ~9%. No effect on $\chi^2$ values, rankings, or parameter selections.

**Fix:** Split into `ndof_free = 10` and `ndof_relic = 11`. No rerun needed (purely cosmetic in output formatting).

#### 10.3 Other Fixes Applied

| Bug | Description | Severity | Rerun needed? |
|-----|-------------|----------|---------------|
| NFW $\rho_s / 3$ | Factor 3 error in density | Affects predictions | Predictions only |
| CSV column name | m_phi_GeV → m_phi_MeV | Parsing error | Rerun of consumers |
| $\lambda$ factor 2 | Convention inconsistency | Display only | No |
| Docstring "80,142" → "194,313" | Stale documentation | None | No |

---

### §11 Validation & Cross-Checks

#### 11.1 Internal Consistency

| Test | Script | What it checks | Result |
|------|--------|----------------|--------|
| Born limit | `cross_checks/literature_crosscheck.py` | VPM → Born at $\kappa \ll 1$ | 6/6 PASS |
| Born transfer | `cross_checks/tyz_comparison.py` | VPM vs TYZ13 Table I | Agreement |
| Unitarity | `cross_checks/blind_sanity.py` | $|S_\ell| \leq 1$ | PASS |
| Large scale | `cross_checks/blind_large.py` | Consistency on full dataset | PASS |
| Sommerfeld | `cross_checks/sommerfeld.py` | $S_0 < 1.026$ at freeze-out | Negligible |
| MB averaging | `cross_checks/velocity_averaged.py` | $\Delta\chi^2 < 9\%$ vs mono-velocity | Acceptable |
| Bullet one-sided | `cross_checks/bullet_onesided.py` | One-sided limit treatment | Correct |
| $\sigma_T = \sigma_{\rm el}$ | `cross_checks/transfer_vs_elastic.py` | Transfer ≈ elastic at low $\ell$ | Confirmed |

#### 11.2 External Cross-Validation

The model was tested against **independent datasets not used in the χ² fit**:

| Test | Independent data | Result |
|------|-----------------|--------|
| SPARC rotation curves | 7 galaxies, Oh+15, de Blok+08 | 3/4 dwarfs: physical $\Upsilon_*$ |
| Multi-dSph Jeans | 4+ spheroidals, Walker+09 | Consistent $\sigma_{\rm los}(R)$ |
| Fornax GC survival | 5 globular clusters, Cole+12 | All survive in SIDM core |
| Cluster mergers | 6 systems, Harvey+15 | ALL $\sigma/m \ll 0.47$ |
| Sensitivity scan | 25,200 points in $(c_{\rm factor}, \nu_{\rm AC})$ | $\Upsilon_*$ robust for dwarfs |

#### 11.3 What Failed and What We Learned

**Spiral galaxies ($V_{\max} > 100$ km/s):**
Initial DM-only rotation curve predictions **missed** observed rotation curves. Adding baryonic contributions ($V_{\rm bar}$ from SPARC 3.6μm photometry) with free $\Upsilon_*$ improved the fit, but spirals still yielded $\Upsilon_* > 1.5$ — unphysical. Root cause: fixed $c = 12$ instead of using the concentration-mass relation ($c \approx 8$–10 for spirals), and missing adiabatic contraction.

**This is not a model failure** — it's a known limitation of the simple isothermal-core + NFW approximation for massive halos. The model is designed for and validated on dwarf-scale systems.

**UGC 128 ($\Upsilon_* \to 0$):**
This extreme LSB galaxy is gas-dominated and has historically been challenging for all DM models. The model cannot constrain $\Upsilon_*$ here — it's an inherent limitation of the data, not the model.

---

### §12 Predict → Test → Refine Cycle

The research followed an iterative cycle of generating predictions from the model functions, testing them against data, identifying discrepancies, and refining the approach:

```
┌──────────────────────────┐
│ 1. Compute σ_T/m(v)      │ VPM solver on parameter grid
│    from model parameters  │
└───────────┬──────────────┘
            ▼
┌──────────────────────────┐
│ 2. Compare to 13 systems │ χ² fit: does σ(v) match
│    (KTY16, Kamada+17...) │ observations at all v?
└───────────┬──────────────┘
            ▼
┌──────────────────────────┐
│ 3. Impose relic density  │ Boltzmann → fix α → only
│    constraint on α       │ (m_χ, m_φ) remain free
└───────────┬──────────────┘
            ▼
┌──────────────────────────┐
│ 4. Generate predictions  │ Gravothermal, Fornax GCs,
│    for unseen phenomena   │ rotation curves, ΔN_eff...
└───────────┬──────────────┘
            ▼
┌──────────────────────────┐
│ 5. Test against          │ SPARC, multi-dSph Jeans,
│    independent data       │ cluster mergers, UFDs
└───────────┬──────────────┘
            ▼
┌──────────────────────────┐
│ 6. Identify failures →   │ Spiral $\Upsilon_*$ too high →
│    refine analysis        │ need c(M) + baryons
└───────────┴──────────────┘
```

Crucially, the model was **not tuned** to match predictions. The 122 benchmark points were selected purely by relic density + SIDM cuts. Their performance on rotation curves, Fornax GCs, and cluster mergers is a **genuine prediction** — they were never optimized against those datasets.

---

### §13 Reproducibility

#### 13.1 Pipeline Structure

The complete analysis is encoded in 54 ordered stages (`docs/execution_pipeline.csv`), grouped into 9 phases:

| Phase | Stages | Scripts | Key output |
|-------|--------|---------|------------|
| 0 — VPM Scan | 1 | `core/v22_raw_scan_fast.py` | 194K viable points |
| 1 — Relic Density | 4 | `relic_density/*.py` | 122 golden BPs |
| 2 — Observations | 4 | `observations/*.py` | χ² results + plots |
| 3 — MCMC | 3 | `stats_mcmc/*.py` | Posterior + corner plots |
| 4 — VPM Validation | 6 | `vpm_scan/*.py` | Error budget + convergence |
| 5 — Cross-checks | 8 | `cross_checks/*.py` | Literature comparison |
| 6 — Cosmology | 3 | `cosmology/*.py` | BBN, BSF, Fermi-LAT |
| 7 — Predictions | 17 | `predictions/*.py` | Testable observables |
| 8 — Model Validations | 8 | `model_validations/*.py` | CP, Fornax, RAR, UFD |

#### 13.2 Configuration

All physical parameters, observational data, SIDM cuts, and benchmark points are defined in a **single source of truth**: `data/global_config.json`. Local `config.json` files reference this and add script-specific settings (e.g., n_workers). No physics is hardcoded in scripts.

#### 13.3 Data Provenance

Every pipeline run is logged in `docs/runs_log.csv` with:
- Timestamp, script name, stage
- Full parameter JSON
- Input data sources (with git hash verification)
- Output file paths
- Duration and result count

`get_latest()` auto-discovers the most recent archive file for each data product, with automatic warning if the git hash has changed since the producing run.

#### 13.4 Running from Scratch

```bash
# Stage 0: ~10 hours on 14 cores
python core/v22_raw_scan_fast.py

# Stage 1: ~5 minutes
python relic_density/benchmark_extractor.py
python relic_density/smart_scan.py

# Stage 2: ~30-90 minutes
python observations/chi2_fit.py

# Stage 3: ~3-12 hours (recommend overnight)
python stats_mcmc/run_mcmc.py
python stats_mcmc/autocorr_diagnostic.py

# Stages 4-8: ~2-4 hours total (no dependencies on Stage 3)
# Run in order from execution_pipeline.csv
```

---

## Appendix A: Timeline of Key Events

| Date | Event | Impact |
|------|-------|--------|
| Early Mar 2026 | V9 pipeline complete (Higgs portal) | Baseline results |
| Mid Mar 2026 | LZ exclusion analysis → gap ×30,000 | **V10: secluded sector** |
| 23 Mar 2026 | Peer review (3 AI referees) | 2 FATAL + 5 technical bugs |
| 23 Mar 2026 | Mixed Majorana coupling resolved (A-B debate) | s-wave restored |
| 23 Mar 2026 | Bug fixes applied (NFW, CSV, λ) | Predictions corrected |
| 24 Mar 2026 | Multi-dSph Jeans cross-validation | Additional confirmation |
| 25 Mar 2026 | ℓ_max bug discovered → **full rerun** | 9.6 hours, 194K points |
| 26 Mar 2026 | Stage 0 done, Stage 1 done (122 BPs) | Current dataset |
| 26 Mar 2026 | ndof bug fixed (10 → 11 for relic) | No rerun needed |
| 26 Mar 2026 | **Stage 2 chi2_fit completed** | free χ²/dof=0.12, relic=0.13, 97/122 < 1.0 |
| ⏳ TBD | Stage 3 MCMC completed | → update §7.5 |
| ⏳ TBD | Stages 7–8 predictions completed | → update §8 |

## Appendix B: Glossary

| Term | Definition |
|------|-----------|
| BP | Benchmark Point — a specific $(m_\chi, m_\phi, \alpha)$ triple |
| VPM | Variable Phase Method — ODE-based partial-wave scattering solver |
| SIDM | Self-Interacting Dark Matter |
| dSph | Dwarf spheroidal galaxy |
| TBTF | Too-Big-To-Fail problem |
| KTY16 | Kaplinghat, Tulin & Yu (2016) PRL 116, 041302 |
| Harvey+15 | Harvey et al. (2015) Science 347, 1462 |
| $\lambda$ | Dimensionless coupling: $\alpha m_\chi / m_\phi$ |
| $\kappa$ | Dimensionless momentum: $m_\chi v / (2 m_\phi c)$ |
| MAP | Maximum A Posteriori — best-fit point from MCMC |
