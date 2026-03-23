# Peer Review — "Self-Interacting Dark Matter from a Majorana Fermion with Light Scalar Mediator"

**Reviewer:** Claude Opus 4.6 (Anthropic)
**Date:** March 23, 2026
**Draft:** V10 — Secluded Scalar Mediator
**Scope:** Full review of preprint, codebase, data files, and cross-checks
**Recommendation:** **Major Revision Required**

---

## Executive Summary

This manuscript presents a secluded Majorana SIDM model with a light scalar mediator, computing velocity-dependent cross sections via the Variable Phase Method (VPM) and identifying 17 cosmological benchmark points forming an "island of viability." The numerical pipeline is extensive and largely well-implemented. However, the work contains several physics issues — one potentially fatal — along with code bugs and internal inconsistencies that must be addressed before the results can be considered reliable.

The **most critical issue** is the claim that the annihilation $\chi\chi \to \phi\phi$ is s-wave for a Majorana fermion with a CP-even scalar coupling, which is contested by standard QFT selection rules (§1 below). The cosmological treatment of the stable mediator also requires significant correction (§2). Beyond these physics issues, I have identified concrete code bugs, data inconsistencies, and λ-convention confusion that undermine the presentation (§§3–8).

---

## PART I: PHYSICS ISSUES (CRITICAL)

### 1. FATAL — The s-wave annihilation claim requires explicit derivation and is likely incorrect

**Location:** §6.1, Eq. for $\langle\sigma v\rangle$; `core/v27_boltzmann_relic.py` line 88; `relic_density/smart_scan.py` line 71

**The claim:**
$$\langle\sigma v\rangle = \frac{y^4}{64\pi m_\chi^2} = \frac{\pi\alpha^2}{4 m_\chi^2} \quad \text{(s-wave, velocity-independent)}$$

**The problem:**
For identical Majorana fermions annihilating into two identical real scalars through a CP-even Yukawa coupling $\mathcal{L} = -(y/2)\bar{\chi}\chi\phi$, standard selection-rule arguments strongly indicate the s-wave amplitude vanishes:

1. **Parity argument:** An $L=0$ state of two identical Majorana fermions must have $S=0$ (Pauli exclusion + antisymmetric spatial wavefunction). The quantum numbers are $J^{PC} = 0^{-+}$. Two identical real scalars in an s-wave (even $L_f$) have $J^{PC} = 0^{++}$. The transition $0^{-+} \to 0^{++}$ violates parity conservation of the scalar interaction.

2. **Explicit amplitude:** In the $v \to 0$ limit, the t/u-channel Majorana amplitudes involve $\bar{v}(p_2)u(p_1)$ which vanishes for zero relative momentum.

If the leading annihilation is p-wave ($\langle\sigma v\rangle \propto v^2$), then:
- The relic-density coupling $\alpha_{\rm relic}$ shifts upward by roughly $\sqrt{x_f/6} \sim 3$–$5\times$
- The "exact overlap" between SIDM-viable and relic-viable parameter space collapses
- All 17 benchmark points produce too much dark matter (overclose the universe)
- The entire "island of viability" is invalidated

**Note:** Both previous reviewers (GPT-5.4, Gemini 3.1 Pro) independently flagged this same issue as fatal. The fact that `cosmology/bsf_estimate.py` internally uses a p-wave formula $a_2 = 3\alpha^2/(16m_\chi^2)$ — *in the same codebase* — adds further evidence that the s-wave claim is incorrect.

**Required action:** Provide an explicit threshold expansion of the $\chi\chi \to \phi\phi$ amplitude at the diagram level, showing whether the $v^0$ term survives. If it does not, the entire relic-density analysis must be redone with the p-wave formula.

---

### 2. FATAL — The cosmology of a stable ~10 MeV scalar is fundamentally mischaracterized

**Location:** §5.3 (Cosmological Safety); `predictions/delta_neff/predict_neff.py`

**Three distinct sub-problems:**

#### 2a. The paper claims $\Delta N_{\rm eff} \approx 0.027$, but the code computes $\Delta N_{\rm eff} \approx 0$

The `predict_neff.py` script correctly identifies that for $m_\phi \sim 10$–14 MeV with $T_\phi \lesssim 0.5$ MeV at BBN, the Boltzmann factor $e^{-m_\phi/T_\phi}$ is catastrophically small ($\ll 10^{-8}$). The code's conclusion: **$\Delta N_{\rm eff} \approx 0$** (the mediator is non-relativistic and Boltzmann-suppressed at all relevant epochs).

The paper, however, quotes $\Delta N_{\rm eff} \approx 0.027$ — a massless-limit formula applied to a *massive* particle. These are contradictory. Moreover:

- If $g_{*S}(T_{\rm dec}) = 86.25$ (correct for $T \sim 20$–90 GeV, below $m_t$), the massless limit gives $\Delta N_{\rm eff} \approx 0.036$, not 0.027
- The value 0.027 apparently uses $g_{*S} = 106.75$ (above top quark threshold), which is inconsistent with the BP1 freeze-out temperature

#### 2b. A stable 10 MeV scalar is non-relativistic matter, not dark radiation

At $T \ll m_\phi$, the scalar is non-relativistic. Its energy density is $\rho_\phi \sim m_\phi n_\phi$, contributing to the **matter budget**, not radiation. Treating it as dark radiation (as §5.3 does) is conceptually wrong.

#### 2c. A thermalized stable 10 MeV boson overclosure problem

If $\phi$ was ever in thermal equilibrium with $\xi \approx 1$, a stable scalar with $m_\phi \sim 10$ MeV that becomes non-relativistic at $T \sim m_\phi$ maintains a relic number density comparable to neutrinos. The resulting energy density:

$$\Omega_\phi h^2 \sim \frac{m_\phi}{93\,\text{eV}} \left(\frac{T_\phi}{T_\nu}\right)^3 \sim \frac{10^7\,\text{eV}}{93\,\text{eV}} \times 0.5 \sim 5 \times 10^4$$

This is **catastrophic overclosure** by $\sim 4$ orders of magnitude. The claim $\Omega_\phi/\Omega_\chi \sim m_\phi/m_\chi \sim 5 \times 10^{-4}$ in §5.3 has no valid derivation — it appears to assume the $\phi$ number density tracks the $\chi$ freeze-out abundance, which is physically wrong for a bosonic species with its own thermal distribution.

**Required action:** Either (a) demonstrate that $\phi$ is efficiently depleted by $\phi\phi \to \chi\chi$ or $\phi\phi \to \phi\phi$ number-changing processes before BBN, (b) solve the coupled Boltzmann equations for both $\chi$ and $\phi$, or (c) acknowledge that an additional mechanism (e.g., a small $\phi^3$ or $\phi^4$ self-coupling) is needed to deplete the scalar abundance.

---

### 3. MAJOR — The single-species Boltzmann equation is inappropriate for the secluded model

**Location:** §6.2; `core/v27_boltzmann_relic.py`

The relic density is computed using a standard single-species WIMP Boltzmann equation:

$$\frac{dY}{dx} = -\sqrt{\frac{\pi}{45}} \frac{g_{*S}}{\sqrt{g_*}} M_{\rm Pl} m_\chi \frac{\langle\sigma v\rangle_0}{x^2}(Y^2 - Y_{\rm eq}^2)$$

This is the standard textbook treatment for visible-sector WIMPs annihilating into SM particles. In the secluded model, $\chi\chi \to \phi\phi$ annihilates into a dark-sector bath. The correct treatment requires:

1. Tracking the dark-sector temperature separately from the SM temperature
2. Including the $\phi$ abundance and its back-reaction on $\chi$ freeze-out
3. Accounting for entropy transfer within the dark sector ($\phi\phi \to$ heating the dark bath)
4. Justifying why $\xi \equiv T_{\rm dark}/T_{\rm SM} \approx 1$ is maintained through freeze-out

The manuscript hand-waves this by invoking a heavy UV mediator $\Sigma$ to thermalize the sectors, but does not parameterize the decoupling epoch, the subsequent $\xi$ evolution, or the dark-sector entropy budget.

---

### 4. MAJOR — σ_T is an elastic cross section, not the transport observable

**Location:** §3.2; `core/v22_raw_scan.py` function `sigma_T_vpm`

The formula used:
$$\sigma_T = \frac{2\pi}{k^2}\left[\sum_{l\,\text{even}}(2l+1)\sin^2\delta_l + 3\sum_{l\,\text{odd}}(2l+1)\sin^2\delta_l\right]$$

This is a spin-weighted **total elastic** cross section, not the momentum-transfer (transport) cross section $\sigma_T^{\rm tr} = \int(1-\cos\theta)\,d\sigma/d\Omega\,d\Omega$ relevant for SIDM halo phenomenology. The appendix (§C) acknowledges a ~20% difference in the Born regime, but this discrepancy can grow in the resonant regime where forward-scattering dominance becomes significant.

The preprint note in §3.2 claims agreement "to within ~20%," citing Appendix A.5. This is adequate for order-of-magnitude studies but marginal for a quantitative $\chi^2$ fit. The magnitude and direction of the systematic (overestimating vs underestimating) should be reported for the benchmark points.

---

### 5. MAJOR — The Higgs portal exclusion argument uses an incomplete operator basis

**Location:** §5.1

The portal is identified as $\lambda_{H\phi}|H|^2\phi^2$ (dimension-4 quartic). However, for a real scalar singlet, the lowest-dimension portal operator is the super-renormalizable trilinear:

$$\mu_{H\phi}\,\phi |H|^2$$

This operator induces mixing without requiring a singlet VEV. The quartic $|H|^2\phi^2$ by itself does **not** induce a mixing angle $\sin\theta$ unless additional terms (like a singlet VEV or $Z_2$-breaking) are present. The derivation in §5.1 implicitly assumes $\phi$ develops a VEV or that a trilinear is present, but this is not stated.

The exclusion conclusion is likely still valid (the $1/m_\phi^4$ enhancement is generic), but the operator analysis should be corrected for precision.

---

## PART II: CODE BUGS AND INCONSISTENCIES

### 6. BUG — NFW ρ_s has an extra factor of 1/3 (3× too low)

**Location:** `predictions/rotation_curves/predict_core_sizes.py` lines 81–82; also in `sensitivity_analysis.py` and `fit_sparc_baryons.py`

```python
delta_c = (200.0 / 3.0) * c ** 3 / f_c    # correct
rho_s = rho_crit * delta_c / 3.0            # BUG: extra /3
```

The standard NFW definition is $\rho_s = \delta_c \cdot \rho_{\rm crit}$ where $\delta_c = (200/3) c^3 / f(c)$. The extra `/3.0` makes $\rho_s$ 3× too low:

- $r_s$ is derived from $V_{\rm max}$, so it compensates ($r_s \to \sqrt{3}\,r_s$), preserving the peak velocity
- But the **density profile shape is distorted**: the inner density is $\sim 1/\sqrt{3}$ too low and the scale radius $\sim\sqrt{3}$ too large
- All predicted SIDM core sizes ($r_1$) and $V(2\,\text{kpc})$ values are affected
- This bug is replicated in 3 files (predict_core_sizes.py, sensitivity_analysis.py, fit_sparc_baryons.py)

**Impact:** Moderate — affects the rotation curve diversity predictions quantitatively but not the qualitative conclusions, since the SIDM thermalization estimate is order-of-magnitude anyway.

---

### 7. BUG — CSV column name mismatch between smart_scan.py and run_mcmc.py

**Location:** `relic_density/smart_scan.py` line 152 writes `m_phi_GeV`; `stats_mcmc/run_mcmc.py` reads `m_phi_MeV`

```python
# smart_scan.py:
writer.writerow(['m_chi_GeV', 'm_phi_GeV', 'alpha', ...])

# run_mcmc.py:
mp = float(row['m_phi_MeV'])    # KeyError if reading smart_scan output
```

The research journal acknowledges a "units bug" that was fixed, but the current code still shows the mismatch. Either the CSV was manually edited post-generation or an intermediate script performs the conversion. This is a reproducibility concern.

---

### 8. λ convention inconsistency — factor of 2 difference across codebase

**Files affected:** All prediction scripts, MCMC, and paper text

| Context | Definition | BP1 value |
|---------|-----------|-----------|
| VPM solver (`v22_raw_scan.py`) | $\lambda = \alpha m_\chi / m_\phi$ | 1.91 |
| Smart scan CSV | $\lambda = \alpha m_\chi / m_\phi$ | 1.91 |
| Paper body (Tables 1–2) | $\lambda = \alpha m_\chi / m_\phi$ | 1.91 |
| Predictions (gravothermal, rotation, clusters) | $\lambda = 2\alpha m_\chi / m_\phi$ | **3.82** |
| MCMC display + Paper §4.7 | $\lambda = 2\alpha m_\chi / m_\phi$ | **3.82** |
| Journal entry: MAP λ | — | 333 (= 2× the solver value ~167) |

This is not a computational bug — the VPM solver always receives the correct internal λ. But it is a **documentation and consistency failure** that makes the paper confusing:

- Table 2 reports $\lambda \in [0.73, 32.4]$ (factor-1 convention)
- §4.7 and Figure 6 report median $\lambda = 59$ (factor-2 convention)
- These appear contradictory without a convention statement

**Required action:** Choose one convention and apply it uniformly across all code outputs, CSV files, and the manuscript. Document the choice explicitly in §3.1.

---

## PART III: DATA AND NUMERICAL ISSUES

### 9. The "island of viability" sits on the numerical error boundary

**Location:** §3.4 Error Budget; data in `v31_true_viable_points.csv`

The stated numerical systematics are:
- **30 km/s:** 2–7% total
- **1000 km/s:** 27–40% total

BP1 values: $\sigma/m(30) = 0.516$ cm²/g, $\sigma/m(1000) = 0.072$ cm²/g.

- A 7% downward shift at 30 km/s → 0.48, below the 0.5 cm²/g threshold
- A 40% upward shift at 1000 km/s → 0.101, above the 0.1 cm²/g cut
- The full 17-point range goes up to $\sigma/m(1000) = 0.099$ — right at the boundary

The island of viability is **not robust** under the paper's own error estimates. This does not invalidate the existence of viable parameter space, but it means the precise boundaries and number of viable points are uncertain.

---

### 10. SIDM "sweet spot" was silently relaxed from [1, 10] to ≥ 0.5 cm²/g

**Location:** §4.2 vs §4.4; Abstract

The original scan criterion (§4.2) is $\sigma/m(30) \in [1, 10]$ cm²/g. The cosmological island (§4.4) uses $\sigma/m(30) \geq 0.5$ cm²/g. All 17 benchmark points have $\sigma/m(30) \in [0.50, 0.79]$ — **none reaches 1 cm²/g**.

The abstract and §8.2 claim "exact overlap" between relic and SIDM regions. This is only true under the relaxed criterion. While 0.5 cm²/g is physically motivated (Kamada+2017), the paper should clearly state that the "overlap" holds only for the relaxed threshold and that the original [1, 10] criterion is not met.

---

### 11. Best-fit at scan boundary suggests incomplete exploration

**Location:** §4.6; `data/v34_results.csv`

The unconstrained $\chi^2$ best-fit has $m_\chi = 100$ GeV — the **maximum value in the scan grid**. This suggests the true statistical optimum lies beyond the scanned range. The paper acknowledges this implicitly ("the upper bound $m_\chi \leq 100$ GeV is a practical choice") but does not quantify the potential improvement from extending the scan.

---

### 12. Observational upper limits treated as measurements

**Location:** §4.6, `stats_mcmc/run_mcmc.py` chi2 function

The Bullet Cluster and Harvey+15 data points have $\sigma/m_{\rm lo} = 0$ (i.e., they are upper limits, consistent with zero). The $\chi^2$ function treats them as Gaussian measurements with an ad-hoc lower error $\sigma^{-} = 0.5 \times \text{central}$ when theory falls below the central value. This inflates $\chi^2$ for these points when the model predicts very low $\sigma/m$.

A more principled treatment would use a one-sided likelihood for upper limits. The impact is probably small given the large uncertainties, but it should be acknowledged.

---

### 13. Benchmark point labeling is inconsistent across modules

| Module | "BP1" | $m_\chi$ | $m_\phi$ | $\alpha$ |
|--------|-------|----------|----------|----------|
| Paper (Table 1) | Primary | 20.69 GeV | 11.34 MeV | 1.048×10⁻³ |
| Chi2 best relic | Rank 1 | 20.69 GeV | **9.91 MeV** | 1.048×10⁻³ |
| Cosmology config | Benchmark | **42.919 GeV** | **4.233 MeV** | 6.172×10⁻⁴ |
| boltzmann_correction.py | "BP1"/"BP2" | **18.421 GeV** | **8.20 MeV** | 8.19×10⁻⁴ |
| Error budget | BM1 | **42.919 GeV** | **4.233 MeV** | 6.172×10⁻⁴ |
| tyz_comparison.py | "BP1 relic" | 20.69 GeV | **9.91 MeV** | 1.048×10⁻³ |

At least 4 different parameter sets are labeled "BP1" or used as the primary benchmark across different modules. The cosmology and error budget analyses were run on a **completely different benchmark point** than the one featured in the paper. While each module is internally consistent, a reader trying to trace the analysis end-to-end will encounter confusion.

---

### 14. Harvey+15 velocity inconsistency

| File | System | v_infall |
|------|--------|----------|
| `predictions/cluster_offsets/clusters_data.csv` | Harvey15 combined | **1000 km/s** |
| `stats_mcmc/config.json` | 72 cluster mergers | **1500 km/s** |
| `observations/config.json` | 72 cluster mergers | **1500 km/s** |

The same measurement is evaluated at different velocities in different analyses. The MCMC and chi2 use 1500 km/s; the cluster offset prediction uses 1000 km/s. This means $\sigma/m$ for this system differs between the observational fit and the prediction module.

---

## PART IV: MINOR ISSUES AND CODE QUALITY

### 15. Dead code in predict_gravothermal.py

Lines 93–94: The first density conversion is immediately overwritten:
```python
rho_g_cm3 = rho_s * MSUN_PER_KPC3_TO_GEV_PER_CM3 * 1.783e-24  # dead code
rho_g_cm3 = rho_s * 1.989e33 / (3.086e21)**3                     # overwrites
```

### 16. Missing resonance filter in benchmark_extractor.py

The docstring describes "Step 3: Quantum resonance filter" checking non-monotonicity of $\sigma/m$ at low velocities, but this filter is **never implemented** in the code. The values at v=5 and v=10 are computed but never used as a selection criterion.

### 17. Freeze-out velocity inconsistency

- `cross_checks/sommerfeld.py`: uses $v_{\rm fo} = 0.3c$
- `cosmology/bsf_estimate.py`: uses $v_{\rm fo} = \sqrt{6/x_{\rm fo}} \approx 0.55c$

Both are common conventions but should be harmonized, especially since they give factors of ~2 difference in $S_0$.

### 18. Hardcoded paths in cross-check scripts

`literature_crosscheck.py`, `blind_sanity.py`, `blind_large.py`, `tyz_comparison.py` all use hardcoded file paths rather than the config.json system used elsewhere. Minor maintainability issue.

---

## PART V: POSITIVE ASPECTS

Despite the issues above, the work has substantial merit:

1. **The VPM solver is correctly implemented and well-validated.** The Numba-accelerated RK4 integrator, validated against scipy (< 0.001% error), Born analytics, and literature values, is a solid numerical tool. The blind testing (v28/v29) adds further confidence.

2. **The Higgs portal exclusion argument is compelling.** The $1/m_\phi^4$ enhancement creating a $3 \times 10^4$ gap between BBN and LZ bounds is a strong, original result (modulo the operator basis issue in §5 above).

3. **The VPM cross-section formula for Majorana fermions is correct.** The spin-weighted partial-wave sum with $w_l = \{1, 3, 1, 3, ...\}$ and the $2\pi/k^2$ identical-particle prefactor are verified against the Born amplitude with Majorana symmetrization in `tyz_comparison.py`.

4. **The chi-squared framework is reasonable.** The 13-system fit spanning 4 decades in velocity is a non-trivial consistency check, even with the caveats about large observational uncertainties and the upper-limit treatment.

5. **The shift from axial-vector to scalar mediator** resolves the spin-dependent potential, gauge-anomaly, and C-parity issues present in earlier versions. This represents genuine intellectual progress.

6. **The code is generally well-organized** with separation of concerns (core solvers, analysis modules, predictions, cross-checks). The config.json system (where used) enables reproducibility.

---

## SUMMARY TABLE

| # | Finding | Severity | Section |
|---|---------|----------|---------|
| 1 | s-wave annihilation claim not established | **FATAL** | §6.1, v27, smart_scan |
| 2 | Stable φ cosmology mischaracterized | **FATAL** | §5.3, predict_neff |
| 3 | Single-species Boltzmann inappropriate | **MAJOR** | §6.2, v27 |
| 4 | σ_T is elastic, not transport | **MAJOR** | §3.2, v22 |
| 5 | Higgs portal operator basis incomplete | **MAJOR** | §5.1 |
| 6 | NFW ρ_s has extra /3 (bug) | **MODERATE** | predict_core_sizes, 3 files |
| 7 | CSV column name mismatch | **MODERATE** | smart_scan ↔ run_mcmc |
| 8 | λ convention (factor 2) inconsistency | **MODERATE** | All predictions, MCMC, paper |
| 9 | Viable island on numerical error boundary | **MODERATE** | §3.4, v31 data |
| 10 | SIDM threshold silently relaxed | **MODERATE** | §4.2 vs §4.4, abstract |
| 11 | Best-fit at scan boundary | **MINOR** | §4.6, v34 data |
| 12 | Upper limits treated as measurements | **MINOR** | §4.6, MCMC |
| 13 | Benchmark point labeling confusion | **MINOR** | Multiple modules |
| 14 | Harvey+15 velocity inconsistency | **MINOR** | Cluster offsets vs MCMC |
| 15 | Dead code in gravothermal | **TRIVIAL** | predict_gravothermal |
| 16 | Missing resonance filter | **TRIVIAL** | benchmark_extractor |
| 17 | Freeze-out velocity inconsistency | **TRIVIAL** | sommerfeld vs bsf |
| 18 | Hardcoded paths in cross-checks | **TRIVIAL** | 4 cross-check scripts |

---

## RECOMMENDATION

**Major Revision Required.** The paper cannot be accepted in its present form due to the unresolved s-wave annihilation claim (§1) and the incorrect stable-mediator cosmology (§2). These are not minor corrections — they strike at the paper's central result (the island of viability and SIDM–relic overlap).

**Priority actions:**
1. Provide an explicit QFT derivation of the $\chi\chi \to \phi\phi$ threshold amplitude. If the leading term is p-wave, redo the relic density analysis.
2. Solve the $\phi$ overclosure problem — either via coupled Boltzmann equations, a depletion mechanism, or revised thermal history.
3. Fix the NFW $\rho_s$ bug and rerun rotation curve predictions.
4. Harmonize the λ convention across all outputs.
5. Correct the $\Delta N_{\rm eff}$ discussion to match the code's correct finding of Boltzmann suppression.

If the s-wave claim holds (with explicit derivation), and the $\phi$ abundance can be tamed, the remaining issues are tractable and the paper's core contribution — velocity-dependent SIDM from a minimal secluded Majorana model via full VPM analysis — is genuinely interesting and publishable.
