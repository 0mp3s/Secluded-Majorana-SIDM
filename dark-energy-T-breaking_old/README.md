# Dark Energy as the T-Violating Component of Dark Matter: A Duality Conjecture

**Author**: Omer P.  
**Date**: March 28, 2026  
**Status**: Preprint / Work in Progress  
**Companion paper**: *Secluded Majorana SIDM: Constraints and Predictions* — [Zenodo DOI: 10.5281/zenodo.19225823](https://zenodo.org/records/19225823)

---

## Core Claim

We identify a universal angle in the Majorana dark matter parameter space:

$$\theta_{\text{relic}} = \arctan\!\left(\frac{1}{\sqrt{8}}\right) = 19.47°$$

This angle satisfies **both** the SIDM self-interaction constraint and the cosmological relic density simultaneously — and is independent of $m_\chi$, $m_\phi$, or $\alpha$.

The conjecture: this angle is not a coincidence. It is the point where the scalar ("electric") and pseudoscalar ("magnetic") components of the dark Yukawa coupling are in exact duality — and the T-violating pseudoscalar component generates the observed dark energy density through a Coleman-Weinberg mechanism.

---

## The Analogy

| Electromagnetism | Dark Sector |
|---|---|
| $\vec{E}$ (electric field) | $y_s = y\cos\theta$ (T-even Yukawa) |
| $\vec{B}$ (magnetic field) | $y_p = y\sin\theta$ (T-odd Yukawa) |
| Lorentz boost mixes $E \leftrightarrow B$ | $\sigma$ field rotates $y_s \leftrightarrow y_p$ |
| $\vec{B}$ does no work | $y_p$ doesn't interfere with mass term ($\gamma^5$ flips sign) |
| Moving charges repel/attract with sign flip | T-odd sector contributes **positive** vacuum energy → DE repulsion |

**Thesis**: Dark energy is not a separate field. It is the T-violating component of the dark matter interaction, made manifest through CP violation in the Majorana sector — exactly as magnetism is the relativistic component of electric fields.

---

## The Lagrangian

Starting from the secluded Majorana SIDM Lagrangian (Paper 1), we promote the CP phase to a dynamical field $\sigma$:

$$y_s + iy_p \;\longrightarrow\; y\cos\!\left(\frac{\sigma}{f}\right) + iy\sin\!\left(\frac{\sigma}{f}\right)\gamma^5 = y\,e^{i\gamma^5 \sigma/f}$$

where $\sigma$ is a dark axion with decay constant $f \sim M_{\text{Pl}}$. The SIDM + relic constraints fix:

$$\alpha_s = \alpha, \quad \alpha_p = \frac{\alpha}{8}, \quad \theta = 19.47°$$

The Coleman-Weinberg 1-loop potential of $\chi$ running in a loop generates $\sigma$-dependent vacuum energy. For the MAP benchmark ($m_\chi = 94$ MeV), with $f \approx 0.27\,M_{\text{Pl}}$:

$$V_\sigma \approx \rho_\Lambda = 2.58 \times 10^{-47} \text{ GeV}^4$$

---

## Repository Contents

### `dark_energy_exploration/`

| File | Description |
|---|---|
| `vev_scan.py` | Scans Coleman-Weinberg effective potential across 5 benchmark points and the $r = \alpha_s/\alpha$ parameter space. Corrected V2: uses hyperbola $\alpha_s \alpha_p = \alpha^2/8$ (not circle). |
| `vev_scan_results.csv` | Output: 5 BPs × 30 $r$-values × grid of ($\mu_3$, $\lambda_4$). Flags: nontrivial VEV, stability, closeness to $\rho_\Lambda$. |
| `analyze_r_minimum.py` | Fine-grained analysis around $r \approx 1.17$ minimum for BP1. Shows CW loop grows with coupling sum $\alpha_s + \alpha_p$. |
| `analyze_vev_results.py` | Statistical summary and best-parameter extraction from the CSV scan. |
| `gravitational_bridge.py` | Tests whether gravitational suppression ($m_\chi/M_{\text{Pl}}$, non-minimal coupling $\xi R \phi^2$) can bridge the 33-order gap $\Delta V \to \rho_\Lambda$. Conclusion: YES for MAP-like benchmarks with $f \sim 0.27\,M_{\text{Pl}}$. |

---

## Key Numerical Results

| Benchmark | $m_\chi$ (MeV) | $\alpha$ | $V_\sigma / \rho_\Lambda$ at $f=0.2 M_{\text{Pl}}$ | $f$ for exact match |
|---|---|---|---|---|
| BP1 | 20.7 | $1.05\times10^{-3}$ | $10^{-6.6}$ | $1.3\times10^{15}$ GeV |
| BP9 | 42.5 | $2.17\times10^{-3}$ | $10^{-4.0}$ | $2.3\times10^{16}$ GeV |
| MAP | 94.1 | $5.73\times10^{-3}$ | $10^{-1.1}$ | $6.6\times10^{17}$ GeV = $0.27\,M_{\text{Pl}}$ |
| MAP\_relic | 85.8 | $5.52\times10^{-3}$ | $10^{-1.4}$ | $5.3\times10^{17}$ GeV = $0.22\,M_{\text{Pl}}$ |

The MAP benchmark (best-fit point from MCMC in Paper 1) reproduces $\rho_\Lambda$ with a sub-Planckian decay constant $f \approx 0.27\,M_{\text{Pl}}$ — no fine-tuning.

---

## Status and Open Questions

- [x] θ=19.47° derived from first principles (SIDM + relic constraints)
- [x] Coleman-Weinberg mechanism identified
- [x] VEV scan across all 5 benchmark points
- [x] Gravitational bridge analysis
- [ ] σ thermalization in early universe (test22: ΔN_eff ~ 0 unless μ₃ > 4×10⁻¹³ GeV)
- [ ] Full CMB constraints on σ as dark axion
- [ ] Connection to gravitational wave background

**This preprint establishes priority for the θ=19.47° duality conjecture as of March 28, 2026.**

---

## Citation

If you use this work, please cite:

> Omer P. (2026). *Dark Energy as the T-Violating Component of Dark Matter: A Duality Conjecture*. Zenodo. [DOI pending]

Companion paper:
> Omer P. (2026). *Secluded Majorana SIDM: Constraints and Predictions*. Zenodo. https://zenodo.org/records/19225823
