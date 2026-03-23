# Peer Review Report
**Reviewer:** GPT-5.4
**Draft Version:** V10
**Recommendation:** **Reject / Major Revision Required**

## Scope of This Review
This review is based only on the contents of the V10 folder. In the submitted folder I found the manuscript and two figures, but no supporting scan scripts, Boltzmann code, or benchmark tables beyond what is written in the draft itself. Accordingly, I can assess the internal physics consistency of the manuscript and the plausibility of its claims, but I cannot independently reproduce the numerical statements from the V10 folder alone.

## Summary
V10 is a conceptually cleaner draft than the axial-vector versions. Moving to a real scalar mediator removes the axial-potential problem, removes gauge-anomaly concerns, and makes the self-scattering potential genuinely Yukawa-like at the operator level. Those are significant improvements.

However, the manuscript introduces a new and more serious cosmological inconsistency by making the mediator stable while simultaneously assuming a standard thermal relic calculation with an unsuppressed annihilation channel into that stable mediator. The treatment of the stable scalar as dark radiation is physically incorrect for $m_\phi \sim 8$--15 MeV, and the relic-density calculation ignores the mediator abundance and hidden-sector thermodynamics. In addition, the claimed s-wave annihilation formula for Majorana dark matter coupled to a CP-even real scalar is not derived and is very likely wrong: for identical Majorana fermions, the threshold $t$- and $u$-channel amplitudes are expected to cancel, leaving a p-wave-suppressed rate.

There is also a persistent scattering-formalism problem: the quantity labeled as the transfer cross section is not the correct transport observable for identical particles. Because the main relic-island claim depends simultaneously on the annihilation formula, the hidden-sector thermal history, and the SIDM observable, the central result of the paper is not established.

---

## 1. Fatal Issue: Stable $\phi$ Is Not Dark Radiation at $T \ll m_\phi$
**Location:** Sections 5.2--5.4 and Summary Table

**Problematic claims:**
1. The model sets the portal to zero and declares $\phi$ stable.
2. The manuscript then states that a stable $\phi$ with $m_\phi \sim 8$--15 MeV contributes as dark radiation with
$$
\Delta N_{\rm eff} = \frac{4}{7}\left(\frac{T_\phi}{T_\nu}\right)^4 \approx 0.027.
$$
3. It further claims
$$
\frac{\Omega_\phi}{\Omega_\chi} \sim \frac{m_\phi}{m_\chi} \sim 5 \times 10^{-4}.
$$

**Why this is fatal:**
This treatment is not physically correct. A particle with mass $m_\phi \sim 10$ MeV is not radiation at temperatures $T \ll m_\phi$; it is non-relativistic matter. The quoted $\Delta N_{\rm eff}$ formula applies only to a relativistic decoupled species. It cannot be used as the late-time cosmological description of a stable 10 MeV scalar.

If $\phi$ was ever a thermal species with $\xi \equiv T_{\rm dark}/T_{\rm SM} \sim 1$, then a stable MeV-scale scalar generally carries a huge relic abundance once it becomes non-relativistic. For a relativistically decoupled boson with $T_\phi \sim T_\nu$, one expects parametrically
$$
\Omega_\phi h^2 \sim \frac{2}{3}\frac{m_\phi}{93\,\mathrm{eV}}\left(\frac{T_\phi}{T_\nu}\right)^3,
$$
which for $m_\phi \sim 10$ MeV is of order $10^4$--$10^5$, not $5 \times 10^{-4}$.

So the manuscript is not facing a small dark-radiation correction. It is facing catastrophic overclosure unless an additional depletion or decay mechanism is introduced. As written, the secluded stable-mediator cosmology is not viable.

---

## 2. Fatal Issue: The Relic-Density Calculation Uses the Wrong Thermal System
**Location:** Sections 5.2, 5.3, and 6.2

**Problematic claim:**
The draft assumes that after decoupling from the Standard Model, the dark sector tracks $T_{\rm dark} = \xi T_{\rm SM}$ with $\xi \approx 1$, and then states that the standard one-species Boltzmann equation for $\chi$ applies without modification.

**Why this is fatal:**
Once the mediator is stable and the sector is secluded, the freeze-out problem is no longer the standard visible-sector WIMP problem. The abundance of $\chi$ depends on the dark-sector temperature, the mediator abundance, entropy transfer within the dark sector, and potentially coupled Boltzmann equations for both $\chi$ and $\phi$. The quantity controlling freeze-out is not simply the SM temperature with the standard WIMP equation copied over.

In particular:
1. $\phi$ is part of the thermal bath into which $\chi$ annihilates.
2. $\phi$ is stable, so its number density and entropy do not disappear.
3. The dark-sector temperature ratio $\xi$ enters both the equilibrium densities and the Hubble rate budget.

Therefore the quoted relic island is not derived from the correct cosmological system. The manuscript would need at least a coupled hidden-sector Boltzmann treatment before any claim of 17 viable cosmological benchmark points can be taken seriously.

---

## 3. Fatal Issue: The Claimed s-Wave Formula for $\chi\chi \to \phi\phi$ Is Not Established and Is Likely Incorrect
**Location:** Abstract, Introduction, and Section 6.1

**Problematic claim:**
The paper asserts that for Majorana dark matter with a CP-even Yukawa coupling,
$$
\langle \sigma v \rangle = \frac{y^4}{64\pi m_\chi^2} = \frac{\pi \alpha^2}{4 m_\chi^2},
$$
and explicitly labels this as an s-wave, velocity-independent annihilation channel.

**Why this is fatal:**
This is the crucial step that creates the claimed broad SIDM-relic overlap, but it is not derived in the paper. For identical Majorana fermions annihilating into two identical CP-even scalars through $t$- and $u$-channel exchange, the standard threshold expectation is that the $s$-wave amplitude cancels due to antisymmetrization, leaving the leading rate p-wave suppressed.

If that standard expectation is correct here, then the main physics message of V10 collapses:
1. The relic-density coupling is no longer in broad overlap with SIDM.
2. The 17-point "island of viability" is an artifact of using the wrong annihilation rate.
3. The comparison with the earlier axial-vector case is no longer valid.

At minimum, the manuscript must provide an explicit amplitude-level derivation showing why the threshold Majorana cancellation does not occur. Without that derivation, the relic-density section is not reliable.

---

## 4. Fatal Issue: The Quantity Labeled $\sigma_T$ Is Not the Correct Transport Observable for Identical Particles
**Location:** Sections 3.1--3.2

**Problematic claim:**
The manuscript states that it computes the transfer cross section, but then uses
$$
\sigma_T = \frac{2\pi}{k^2}\left[\sum_{l\,\mathrm{even}}(2l+1)\sin^2\delta_l + 3\sum_{l\,\mathrm{odd}}(2l+1)\sin^2\delta_l\right].
$$

**Why this is fatal:**
That expression is not the standard transport cross section. It is structurally a spin-weighted total cross section built out of partial-wave probabilities. The true transfer or viscosity cross section requires angular weighting to suppress forward scattering, and for identical particles the correct transport observable is usually the viscosity cross section, not the simple quantity written above. In partial waves, these transport observables involve interference between neighboring partial waves rather than only $\sin^2\delta_l$ terms.

This matters because the entire scan is filtered against the dwarf/cluster SIDM criteria using this quantity. If the wrong observable is being used, then the reported 80,142 viable points and the 17-point relic island are not physically established.

---

## 5. Major Issue: The Higgs-Portal Analysis Is Built on an Incorrect Operator Statement
**Location:** Sections 2 and 5.1

**Problematic claim:**
The manuscript identifies the Higgs portal as
$$
\lambda_{H\phi}|H|^2 \phi^2,
$$
and states that this induces a mixing angle $\sin\theta$ after electroweak symmetry breaking.

**Why this is major:**
This is not the correct general statement for a real singlet scalar. The lowest-dimension gauge-invariant portal is actually the super-renormalizable operator
$$
\mu_{H\phi}\,\phi |H|^2,
$$
and the quartic portal $|H|^2\phi^2$ by itself does not induce a mixing angle unless additional assumptions are made, such as a nonzero singlet vev or additional terms in the scalar potential. None of that is specified in the manuscript.

As written, Section 5.1 mixes together three different things:
1. a quartic portal,
2. a linear mixing angle description,
3. a decay phenomenology for $\phi$.

Those are not equivalent statements. So the conclusion that the Higgs portal is generically excluded is much less general than the manuscript claims.

---

## 6. Major Issue: The Claimed "Exact Overlap" Is Obtained Only After Relaxing the SIDM Criterion
**Location:** Abstract, Sections 4.2, 4.4, 6.2, 7, and 8.2

**Problematic claim:**
The paper repeatedly says that the relic-density solution lies squarely within the SIDM-viable region and even claims "exact overlap".

**Why this is major:**
The original scan criterion is
$$
\sigma/m(30\,\mathrm{km/s}) \in [1,10]~\mathrm{cm}^2/\mathrm{g},
$$
but the cosmological island in Section 4.4 is defined only after relaxing the dwarf-scale threshold to
$$
\sigma/m(30\,\mathrm{km/s}) \ge 0.5~\mathrm{cm}^2/\mathrm{g}.
$$

Indeed, the headline benchmark has
$$
\sigma/m(30) = 0.516~\mathrm{cm}^2/\mathrm{g},
$$
and the full 17-point island spans only 0.50--0.79 cm$^2$/g. So the paper has not shown exact overlap with its original 1--10 cm$^2$/g SIDM sweet spot. It has shown overlap only with a later, relaxed criterion.

That may still be physically interesting, but it is not what the abstract and summary currently claim.

---

## 7. Major Issue: The Cluster-Scale Benchmarks Are Not Numerically Robust
**Location:** Sections 3.4 and 4.4, Appendix B

**Problematic claim:**
The paper presents 17 viable points with
$$
\sigma/m(1000) = 0.072\text{--}0.099~\mathrm{cm}^2/\mathrm{g},
$$
while quoting a 27--40% numerical systematic uncertainty at 1000 km/s.

**Why this is major:**
None of these points sits comfortably below the cluster cut once the stated numerical uncertainty is accounted for. The worst-case point at 0.099 is obviously non-robust, but even the primary benchmark at 0.072 is only marginally safe when the quoted uncertainty is as large as 40%.

So the claimed island of viability is sitting almost entirely on the edge of the numerical error bar. That weakens the central phenomenological conclusion substantially.

---

## 8. Major Issue: The CMB and Indirect-Detection Section Is Overstated
**Location:** Section 6.3

**Problematic claim:**
The paper states that because $\chi\chi \to \phi\phi$ remains in the dark sector, CMB and indirect-detection constraints simply do not apply.

**Why this is major:**
That conclusion is too strong. In a genuinely secluded model with a stable mediator, the phenomenology is not captured solely by whether SM photons or electrons are produced promptly. One must also understand the cosmological role of the stable mediator itself, including its relic abundance, its equation of state, and whether late-time dark-sector energy density modifies expansion or structure formation.

Because the manuscript has not solved the stable-mediator cosmology correctly, Section 6.3 cannot be treated as closed.

---

## 9. Minor Issue: The V10 Folder Does Not Contain the Supporting Numerical Machinery
**Location:** Submission package, not the prose itself

The V10 folder contains the manuscript and two figures, but not the scan code, the Boltzmann solver, the benchmark tables, or the data behind the figures. That does not by itself invalidate the science, but it prevents an independent audit of the numerical statements within the submitted folder. Given how central the 17-point island is to the paper, that is a meaningful weakness in the present package.

---

## Positive Changes in V10
The following are real improvements compared with the axial-vector direction:
1. A scalar mediator is a much cleaner route to a genuinely central Yukawa potential for self-scattering.
2. The draft correctly avoids the axial-vector C-parity and anomaly problems.
3. The manuscript is more focused and the main phenomenological target is clear.
4. The use of a dedicated relic-density scan, rather than only post-hoc commentary, is a good step in principle.

These strengths make the project worth continuing, but not in its current published form.

---

## What Is Required for a Publishable Revision
1. Replace the stable-mediator cosmology with a physically consistent scenario: either a decaying mediator with a specified portal, or a full hidden-sector cosmological treatment that computes the $\phi$ abundance and temperature evolution correctly.
2. Derive the annihilation amplitude for $\chi\chi \to \phi\phi$ explicitly and settle whether the leading contribution is s-wave or p-wave.
3. Replace the currently labeled $\sigma_T$ with the correct transport observable for identical-particle SIDM.
4. Rewrite the Higgs-portal discussion using a correct scalar potential and correct statement of when mixing occurs.
5. Recompute the relic-island benchmarks after items 1--4 are fixed.

---

## Final Recommendation
**Reject**

V10 is cleaner than the axial-vector drafts, but the cosmology of a stable MeV scalar mediator is treated incorrectly, the relic-density calculation is not derived from the correct hidden-sector system, and the claimed s-wave annihilation formula is not established. Together with the transport-cross-section issue, these problems invalidate the paper's main claim of a robust SIDM-relic island of viability.