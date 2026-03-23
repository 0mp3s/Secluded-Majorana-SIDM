# Referee Report — Manuscript: "Self-Interacting Dark Matter from a Majorana Fermion with Light Scalar Mediator: Velocity-Dependent Cross Sections via Yukawa Partial-Wave Analysis"

**Journal:** Physical Review D  
**Referee:** GPT-5.4  
**Date:** March 22, 2026  
**Recommendation:** **Reject**

---

## Summary

This manuscript is more coherent than the author's earlier dual-mediator construction. It removes the gauge-sector complications, explicitly distinguishes raw from representative scan counts, and presents a more honest numerical error budget for the VPM cross section. Those are real improvements.

However, the present version still does not meet publication standard. The central physics case rests on a relic-density mechanism and a cosmological argument that are, in my view, fundamentally flawed. In addition, the manuscript's own numerical uncertainty estimate is large enough to undermine the claimed "island of viability," and the Higgs-portal exclusion argument is built on an incompletely specified portal sector.

The most serious problem is not stylistic. It is that the paper appears to treat a **stable, thermalized 8–15 MeV scalar** as if it were harmless dark radiation with

$$
\Delta N_{\rm eff} \approx 0.027,
$$

while simultaneously using a standard single-species WIMP freeze-out equation for a secluded sector whose temperature history is not derived. That is not a controlled cosmological framework.

I therefore recommend rejection in the present form.

---

## Major Findings

### 1. The relic-density mechanism is not convincingly established; the claim that $\chi\chi \to \phi\phi$ is s-wave for this Majorana model is highly doubtful

Section 6.1 states that for a Majorana fermion with Yukawa coupling

$$
\mathcal{L}_{\rm int} = -\frac{y}{2}\bar\chi\chi\phi,
$$

the annihilation

$$
\chi\chi \to \phi\phi
$$

has the s-wave cross section

$$
\langle \sigma v \rangle = \frac{y^4}{64\pi m_\chi^2} = \frac{\pi \alpha^2}{4 m_\chi^2}.
$$

This is the central formula on which the relic-density argument, the "island of viability," and the comparison with the earlier axial-vector model all depend.

I do not find this claim credible as written.

For two identical Majorana fermions annihilating into two identical real scalars through a CP-even Yukawa coupling, the standard selection-rule argument strongly suggests that the s-wave amplitude is absent and that the leading contribution is p-wave. The manuscript does not address this issue at all. It simply writes down a velocity-independent rate and builds the rest of the paper on it.

If the leading annihilation is actually p-wave rather than s-wave, the entire phenomenological story changes:

1. the "exact overlap" between relic-density couplings and SIDM couplings disappears,
2. the 17-point cosmological island is not established,
3. and the paper's claimed advantage over the axial-vector model largely collapses.

This is not a minor correction. It strikes at the paper's main result.

At minimum, the authors must supply a correct threshold expansion of the full $t/u$-channel amplitude for identical Majorana fermions and explicitly demonstrate the presence of a non-vanishing s-wave term. In the absence of such a derivation, I regard the relic-density section as unreliable.

---

### 2. The cosmology of a stable 8–15 MeV scalar is treated incorrectly

The manuscript adopts a fully secluded model in which:

1. $\phi$ has no SM coupling,
2. $\phi$ is stable,
3. and the dark sector is assumed to have been in thermal equilibrium with the SM at early times with $\xi \approx 1$.

Section 5.3 then argues that this stable scalar is cosmologically harmless because it contributes only

$$
\Delta N_{\rm eff} \approx 0.027
$$

and because

$$
\frac{\Omega_\phi}{\Omega_\chi} \sim \frac{m_\phi}{m_\chi} \sim 5 \times 10^{-4}.
$$

This discussion is not acceptable.

First, the text explicitly states that the stable scalar contributes as dark radiation at temperatures

$$
T \ll m_\phi.
$$

That is simply wrong. At $T \ll m_\phi$, the scalar is non-relativistic and contributes to the matter density, not to radiation.

Second, even if one temporarily ignored that conceptual mistake and treated $\phi$ as relativistic, the quoted number

$$
\Delta N_{\rm eff} \approx 0.027
$$

does not follow from the formula given in the manuscript. For a single real scalar with the same temperature as the neutrinos one would have

$$
\Delta N_{\rm eff} = \frac{4}{7}\left(\frac{T_\phi}{T_\nu}\right)^4 = \frac{4}{7} \approx 0.57,
$$

not 0.027. To obtain 0.027, one would need a substantially colder dark sector, but the paper assumes $\xi \approx 1$ and does not derive any suppression of $T_\phi / T_\nu$.

Third, the late-time abundance estimate

$$
\Omega_\phi / \Omega_\chi \sim m_\phi / m_\chi
$$

has no clear basis. A stable thermalized 10 MeV boson is not a negligible correction to the DM density. If it ever shared the thermal bath with $\xi \sim 1$, its relic abundance is generically enormous unless a depletion or decay mechanism is provided. In fact, such a species would typically overclose the Universe by many orders of magnitude.

In other words, the cosmology of the secluded stable scalar is not merely incomplete. It is qualitatively mischaracterized.

This issue alone is severe enough to invalidate the paper's claims of cosmological safety.

---

### 3. The relic-density calculation is not "exact" for the secluded model being proposed

Section 6.2 presents a single Boltzmann equation for $Y_\chi$ and describes the result as an "exact numerical Boltzmann solver." That description is much too strong for the model under discussion.

The secluded setup assumes:

1. an early epoch of SM–dark thermal equilibrium via unspecified UV physics,
2. a hidden-sector temperature ratio $\xi \lesssim 1$ but effectively $\xi \approx 1$,
3. a stable mediator $\phi$,
4. and annihilation $\chi\chi \leftrightarrow \phi\phi$ fully within the dark sector.

In such a framework, the standard one-species visible-sector WIMP equation is not, by itself, an exact treatment. One should in general track:

1. the hidden-sector temperature evolution,
2. entropy transfer within the dark sector,
3. the mediator abundance,
4. and the conditions under which $\xi \approx 1$ is actually maintained.

The paper does none of this. Instead it assumes some UV interaction strong enough to thermalize the sectors, but negligible for every other purpose. That may be possible in a complete model, but the manuscript does not present such a model, nor does it quantify the required portal strength or decoupling epoch.

Thus the relic-density calculation is not exact for the secluded framework actually under discussion. It is a simplified one-equation approximation contingent on unstated assumptions about the hidden-sector thermal history.

Since the 17-point viable island is generated entirely from this calculation, the omission is serious.

---

### 4. The paper's own numerical error budget undermines the claimed viable island

One of the better features of this version is that it finally acknowledges sizable numerical uncertainties. Unfortunately, those uncertainties are large precisely where the paper claims success.

Section 3.4 reports a total systematic uncertainty of:

1. **2–7% at 30 km/s**,
2. **27–40% at 1000 km/s**.

The primary benchmark is then quoted as:

$$
\sigma/m(30\,\mathrm{km/s}) = 0.516\;\mathrm{cm}^2/\mathrm{g},
$$

$$
\sigma/m(1000\,\mathrm{km/s}) = 0.072\;\mathrm{cm}^2/\mathrm{g}.
$$

This is not a robust benchmark.

At the low-velocity end, a 7% downward shift gives approximately

$$
0.516 \times 0.93 \approx 0.48,
$$

which falls below the manuscript's relaxed dwarf threshold of 0.5 cm$^2$/g.

At the cluster-velocity end, a 40% upward shift gives approximately

$$
0.072 \times 1.40 \approx 0.101,
$$

which lies slightly above the claimed cluster bound of 0.1 cm$^2$/g.

So the paper's showcase point is not safely viable under its own stated numerical systematics.

The problem is even worse for the 17-point island, whose quoted range extends up to

$$
\sigma/m(1000) = 0.099\;\mathrm{cm}^2/\mathrm{g}.
$$

With a 27–40% uncertainty, such points are plainly not established as cluster-safe.

The authors try to dismiss this by saying that points with $\sigma/m(1000) \ll 0.1$ remain safely viable. But that is not the region they actually highlight. Their island sits directly on top of the region where the numerical uncertainty is largest and the exclusion cut is most consequential.

This is a major internal inconsistency.

---

### 5. The Higgs-portal exclusion argument is not formulated consistently

Section 5.1 states that the "most economical" connection to the Standard Model is the portal

$$
\lambda_{H\phi} |H|^2 \phi^2,
$$

and then immediately proceeds to interpret this in terms of an induced scalar mixing angle $\sin\theta$ that governs both direct detection and the decay $\phi \to e^+e^-$.

As written, this is not consistent.

The quartic operator $|H|^2 \phi^2$ does not by itself generate Higgs–scalar mixing after electroweak symmetry breaking unless additional structure is present, such as:

1. a linear portal term,
2. a nonzero vacuum expectation value for $\phi$,
3. or explicit breaking of a $\phi \to -\phi$ symmetry.

The manuscript does not define such a portal sector. It simply introduces $\sin\theta$ as though it followed automatically from $|H|^2 \phi^2$.

That matters because the paper's argument that the Higgs portal is "generically excluded" relies entirely on comparing:

1. a direct-detection bound on $\sin\theta$,
2. and a BBN decay bound on the same $\sin\theta$.

If the portal itself is not consistently specified, the exclusion claim is not well posed.

I do not object to the general intuition that very light scalar mediators face a serious tension between direct detection and prompt decay. I object to the fact that the paper presents this as a quantitative exclusion theorem while working with an incompletely defined portal model.

---

### 6. The paper shifts its phenomenological standard midstream

The early scan is defined using the classic SIDM sweet-spot criterion:

$$
\sigma/m(30\,\mathrm{km/s}) \in [1,10]\;\mathrm{cm}^2/\mathrm{g},
$$

but the cosmological scan later relaxes this to

$$
\sigma/m(30\,\mathrm{km/s}) \ge 0.5\;\mathrm{cm}^2/\mathrm{g}.
$$

That may be defensible, but the manuscript does not manage the transition cleanly.

The abstract still foregrounds the 1–10 cm$^2$/g sweet spot, yet the primary benchmark quoted in the abstract has

$$
\sigma/m(30) = 0.52\;\mathrm{cm}^2/\mathrm{g},
$$

which is outside that original window.

This is not fatal by itself, but it reinforces the impression that the target is being moved in order to rescue the cosmological benchmarks.

If the true thesis of the paper is that 0.5 cm$^2$/g is sufficient and observationally acceptable, then that should be the central criterion from the start rather than a later relaxation.

---

### 7. The numerical claims are not auditable from the V10 folder as submitted

The V10 folder contains only:

1. the manuscript,
2. one island-of-viability plot,
3. and one benchmark velocity-profile plot.

It does **not** contain:

1. the VPM scan code,
2. the Boltzmann solver,
3. the raw 80,142-point sample,
4. the 17 benchmark-point table in machine-readable form,
5. or the scripts used to produce the figures.

For a paper whose main claims are numerical, this is a serious weakness. A referee cannot independently audit the error budget, the scan counts, the benchmark selection, or the relic-density island from the material in V10 alone.

This problem is not the main reason for rejection, but it is a real one.

---

## Additional Concerns

### 8. The statement that the secluded model is automatically free of CMB constraints is too strong

The paper argues that because the annihilation products remain in the dark sector,

$$
f_{\rm eff} = 0,
$$

and therefore the Planck bound does not apply.

That is too glib. It is true that there is no direct SM energy injection if $\phi$ never couples back to the visible sector. But once the model contains a stable thermalized mediator, the CMB is still sensitive through the **background cosmology**: dark radiation, extra matter, and hidden-sector thermal history. The paper's treatment of the CMB therefore addresses the wrong question. It removes one visible-sector constraint while leaving the dominant hidden-sector cosmology untreated.

### 9. The bound-state-formation discussion is underdeveloped and numerically suspect

Section 6.4 dismisses BSF using a Coulomb estimate

$$
\sigma v_{\rm BSF} / \langle \sigma v \rangle_{\rm ann} \sim \alpha^3 / v^2 \sim 10^{-12}
$$

at freeze-out. For the quoted viable range $\alpha \sim 5 \times 10^{-4}$–$5 \times 10^{-3}$ and typical freeze-out velocities, this estimate does not obviously evaluate to $10^{-12}$. Even if BSF ultimately remains subleading, the manuscript gives no controlled derivation and no parameter-dependent estimate across the claimed viable island.

This is a secondary issue, but the current dismissal is not persuasive.

---

## Recommendation

I recommend **rejection**.

The idea of revisiting Majorana SIDM with a light scalar mediator is worthwhile, and the present manuscript is more focused than the author's earlier attempts. But the current paper still fails at the level of foundational physics.

The main conclusions depend on a chain of assertions that are not adequately established:

1. the annihilation channel is taken to be s-wave without confronting the Majorana selection rules,
2. the stable mediator cosmology is treated incorrectly,
3. the relic-density solver is described as exact when it is not exact for the secluded hidden-sector setup,
4. and the numerical error bars are large enough to erase the claimed safety margin of the preferred benchmarks.

These are not polishing issues. They are central to the paper's viability.

---

## What Would Be Required For A Serious Resubmission

At minimum, a defensible resubmission would require all of the following:

1. **A correct threshold-level derivation of $\chi\chi \to \phi\phi$ for identical Majorana fermions**, explicitly establishing whether the leading term is s-wave or p-wave.
2. **A complete hidden-sector cosmology** for the stable mediator, including $T_{\rm dark}/T_{\rm SM}$ evolution, mediator abundance, and the late-time contribution of $\phi$ to radiation and/or matter.
3. **A relic-density calculation appropriate to the secluded model**, rather than a standard one-species WIMP equation with an assumed $\xi \approx 1$.
4. **A benchmark set robust under the stated numerical systematics**, not one sitting directly on the observational boundaries.
5. **A consistently defined Higgs-portal comparison model**, if the paper wants to claim a generic portal exclusion.
6. **A real computational supplement**, containing the scan code, raw benchmark table, and plotting scripts.

Until those points are addressed, the manuscript is not ready for publication.

---

*Referee: GPT-5.4*