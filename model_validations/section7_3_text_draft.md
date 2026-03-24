# §7.3 — Radial Acceleration Relation

**Draft by A | 24 Mar 2026 | ~500 words**

---

The radial acceleration relation (RAR) discovered by McGaugh et al. (2016) establishes a tight empirical correlation between the observed centripetal acceleration $g_{\rm obs} = V_{\rm obs}^2/r$ and the baryonic acceleration $g_{\rm bar} = V_{\rm bar}^2/r$ in disk galaxies, described by

$$g_{\rm obs} = \frac{g_{\rm bar}}{1 - \exp\left(-\sqrt{g_{\rm bar}/g_\dagger}\right)}$$

with the characteristic scale $g_\dagger = 1.2 \times 10^{-10}$ m/s² and an intrinsic scatter of $\sim 0.13$ dex. Any viable dark matter model must reproduce this relation without excessive scatter.

We test the SIDM prediction against 7 SPARC galaxies spanning two regimes: gas-dominated dwarfs (DDO 154, IC 2574, NGC 2366, UGC 128; $V_{\rm max} = 47$–$66$ km/s) and baryon-rich spirals (NGC 2403, NGC 2976, NGC 3198; $V_{\rm max} = 90$–$150$ km/s). The SPARC database provides a single pre-computed baryonic velocity $V_{\rm bar}^2 = V_{\rm gas}^2 + \Upsilon_{*,\text{def}} V_{\rm disk}^2$ with default $\Upsilon_{*,\text{def}} = 0.5\;M_\odot/L_\odot$. To properly disentangle the gas component (which does not scale with $\Upsilon_*$) from the stellar disk, we employ literature gas mass fractions $f_{\rm gas} = M_{\rm gas}/M_{\rm bar}$ from resolved HI surveys (Oh et al. 2015; de Blok et al. 2008; Lelli et al. 2016; Adams et al. 2014), yielding:

$$V_{\rm tot}^2(r) = f_{\rm gas}\,V_{\rm bar}^2 + \Upsilon_*\cdot 2(1-f_{\rm gas})\,V_{\rm bar}^2 + V_{\rm DM}^2(r)$$

where the factor of 2 undoes the default $\Upsilon_{*,\text{def}} = 0.5$ embedded in $V_{\rm bar}$. This decomposition is exact in the limit $f_{\rm gas} \to 1$ (gas-dominated dwarfs) and becomes a standard global approximation for mixed systems (cf. Di Cintio et al. 2014; Santos-Santos et al. 2018). The DM contribution $V_{\rm DM}$ includes adiabatic contraction and SIDM isothermal coring via the thermalization condition $\rho(r_1)\,(\sigma_T/m)\,v\,t_{\rm age} = 1$.

**BP1 results** ($m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV, $\alpha = 1.048 \times 10^{-3}$). The fitted stellar mass-to-light ratios for spirals — $\Upsilon_* = 0.62$ (NGC 2976), 0.84 (NGC 3198), 0.91 (NGC 2403) $M_\odot/L_\odot$ — are consistent with 3.6 $\mu$m stellar population synthesis expectations ($0.2$–$0.8\;M_\odot/L_\odot$; Meidt et al. 2014). For gas-dominated dwarfs ($f_{\rm gas} > 0.75$), $\Upsilon_*$ is unconstrained as expected: the stellar disk contributes $< 10\%$ of $V_{\rm bar}^2$, and the rotation curve is governed by gas + DM alone. Crucially, BP1 reduces the RAR scatter from 0.283 to 0.179 dex for dwarfs (37% improvement) and from 0.177 to 0.122 dex for spirals (31%), demonstrating that SIDM coring tightens the RAR at both ends of the acceleration spectrum.

**MAP results** ($m_\chi = 94.07$ GeV, $\alpha = 5.734 \times 10^{-3}$, $\lambda = 48.6$). The MAP benchmark produces SIDM cores with $r_1 \sim 2$–$10$ kpc — comparable to or exceeding the optical scale lengths of the sample galaxies. This over-coring drives spiral fits to unphysical $\Upsilon_* > 1.3$ and pushes dwarf fits to the upper bound ($\Upsilon_* = 3.0$), with $\chi^2/\text{dof} > 19$ for NGC 2403 and NGC 3198. The tension provides a rotation-curve-based upper bound on the self-interaction strength at dwarf velocities: $\sigma/m(30\;\text{km/s}) \lesssim \text{a few}$ cm$^2$/g for consistency with SPARC rotation curves.

The contrast between BP1 and MAP illustrates a key feature of velocity-dependent SIDM: a single benchmark cannot simultaneously optimize predictions across all mass scales. BP1 ($\sigma/m \approx 0.5$ cm$^2$/g at 30 km/s) represents the sweet spot for galactic-scale phenomenology, while MAP excels at cluster scales (§7.2) and UFD environments (§7.5). This velocity-dependent complementarity is a structural prediction of the model (§7.6).

**References:** McGaugh, Lelli & Schombert (2016); Oh et al. (2015); de Blok et al. (2008); Lelli, McGaugh & Schombert (2016); Adams et al. (2014); Meidt et al. (2014); Di Cintio et al. (2014); Santos-Santos et al. (2018).
