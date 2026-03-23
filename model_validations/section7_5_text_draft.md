# §7.5 — Ultra-Faint Dwarfs and the Crater II Case Study (Draft Text)

*B — draft for paper. ~600 words. All numbers from `model_validations/ufd_crater/predict_ufd.py`.*

---

### §7.5 Dwarf Spheroidal Core Sizes and Ultra-Faint Dwarfs

A strong test of SIDM models is whether they produce dark matter cores of the correct size in dwarf spheroidal galaxies (dSphs), where the core-cusp problem is most acute. We apply the Kaplinghat, Tulin & Yu (2016) criterion: a halo forms an observable core at radius $r_1$ where the cumulative scattering rate satisfies

$$\rho(r_1) \cdot \frac{\sigma}{m}\bigl(v_{\rm rel}\bigr) \cdot v_{\rm rel} \cdot t_{\rm age} = 1\,,$$

with $v_{\rm rel} = \sqrt{2}\,\sigma_v$ the typical relative velocity. We adopt NFW profiles with halo masses from abundance matching (Read et al. 2017; Errani et al. 2018), concentrations from the Correa et al. (2015) relation at $z = 0$, and a conservative dynamical age $t_{\rm age} = 10$ Gyr. We compute $\sigma/m$ at each galaxy's characteristic velocity using the full partial-wave VPM calculation.

Table~\ref{tab:ufd} presents predictions for 15 galaxies: 8 classical dSphs (Fornax, Sculptor, Draco, Carina, Sextans, Leo I, Leo II, Ursa Minor) and 6 ultra-faint dwarfs (Tucana II, Segue 1, Reticulum II, Tucana III, Carina II, Grus I), plus the extended satellite Crater II as a special case study.

**Classical dSphs.** For the MAP benchmark, 6 of 8 classical dSphs produce dark matter cores ($N_{\rm scatter} > 1$), with core radii $r_{\rm core} = 462$–$985$ pc — consistent with the observed cored profiles of Fornax ($r_{\rm core} \approx 500$–$900$ pc; Walker \& Peñarrubia 2011), Sculptor, Carina, and others. The two exceptions — Sextans ($N = 0.99$) and Leo I ($N = 0.85$) — are formally below threshold but within the systematic uncertainty: adopting $t_{\rm age} = 12$ Gyr (appropriate for stellar populations older than 10 Gyr; Weisz et al. 2014) yields $N > 1$ for all 8 classicals. In contrast, BP1 and BP9 produce no cores in any dSph ($N_{\rm scatter} < 0.65$ everywhere), as their cross sections at $v \lesssim 12$ km/s are too small ($\sigma/m \sim 0.3$–$0.5$ cm$^2$/g for BP1, $\sim 0.3$ cm$^2$/g for BP9).

**Ultra-faint dwarfs.** At UFD velocities ($v \sim 2$–$4$ km/s), the MAP benchmark gives $\sigma/m \approx 1.5$ cm$^2$/g — benefiting from the resonant plateau identified in §7.2. This produces cores in 5 of 6 UFDs (Tucana II, Segue 1, Reticulum II, Carina II, Grus I), with $r_{\rm core} = 104$–$309$ pc. Tucana III ($N = 0.78$, $\sigma_v = 1.5$ km/s) is the smallest system and remains formally cuspy, though its halo mass is highly uncertain.

**Non-universality of $r_{\rm core}/r_{\rm half}$.** A key prediction of velocity-dependent SIDM is that the ratio $r_{\rm core}/r_{\rm half}$ is **not constant** across galaxies. In our model, this ratio varies from 0.19 (Crater II) to 10.6 (Segue 1) for the MAP benchmark, with a mean of $4.8 \pm 3.2$ (67\% scatter). This large scatter arises because $\sigma/m(v)$ varies significantly across the dwarf velocity range, and because $r_{\rm half}$ (a stellar tracer) need not track the dark matter core radius. This prediction distinguishes velocity-dependent SIDM from constant-$\sigma/m$ models, which generically predict $r_{\rm core}/r_{\rm half} \approx$ const (Kamada et al. 2017), and is testable with future kinematic surveys of ultra-faint satellites (e.g., Rubin Observatory LSST; Simon 2019).

### Crater II Case Study

Crater II presents an extreme test: its half-light radius ($r_{\rm half} = 1066$ pc; Torrealba et al. 2016) is among the largest of any Milky Way satellite, yet its velocity dispersion is remarkably low ($\sigma_v = 2.7$ km/s; Caldwell et al. 2017). Even the MAP benchmark, which gives the largest cross section at $v \approx 3.8$ km/s ($\sigma/m = 1.55$ cm$^2$/g), produces an SIDM core of only $r_{\rm core} = 197$ pc — just 19\% of $r_{\rm half}$. The tidal radius ($r_t \approx 3.7$ kpc) is comfortably larger, so Crater II is not tidally truncated, but its extreme stellar extent is best explained by **tidal heating**: interactions with the Milky Way potential inflate the stellar distribution without significantly altering the dark matter core (Fattahi et al. 2018; Fu et al. 2019; Sanders et al. 2018). Crater II therefore probes the joint SIDM+tidal regime rather than SIDM in isolation — a conclusion consistent with N-body simulations that require both mechanisms to reproduce its properties.

---

*References cited: Kaplinghat+2016, Read+2017, Errani+2018, Correa+2015, Walker & Peñarrubia 2011, Weisz+2014, Kamada+2017, Simon 2019, Torrealba+2016, Caldwell+2017, Fattahi+2018, Fu+2019, Sanders+2018.*
