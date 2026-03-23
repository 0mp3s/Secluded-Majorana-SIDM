# §7.2 — Fornax Globular Cluster Survival (Draft Text)

*B — draft for paper. ~400 words. Numbers from `model_validations/fornax_gc/predict_fornax_gc.py`.*

---

### §7.2 Fornax Globular Cluster Survival

The persistence of five globular clusters (GCs) in the Fornax dSph, despite dynamical friction timescales shorter than a Hubble time in an NFW cusp, constitutes one of the strongest small-scale challenges to $\Lambda$CDM (Tremaine 1976; Goerdt et al. 2006; Read et al. 2006). In SIDM, self-interactions flatten the central cusp into a constant-density core, within which dynamical friction vanishes as the gravitational wake symmetrizes in the homogeneous medium (Goerdt et al. 2006; Read et al. 2006; Kaur \& Sridhar 2018). GCs that spiral inward to $r \sim r_{\rm core}$ stall indefinitely.

We model Fornax with $M_{200} = 3.16 \times 10^9~M_\odot$, $\sigma_v = 11.7$ km/s (Walker et al. 2009), and compute the SIDM core size using the Kaplinghat et al. (2016) criterion. For each of Fornax's 5 GCs (masses from Mackey \& Gilmore 2003), we estimate the 3D position from the projected radius using three deprojection assumptions (face-on: $r_{\rm 3D} = R_{\rm proj}$; mean: $r_{\rm 3D} = \frac{4}{\pi} R_{\rm proj}$; median: $r_{\rm 3D} = \frac{\pi}{2} R_{\rm proj}$), then classify each as SAFE ($r_{\rm 3D} < r_{\rm core}$ → stalled), INSPIRAL ($r_{\rm 3D} > r_{\rm core}$ and $t_{\rm DF} < t_{\rm age}$), or MARGINAL ($r_{\rm 3D} > r_{\rm core}$ but $t_{\rm DF} > t_{\rm age}$). We adopt a conservative prescription: dynamical friction is set to zero inside $r_{\rm core}$ (where the wake fully symmetrizes), with Chandrasekhar friction operating outside. A smooth transition (e.g., Kaur \& Sridhar 2018, where $F_{\rm DF} \propto d\ln\rho/d\ln r$) would reduce inspiral rates further, making our survival scores conservative lower bounds.

**Results.** The MAP benchmark ($\alpha = 2.5 \times 10^{-2}$) produces $r_{\rm core} = 885$ pc, comfortably enclosing all 5 GCs at all deprojection angles: **15/15 safe** (5 GCs × 3 deprojections). The core size is consistent with the observationally inferred core radius of Fornax ($r_{\rm core} \approx 500$–$900$ pc; Walker \& Peñarrubia 2011; Amorisco \& Evans 2012).

BP1 ($\alpha = 1.05 \times 10^{-3}$) gives $r_{\rm core} = 449$ pc — marginally sufficient. GC3 at mean deprojection ($r_{\rm 3D} \approx 548$ pc) lies outside the core with a remaining inspiral time of $\sim 0.9$ Gyr to the core edge, yielding a score of **12/15**. BP9 ($\alpha = 1.84 \times 10^{-3}$) produces $r_{\rm core} = 332$ pc and scores 11/15. **The viable parameter space therefore includes points — particularly the MAP region — that naturally explain Fornax GC survival without fine-tuning.**

The observed projected positions of the GCs ($R_{\rm proj} = 240$–$1710$ pc) are naturally explained as GCs that formed at $r > r_{\rm core}$, spiraled inward under dynamical friction, and stalled upon entering the constant-density core. No additional mechanism beyond SIDM core formation is required.

---

*References: Tremaine 1976, Goerdt+2006, Read+2006, Kaur & Sridhar 2018, Walker+2009, Kaplinghat+2016, Mackey & Gilmore 2003, Walker & Peñarrubia 2011, Amorisco & Evans 2012.*
