# §7.1 — CP-Violating Structure of the Coupling Space (Draft Text)

*A — draft for paper. ~550 words. Numbers from `model_validations/cp_separation/cp_separation_table.py` (BP1) and `cp_separation_MAP.py` (MAP).*

---

### §7.1 CP-Violating Structure of the Coupling Space

Our mixed Majorana Lagrangian $\frac{1}{2}\bar{\chi}(y_s + iy_p\gamma_5)\chi\,\phi$ introduces two independent Yukawa couplings, scalar ($y_s$) and pseudoscalar ($y_p$), parametrized through $\alpha_{s,p} = y_{s,p}^2/(4\pi)$. The relic density constraint $\langle\sigma v\rangle_0 = 2\pi\alpha_s\alpha_p/m_\chi^2$ fixes the product $\alpha_s \times \alpha_p$, while self-interactions depend only on $\alpha_s$ through the Yukawa potential $V(r) = -\alpha_s e^{-m_\phi r}/r$. This creates a **CP-separation band**: a continuous family of viable models parametrized by $\alpha_s/\alpha_p$, ranging from the CP-symmetric point ($\alpha_s = \alpha_p$) to highly CP-violating configurations ($\alpha_s \gg \alpha_p$ or vice versa).

To map this band, we fix the mass spectrum at each benchmark point and scan $\alpha_s$ while enforcing three constraints: (i) the relic product $\alpha_s \alpha_p = 1.387 \times 10^{-7}$, (ii) dwarf-scale SIDM: $\sigma/m(30\text{ km/s}) \in [0.1, 10]$ cm$^2$/g, and (iii) cluster safety: $\sigma/m(1000\text{ km/s}) < 1$ cm$^2$/g.

**BP1 masses** ($m_\chi = 20.69$ GeV, $m_\phi = 11.34$ MeV). Here $\lambda = \alpha_s m_\chi/m_\phi$ ranges from 2.4 to 9.9 across the viable band, straddling the first Born resonance at $\lambda = \pi$. The 13 viable points span $\alpha_s/\alpha_p \in [13, 212]$ — a dynamic range of 1.22 decades. Points with $\lambda < \pi$ exhibit an $s$-wave scattering plateau at low velocities, while those with $\lambda > \pi$ transition to the classical regime with mild velocity dependence $\sigma_T \propto (\ln\lambda)^2/v^2$.

**MAP masses** ($m_\chi = 90.64$ GeV, $m_\phi = 13.85$ MeV). The high mass ratio $m_\chi/m_\phi = 6{,}545$ places even modest $\alpha_s$ values deep in the resonant regime ($\lambda \gg \pi$). **All 500 scanned points pass viability**, spanning $\alpha_s/\alpha_p \in [1.8, 11{,}532]$ — a dynamic range of **3.81 decades**. The oscillating VPM cross-section maintains $\sigma/m(30) > 0.1$ cm$^2$/g via resonance peaks even at very small $\alpha_s$, while the $\sigma/m \propto 1/v^4$ suppression at cluster velocities keeps $\sigma/m(1000)$ between 0.003 and 0.19 cm$^2$/g across the entire band — safely below the Bullet Cluster limit.

Table~\ref{tab:cp_bp1} and Table~\ref{tab:cp_map} present the full viable bands. The key conclusion is that **CP violation is a generic prediction** of the model: the relic-SIDM constraints permit coupling asymmetries spanning 1–4 orders of magnitude depending on the mass spectrum, with the MAP region showing the widest viable band. At the MAP benchmark itself ($\alpha_s = 2.546 \times 10^{-2}$, $\alpha_s/\alpha_p = 4{,}674$), the coupling ratio exceeds $10^3$, indicating strong CP violation as the natural state of this parameter space.

The CP-separation band has two important implications: (1) **phenomenological distinguishability** — different points along the band predict identical self-interactions but different annihilation signatures, relic pathways, and (in principle) collider phenomenology in CP-sensitive observables; (2) **robustness** — the model's SIDM predictions are insensitive to the degree of CP violation, depending only on $\alpha_s$ and the mass spectrum. Any future measurement constraining $\alpha_s/\alpha_p$ (e.g., through indirect detection channels sensitive to $s$-wave vs. $p$-wave annihilation) would select a unique point within the band without affecting the astrophysical predictions.

---

*References: Kaplinghat+2016, Tulin+2013 (VPM), Randall+2008 (Bullet Cluster), Harvey+2015.*
