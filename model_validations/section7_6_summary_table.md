# §7.6 — Summary Table (Draft Skeleton)

*B — draft for review by A. Numbers from all model_validations/ scripts.*

---

## Table: Phenomenological Predictions Across Benchmark Points

| Section | Observable | BP1 ($m_\chi = 20.7$ GeV, $\lambda = 1.9$) | BP9 ($m_\chi = 37.9$ GeV, $\lambda = 4.3$) | MAP ($m_\chi = 94.07$ GeV, $\lambda = 48.6$) | Observation / Constraint |
|---|---|---|---|---|---|
| §7.1 | CP band $\alpha_s/\alpha_p$ | 13 — 212 (1.22 dec) | — | **1.8 — 11,532 (3.81 dec)** | Testable at colliders |
| §7.1 | $\lambda$ range in band | 2.4 — 9.9 | — | 3.3 — 261.8 | — |
| §7.2 | Fornax GC survival | ⚠ 12/15 | ⚠ 11/15 | ✅ **14/15** | All 5 GCs survive (observed) |
| §7.2 | SIDM core $r_{\rm core}$ (Fornax) | 449 pc | 332 pc | **829 pc** | $\gtrsim 500$ pc (Walker+2011) |
| §7.3 | RAR $\Upsilon_*$ (spirals) | 1.2 — 1.6 | — | 1.4 | 0.5 — 2.0 (McGaugh+2016) |
| §7.3 | RAR scatter (spirals) | 0.122 dex | — | — | 0.177 dex (observed) |
| §7.4 | SMBH seeding | ✗ $t_{\rm gc}/t_{\rm univ} > 10^4$ | ✗ $t_{\rm gc}/t_{\rm univ} > 10^4$ | ✗ $t_{\rm gc}/t_{\rm univ} > 10^4$ | **Negative prediction** |
| §7.5 | Classical dSph cores ($t_{\rm age} = 10$ Gyr) | ✗ 0/8 | ✗ 0/8 | ✅ **5/8** (7/8 at 12 Gyr) | Observed cores in Fornax, Sculptor, Carina, Sextans |
| §7.5 | UFD cores | ✗ 0/6 | ✗ 0/6 | ✅ **2/6** | Limited data (future surveys) |
| §7.5 | Crater II $r_{\rm core}$ | 27 pc | 32 pc | ⚠ 180 pc | $r_{\rm half} = 1066$ pc (tidal+SIDM) |
| §7.5 | $r_{\rm core}/r_{\rm half}$ universality | NOT universal (86%) | NOT universal (82%) | **NOT universal (67%)** | Testable prediction |
| — | $\sigma/m(30$ km/s) | 0.52 | 0.32 | **1.71** | 0.1 — 10 cm²/g |
| — | $\sigma/m(1000$ km/s)$^\dagger$ | 0.072 | 0.040 | 0.203 | $< 1$ cm²/g (Bullet Cluster) |
| — | Low-$v$ plateau? | ✗ No (92% variation) | — | ✅ Yes (12% variation) | — |
| — | Relic density $\Omega h^2$ | 0.120 | 0.120 | 0.120 | 0.120 ± 0.001 (Planck) |
| — | VPM regime | Born ($\lambda = 1.9$) | Classical ($\lambda = 4.3$) | Resonant ($\lambda = 48.6$) | Determines $\sigma/m(v)$ shape |
| — | $\Delta N_{\rm eff}$ | $\approx 0$ | $\approx 0$ | $\approx 0$ | $< 0.3$ (Planck) |

---

## Narrative Summary

The dual-benchmark strategy is strongly justified by the data:

**MAP** ($\lambda = 48.6$, resonant regime) is the **astrophysical champion**:
- Produces observable cores in most dwarf spheroidals (5/8 classical, 2/6 UFDs)
- Achieves near-perfect Fornax GC survival (14/15)
- Shows true $\sigma/m$ plateau at low velocities (12% variation across $v = 1$–$12$ km/s)
- Spans 3.81 decades of CP violation — the widest viable band

**BP1** ($\lambda = 1.9$, Born-to-classical transition) is the **collider-accessible benchmark**:
- Lower $m_\chi = 20.7$ GeV — potentially accessible to future dark matter searches
- Existing CP-separation band (1.22 decades) demonstrates the effect exists even at low $\lambda$
- Marginal Fornax GC survival motivates the full MAP analysis

**Both benchmarks** satisfy relic density, cluster constraints, and $\Delta N_{\rm eff}$ by construction.

**Structural predictions:**
1. SMBH seeds do NOT form via gravothermal collapse — structural consequence of cluster-safe $\sigma/m(v)$
2. $r_{\rm core}/r_{\rm half}$ is NOT universal in velocity-dependent SIDM — distinguishes from constant-$\sigma/m$ models
3. Crater II requires tidal processing in addition to SIDM — not a pure SIDM test
4. The partial-wave weights for identical Majorana scattering (singlet: even-$\ell$, weight 1; triplet: odd-$\ell$, weight 3) produce a $\sigma_T(v)$ profile quantitatively distinguishable from Dirac+scalar models — measurements of $\sigma/m$ at three or more velocity scales could discriminate between the two

---

$^\dagger$ All viable points satisfy $\sigma/m(1000) < 1$ cm$^2$/g (Bullet Cluster bound); the MAP benchmark gives $\sigma/m(1000) = 0.203$ cm$^2$/g. CP band maximum $\sigma/m(1000) = 0.28$ cm$^2$/g (at resonance peak).

*RAR v2 (gas-fraction fix) validated — see §7.3. BP1 dwarfs: 0.179 dex scatter (37% improvement). MAP overcores at all scales.*
