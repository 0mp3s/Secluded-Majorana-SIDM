# יומן מחקר V10: SIDM Majorana Dark Matter — Secluded Scalar Mediator
**Self-Interacting Dark Matter from a Majorana Fermion with Light Scalar Mediator (Secluded)**

מחבר: עומר פ. | מרץ 2026  
כלים: Python 3.14 + Numba 0.64 + NumPy | ביקורת עמיתים: Claude Opus, Gemini, GPT

---

## סיכום V9 → V10: מעבר ל-Secluded Dark Sector

### הסיבה למעבר

V9 הוגש לביקורת עמיתים. ביקורת זיהתה בעיה **קריטית C1**:

> **נוסחת σ_SI עם Higgs portal** מכילה enhancement של 1/m_φ⁴.
> עבור m_φ = 11 MeV, הדרישה מ-LZ נותנת sinθ < 6×10⁻¹⁰,
> בעוד BBN דורש sinθ > 2×10⁻⁵.
> **פער של ×30,000** — Higgs portal סגור לחלוטין.

**הפתרון:** אימוץ **secluded dark sector** (sinθ = 0) — זה הפך ל-V10.

### מה השתנה ב-V10

| Feature | V9 (Higgs Portal) | V10 (Secluded) |
|---------|-------------------|----------------|
| sinθ | Tunable | **= 0 exactly** |
| σ_SI | Testable via DD | **= 0** (prediction) |
| φ stability | Decays φ→e⁺e⁻ | **Stable** |
| BBN filter | m_φ > 2m_e required | Not needed |
| ΔN_eff | From decay products | 0.027 (stable φ) |
| CMB constraint | Applies (s-wave) | **Evaded** (χχ→φφ) |
| Parameters | (m_χ, m_φ, α, sinθ) | **(m_χ, m_φ, α)** only |

### מה נשאר זהה

- **VPM Solver**: ללא שינוי — אותו V(r) = -α e^{-m_φ r}/r
- **Parameter scan**: 80,142 raw viable, 646 representative
- **17 cosmological BPs**: אותה "island of viability"
- **Error budget**: 2-7% at 30 km/s
- **Blind tests**: כולם עוברים

---

## v32: Literature Cross-Check (Born Regime)

### מטרה
אימות VPM solver מול Born phase shifts מחושבים באופן עצמאי (scipy quadrature).

### תוצאות: 6/6 PASS

| Test | תיאור | תוצאה |
|------|--------|--------|
| T1 Born | VPM vs Born partial-wave sum (κ∈[0.1,0.5]) | **PASS** — 85-97% agreement |
| T2 Classical | VPM in classical regime (κ>10) | **PASS** — median ratio ~0.75 |
| T3 Velocity | σ/m decreases with velocity | **PASS** — monotonic decrease |
| T4 Unitarity | |S_l| ≤ 1 for all partial waves | **PASS** |
| T5 Maj/Dirac | σ_Maj/σ_Dirac = 4 | **PASS** |
| T6 BP1 | Consistency with preprint values | **PASS** — 0.9% at 30 km/s, 0.4% at 1000 km/s |

### הערות
- T1 Born: הסטייה של ~15% ב-κ=0.1 נובעת מ-barrier cutoff ביוקאווה שנעדר ב-Born
- תיקון: גרסה ראשונה השוותה total elastic (VPM) מול transfer cross section (TYZ Eq. 11) — תפוזים מול תפוחים. תוקן לשימוש באותה נוסחת partial-wave sum

---

## v33: השוואה לנתונים תצפיתיים

### מטרה
להשוות עקומות σ/m(v) תיאורטיות מול נתוני תצפית אמיתיים מ-5 מאמרים מפורסמים.

### מקורות נתונים

| מקור | תצפיות | v range |
|------|---------|---------|
| Kaplinghat, Tulin & Yu (2016) PRL 116, 041302 | 8 systems: dSphs, spirals, groups, clusters | 12–1200 km/s |
| Kamada et al. (2017) PRL 119, 111102 | Diverse rotation curves | ~40 km/s |
| Randall et al. (2008) ApJ 679 | Bullet Cluster | 4700 km/s |
| Harvey et al. (2015) Science 347 | 72 cluster mergers | 1500 km/s |
| Elbert et al. (2015) MNRAS 453 | TBTF dwarfs | 30 km/s |

### תוצאות: BP1 תואם 11/13 תצפיות

| System | v [km/s] | Theory σ/m | Obs range | Status |
|--------|----------|-----------|-----------|--------|
| Draco dSph | 12 | 0.430 | [0.1, 2.0] | ✅ |
| Fornax dSph | 12 | 0.430 | [0.2, 3.0] | ✅ |
| NGC 2976 | 60 | 0.498 | [0.5, 5.0] | ⚠️ barely below (0.4%) |
| NGC 1560 | 55 | 0.504 | [1.0, 8.0] | ❌ below lower bound |
| IC 2574 | 50 | 0.509 | [0.3, 5.0] | ✅ |
| NGC 720 (group) | 250 | 0.325 | [0.1, 1.5] | ✅ |
| NGC 1332 (group) | 280 | 0.307 | [0.05, 1.0] | ✅ |
| Abell 611 | 1200 | 0.052 | [0.02, 0.3] | ✅ |
| Abell 2537 | 1100 | 0.061 | [0.03, 0.4] | ✅ |
| Diverse RC band | 40 | 0.516 | [0.5, 10.0] | ✅ |
| Bullet Cluster | 4700 | 0.004 | [0.0, 1.25] | ✅ |
| 72 cluster mergers | 1500 | 0.035 | [0.0, 0.47] | ✅ |
| TBTF dwarfs | 30 | 0.515 | [0.5, 5.0] | ✅ |

### ניתוח
- BP1 מייצר **velocity dependence נכון**: σ/m ≈ 0.5 cm²/g ב-dwarfs, יורד דרך groups ל-≈0.05 ב-clusters
- NGC 1560 דורש σ/m ≥ 1.0 — BP1 נותן 0.504. זה system בעל אי-ודאויות גדולות (~0.5 dex) מ-halo modeling
- NGC 2976 shortfall של 0.4% בלבד — בתוך שגיאת המדידה
- BP5 ו-BP17 מראים velocity profile דומה — לא ייחודי ל-BP1

### Figure 3
נוצר: v33_observational_comparison.png — 3 עקומות BP + 13 נקודות data עם error bars

---

## v34: התאמת χ² לנתונים תצפיתיים

### מטרה
מעבר להשוואה איכותית (§4.5) — ביצוע התאמת χ² מלאה ל-13 מערכות אסטרופיזיקליות עם full VPM solver לכל נקודה.

### מתודולוגיה
- 5,009 נקודות מדגמות מ-all_viable_raw_v8.csv (80,142) + 17 relic BPs = 5,026 נקודות
- 12 מהירויות ייחודיות = 60,312 קריאות VPM מלאות
- מקבילות: 12 עובדים מתוך 14 ליבות (`multiprocessing.Pool`, `imap_unordered`, `chunksize=32`)
- Fine-scan סביב 50 הנקודות הטובות ביותר — 400 הערכות נוספות
- זמן ריצה: 697 שניות (~11.5 דקות)
- שגיאות אסימטריות: σ⁺ ו-σ⁻ לפי סימן השארית (residual)

### באג יחידות שתוקן
`sigma_T_vpm` מקבל m_φ ב-**GeV**. ה-CSV של raw scan (קובץ all_viable_raw_v8.csv) שומר m_φ ב-GeV, אבל ה-CSV של relic BPs (v31_true_viable_points.csv) שומר ב-MeV.
גרסאות קודם המירו MeV ל-VPM — החזירו σ/m = 0 לכל הנקודות.
תיקון: raw CSV נשאר כמות שהוא (GeV); relic CSV מחולק ב-1000 (MeV→GeV).
אימות: `sigma_T_vpm(20.69, 11.34e-3, 1.048e-3, 30.0)` = 0.515 cm²/g ✓

### תוצאות: התאמה מלאה

**Free best-fit (ללא אילוץ relic):**

| פרמטר | ערך |
|---------|------|
| χ²/dof | **2.59/10 = 0.26** |
| m_χ | 100 GeV |
| m_φ | 9.15 MeV |
| α | 2.84×10⁻³ |
| λ | 31.04 |
| σ/m(30) | 1.91 cm²/g |
| σ/m(1100) | 0.11 cm²/g |
| גדול ביותר pull | Bullet Cluster: −0.99σ |

**Relic-constrained best-fit (BP1):**

| פרמטר | ערך |
|---------|------|
| χ²/dof | **5.40/10 = 0.54** |
| m_χ | 20.69 GeV |
| m_φ | 9.91 MeV |
| α | 1.048×10⁻³ |
| λ | 2.19 |
| σ/m(30) | 0.79 cm²/g |
| גדול ביותר pull | NGC 1560: −1.12σ |

**כל 17 נקודות relic:** χ²/dof בין 0.54 ל-0.84 — **כולן** מתחת 1.0

### ניתוח
- χ²/dof < 1 משקף את אי-הוודאויות הגדולות בנתונים האסטרופיזיקליים — המודל לא "over-fits"
- התוצאה העיקרית: סט אחד של פרמטרים מתאים **בו-זמנית** לכל 13 המערכות (12–4700 km/s)
- אילוץ relic density לא פוגע באיכות ההתאמה: ה-BPs עדיין נותנים χ²/dof < 0.85
- NGC 1560 ה-pull הגדול ביותר (−1.12σ) — עדיין מקובל לחלוטין

### Figure 4
נוצר: v34_chi2_fit.png — עקומות best-fit + 13 נקודות data עם error bars

### קבצים
- v34_results.csv (5,426 שורות): כל הנקודות עם χ², σ/m בכל 12 מהירויות
- v34_chi2_fit.png: גרף
- v34_output.txt: פלט מלא

---

## סטטוס המאמר

### מבנה מעודכן
- §1–§3: Model + VPM (ללא שינוי)
- **§4.5**: Comparison with Observational Data — 13 data points מ-5 מאמרים
- **§4.6 חדש**: Quantitative χ² Fit — χ²/dof = 0.26 (free), 0.54 (relic BP1), Table 4 עם 17 BPs
- §5: Secluded Dark Sector (ללא שינוי)
- §6–§7: Relic Density + Predictions (ללא שינוי)
- **§8.2 מעודכן**: נקודה 5 מעודכנת עם תוצאות χ²
- **Appendix A.4**: Literature cross-check (Born regime, 6/6 PASS)
- **References [20]–[22]**: Randall+08, Harvey+15, Elbert+15

### Figures
1. v31_island_of_viability.png — Island of 17 BPs
2. v31_bp1_velocity_profile.png — BP1 σ/m(v)
3. v33_observational_comparison.png — BP1/BP5/BP17 vs real data
4. **v34_chi2_fit.png** — Best-fit σ/m(v) עם 13 data points + error bars (NEW)

### Validation Summary
| Script | Purpose | Result |
|--------|---------|--------|
| v21 | Born validation | PASS |
| v23 | Error budget | 2-7% at 30 km/s |
| v25 | Mediator cosmology | 5/5 PASS |
| v26 | BSF | 6/6 PASS |
| v27 | Boltzmann solver | 6/6 PASS |
| v28/v29 | Blind sanity | 3/4 + 3/3 PASS |
| v30 | Benchmark extraction | PASS |
| v31 | Cosmological scan | 17 viable BPs |
| **v32** | **Literature cross-check** | **6/6 PASS** |
| **v33** | **Observational comparison** | **11/13 compatible** |
| **v34** | **χ² fit to 13 systems** | **χ²/dof = 0.26 (free), 0.54 (relic)** |
| **v35** | **VPM vs Born transfer cross-check** | **VPM/Born ≈ 0.8–0.9 in Born regime (κ > 0.1)** |
| **v36** | **Sommerfeld enhancement** | **S < 1.026 at freeze-out, tree-level valid** |

---

## v35: השוואת VPM מול Born transfer cross section

### מטרה
בדיקת תקפות סולבר ה-VPM ע"י השוואה לחתך פעולה אנליטי Born transfer עבור פרמיוני מיורנה זהים.

### שיטה
חישוב אינטגרל Born מדויק עם:
- אמפליטודה מסומטרת של מיורנה: $\frac{1}{4}|f+f'|^2 + \frac{3}{4}|f-f'|^2$
- משקל momentum transfer: $(1-\cos\theta)$
- מקדם $1/2$ לחלקיקים זהים
- אינטגרציה נומרית ב-scipy.integrate.quad

### תוצאות

**Born regime ($\lambda < 1$, $\kappa > 0.1$):**
| נקודה | $\lambda$ | $\kappa$ | VPM/Born |
|:-------|:---------:|:--------:|:--------:|
| Born-4 ($m_\chi=100, v=200$) | 0.02 | 0.17 | 0.91 |
| Born-1 ($m_\chi=50, v=100$) | 0.05 | 0.08 | 0.78 |

**Resonant regime ($\lambda > 1$):**
| נקודה | $\lambda$ | $\kappa$ | VPM/Born |
|:-------|:---------:|:--------:|:--------:|
| BP1 ($v=1000$) | 2.19 | 0.35 | 0.80 |
| BP1 ($v=30$) | 2.19 | 0.10 | 0.26 |
| Free-best ($v=30$) | 31.04 | 0.55 | 0.02 |

### מסקנות
1. ב-Born regime עם $\kappa > 0.1$: VPM/Born ≈ 0.8–0.9. הפער ~20% הוא ההבדל הידוע בין elastic ל-transfer cross section
2. ב-resonant regime ($\lambda \gg 1$): Born נכשל לחלוטין — VPM/Born ~0.02. מאשר הצורך ב-VPM מלא
3. הוספנו הערה ב-§3.2 של ה-preprint + Appendix C עם תוצאות מפורטות

### קבצים
- v35_tyz_comparison.py: סקריפט ההשוואה
- v35_tyz_comparison.png: גרף דו-פאנלי (scatter + velocity scan)

---

## v36: Sommerfeld Enhancement — חישוב נומרי

### מטרה
חישוב נומרי מלא של ה-Sommerfeld enhancement factor $S_0$ לאניהילציה s-wave, עבור פוטנציאל Yukawa אטרקטיבי, בכל 17 benchmark points.

### שיטה
פתרון משוואת שרדינגר רדיאלית $l=0$ ביחידות חסרות מימד $x = m_\phi r$:
$$u''(x) + [\kappa^2 + \lambda e^{-x}/x] u(x) = 0$$
- RK4 integration עבור $\kappa < 10$
- Coulomb analytic upper bound עבור $\kappa \geq 10$ (tight bound כש-$\kappa \gg \lambda$)
- חילוץ אמפליטודה: $S_0 = 1/(\kappa^2 A^2)$ עם $A = \sqrt{u^2 + (u'/\kappa)^2}$

### Validation
1. **Free particle** ($\alpha \to 0$): $S = 1.000$ בכל מהירויות ✓
2. **Perturbative limit** (high v): $S \approx 1 + \pi\alpha/v$ ✓
3. **Yukawa vs Coulomb** ($\kappa \gg 1$): ratio = 1.0000 ✓

### תוצאות

**Freeze-out ($v \sim 0.3c$):**
- $S_0$ range: 1.003–1.025
- Maximum relic correction: $|\Delta\Omega h^2 / \Omega h^2| < 2.5\%$
- Tree-level Boltzmann solver תקף ✓

**Late-time ($v \sim 30$ km/s):**
- $S_0$ range: 5.7–417
- Depletion timescale: $\tau \sim 10^7$–$10^8\,t_{\rm Hubble}$ — אירלוונטי
- Secluded model: אין CMB/Fermi-LAT constraints

### עדכונים ל-preprint
- §6.5 חדש: "Sommerfeld Enhancement" — חישוב מלא + טבלה + ניתוח
- §8.2 point 2: הוספת "Sommerfeld enhancement is negligible ($S_0 < 1.026$)"

### קבצים
- v36_sommerfeld.py: סקריפט החישוב (RK4 + Coulomb hybrid)
- v36_sommerfeld.png: גרף דו-פאנלי ($S(v)$ + $S$ vs $\lambda$ for all BPs)
- v36_output3.txt: פלט מלא

---

## v37: מיצוע מהירויות Maxwell–Boltzmann

### מטרה
בדיקה: כמה משנה מיצוע על התפלגות מהירויות MB במקום הערכה במהירות בודדת $v_{\rm char}$?

### שיטה
עבור כל מערכת תצפיתית ($v_{\rm char}$), חישוב:
$$\langle\sigma/m\rangle_{\rm MB} = \frac{\int_0^\infty (\sigma/m)(v)\;v^3\,e^{-v^2/2v_0^2}\,dv}{\int_0^\infty v^3\,e^{-v^2/2v_0^2}\,dv}$$
עם $v_0 = v_{\rm char}/\sqrt{2}$ (dispersal 1D → most probable relative speed).
אינטגרציה: Gauss–Legendre quadrature, $n=30$ nodes על $[1, 5v_{\rm char}]$.

### תוצאות

**סטייה per-system (BP1, $\lambda=1.91$):**
- גמדים ($v \lesssim 60$ km/s): $\lesssim 5\%$ — cross section שטוח, מיצוע כמעט לא משנה
- קבוצות/clusters ($v \sim 250$–$1500$ km/s): $\sim$12–17% reduction — הזנב של MB דוגם $\sigma/m$ נמוך יותר
- Bullet Cluster ($v = 4700$ km/s): $\sim$2.5% — כבר ב-power-law regime

**השפעה על $\chi^2$ (כל 17 BPs):**
- שינוי $\chi^2$: 2–9%, ממוצע 6%
- Fixed $\chi^2/\nu$: min=0.42, max=0.65, mean=0.53
- Averaged $\chi^2/\nu$: min=0.45, max=0.69, mean=0.56
- **17/17 BPs נשארים עם $\chi^2/\nu < 1$** — אף BP לא משנה סטטוס

### מסקנה
הקירוב של מהירות בודדת מספיק ברמת הדיוק הנוכחית של אילוצי SIDM ($\sim$0.5 dex).
שינוי $\chi^2$ של $\sim$6% הוא זניח לעומת אי-ודאויות תצפיתיות.

### עדכונים ל-preprint
- §4.6: פסקה חדשה "Effect of velocity averaging"
- §8.2 point 5: הוספת "MB velocity averaging shifts χ² by ~6%"
- Appendix D חדש: שיטה + טבלה + תוצאות $\chi^2$ + Figure D1

### קבצים
- v37_velocity_averaged.py: סקריפט (Gauss–Legendre MB averaging)
- v37_velocity_averaged.png: Figure D1 (BP1 fixed vs averaged + ratio)
- v37_output2.txt: פלט מלא

---

## v38: MCMC Posterior Sampling — Bayesian Constraints

### מטרה
דגימת פוסטריור בייסיאנית של מרחב הפרמטרים $(m_\chi, m_\phi, \alpha)$ תוך שימוש ב-13 אילוצי תצפית. 
מטרות: (1) corner plot עם contours של 68% ו-95%; (2) אילוצים כמותיים על הפרמטרים; (3) אימות שכל 17 relic BPs נמצאים בתוך ה-posterior.

### שיטה
- **Sampler:** emcee (affine-invariant ensemble sampler)
- **Walkers:** 32
- **Burn-in:** 300 steps, **Production:** 2,000 steps
- **Total evaluations:** 73,600 likelihood calls
- **Workers:** 12 processes מתוך 14 ליבות (multiprocessing.Pool)
- **פריורים:** uniform ב-log₁₀ space:
  - $\log_{10}(m_\chi/\text{GeV}) \in [\log_{10}(5), \log_{10}(200)]$
  - $\log_{10}(m_\phi/\text{MeV}) \in [\log_{10}(3), \log_{10}(30)]$
  - $\log_{10}\alpha \in [\log_{10}(10^{-5}), \log_{10}(0.05)]$
- **Likelihood:** $\ln\mathcal{L} = -\chi^2/2$ עם asymmetric σ (same as v34)
- **Seeds:** 17 relic BPs + v34 unconstrained best fit, perturbation $\mathcal{N}(0, 0.05)$ ב-log space
- **Seed:** `np.random.default_rng(42)` — reproducible

### אבחון התכנסות
- **Acceptance fraction:** 0.522 (optimal range 0.2–0.5, מעט גבוה אבל תקין)
- **Max autocorrelation time:** 73.3 steps
- **Effective samples:** ~874 (64,000 / 73.3)
- **Chain traces:** רואים mixing טוב ב-chains plot — אין drift מתמשך

### תוצאות

**Best-fit (MAP):**

| פרמטר | ערך |
|---------|------|
| $m_\chi$ | 90.64 GeV |
| $m_\phi$ | 13.85 MeV |
| $\alpha$ | 2.546×10⁻² |
| $\lambda = 2\alpha m_\chi / m_\phi$ | 333.37 |
| $\chi^2$ | 1.575 |
| $\chi^2/\text{dof}$ | **0.1575** (dof = 10) |

**אילוצים (median ± 68% CI):**

| פרמטר | Median | 16th %ile | 84th %ile |
|---------|--------|-----------|-----------|
| $m_\chi$ [GeV] | 72.98 | 16.97 | 115.91 |
| $m_\phi$ [MeV] | 10.83 | 6.51 | 15.59 |
| $\alpha$ | 0.004 | 0.001 | 0.022 |
| $\lambda$ (derived) | 59.0 | 3.6 | 272.2 |

**Relic BPs vs posterior — 17/17 within 95% CI:**

| BP | $m_\chi$ | $m_\phi$ [MeV] | $\alpha$ | $\lambda$ | $\chi^2$ | within 95% |
|----|----------|-----------------|----------|-----------|----------|------------|
| 1 | 20.7 | 11.34 | 1.048e-03 | 3.83 | 8.418 | ✅ |
| 2 | 29.8 | 12.98 | 1.473e-03 | 6.76 | 8.282 | ✅ |
| 3 | 16.2 | 9.91 | 8.414e-04 | 2.76 | 7.796 | ✅ |
| 4 | 23.4 | 11.34 | 1.173e-03 | 4.83 | 7.203 | ✅ |
| 5 | 33.6 | 12.98 | 1.654e-03 | 8.56 | 7.507 | ✅ |
| 6 | 12.7 | 8.66 | 6.796e-04 | 2.00 | 7.448 | ✅ |
| 7 | 18.3 | 9.91 | 9.386e-04 | 3.47 | 6.389 | ✅ |
| 8 | 26.4 | 11.34 | 1.314e-03 | 6.11 | 6.304 | ✅ |
| 9 | 37.9 | 12.98 | 1.858e-03 | 10.86 | 6.951 | ✅ |
| 10 | 88.6 | 14.85 | 4.262e-03 | 50.85 | 6.562 | ✅ |
| 11 | 100.0 | 14.85 | 4.806e-03 | 64.73 | 6.147 | ✅ |
| 12 | 10.0 | 7.56 | 5.538e-04 | 1.46 | 7.286 | ✅ |
| 13 | 78.5 | 14.85 | 3.781e-03 | 39.96 | 6.979 | ✅ |
| 14 | 14.4 | 8.66 | 7.555e-04 | 2.51 | 5.814 | ✅ |
| 15 | 42.8 | 12.98 | 2.089e-03 | 13.78 | 6.495 | ✅ |
| 16 | 20.7 | 9.91 | 1.048e-03 | 4.38 | **5.396** | ✅ |
| 17 | 29.8 | 11.34 | 1.473e-03 | 7.73 | 5.692 | ✅ |

### ניתוח
1. **MAP vs relic BPs:** ה-MAP ($m_\chi \approx 91$ GeV) שונה מ-relic BPs ($m_\chi \sim 10$–$100$ GeV), אבל שניהם בתוך ה-posterior — ללא relic constraint, מסות גבוהות יותר מועדפות
2. **broad posterior:** ה-CI רחב ($m_\chi = 17$–$116$ GeV ב-68%) — הנתונים התצפיתיים עם אי-ודאויות $\sim$0.5 dex לא מצליחים לאלץ חזק
3. **λ posterior:** האטם bimodal/long-tailed, median ≈ 59, range 3.6–272 — resonant regime ($\lambda > 1$) מועדף
4. **17/17 BPs within 95%:** תוצאה מרכזית — ה-island of viability שלנו תואמת את ה-posterior הבייסיאני
5. **BP16 ($\chi^2 = 5.4$) הטוב ביותר מבין relic BPs** — $(m_\chi, m_\phi) = (20.7, 9.91)$ ב-center של ה-island

### הערות טכניות
- 874 effective samples מספיקים לאמידת contours ותיעוד credible intervals, אבל לא לאמידת tails באופן מדויק
- אם נרצה tails מדויקים יותר: $N_{\text{steps}} \times 3$–$5$ (6,000–10,000 production)
- ה-acceptance fraction (0.52) מעט גבוה — ייתכן שה-proposal stretch parameter (a=2 default) צריך עדכון, אבל לא קריטי

### עדכונים ל-preprint
- **§4.7 חדש:** "Bayesian Posterior Constraints" — MCMC methodology + corner plot + credible intervals
- **Table 5 חדש:** Parameter constraints table (median ± CI)
- **Figure 5 חדש:** v38_corner.png — corner plot with 68%/95% contours
- **Figure 6 חדש:** v38_lambda_posterior.png — posterior density of derived λ
- **§8.2 point 6:** הוספת "MCMC posterior sampling confirms all 17 BPs within 95% CI"

### קבצים
- run_mcmc.py (stats_mcmc/): סקריפט MCMC (emcee + corner + multiprocessing)
- output/v38_corner.png: Corner plot 3×3
- output/v38_mcmc_chains.png: Chain trace plot (convergence diagnostic)
- output/v38_lambda_posterior.png: Posterior density of λ
- output/v38_mcmc_samples.csv: Full chain (64,000 rows × 8 columns)

---

## סיכום סטטוס מעודכן (אחרי v38)

### Validation Summary

| Script | Purpose | Result |
|--------|---------|--------|
| v21 | Born validation | PASS |
| v23 | Error budget | 2-7% at 30 km/s |
| v25 | Mediator cosmology | 5/5 PASS |
| v26 | BSF | 6/6 PASS |
| v27 | Boltzmann solver | 6/6 PASS |
| v28/v29 | Blind sanity | 3/4 + 3/3 PASS |
| v30 | Benchmark extraction | PASS |
| v31 | Cosmological scan | 17 viable BPs |
| **v32** | **Literature cross-check** | **6/6 PASS** |
| **v33** | **Observational comparison** | **11/13 compatible** |
| **v34** | **χ² fit to 13 systems** | **χ²/dof = 0.26 (free), 0.54 (relic)** |
| **v35** | **VPM vs Born transfer** | **VPM/Born ≈ 0.8–0.9 in Born regime** |
| **v36** | **Sommerfeld enhancement** | **S < 1.026 at freeze-out** |
| **v37** | **MB velocity averaging** | **Δχ² < 9%** |
| **v38** | **MCMC posterior** | **χ²/dof = 0.16 (MAP), 17/17 BPs in 95% CI** |

### Figures עבור ה-preprint

1. v31_island_of_viability.png — Island of 17 BPs
2. v31_bp1_velocity_profile.png — BP1 σ/m(v)
3. v33_observational_comparison.png — BP1/BP5/BP17 vs real data
4. v34_chi2_fit.png — Best-fit σ/m(v) עם 13 data points + error bars
5. **v38_corner.png** — MCMC corner plot עם 68%/95% contours (NEW)
6. **v38_lambda_posterior.png** — Posterior density of λ (NEW)
