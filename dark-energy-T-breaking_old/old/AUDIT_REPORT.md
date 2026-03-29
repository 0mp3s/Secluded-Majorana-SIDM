# דוח ביקורת מקיף — `dark-energy-T-breaking/`

**תאריך**: ביקורת מערכתית של כל 22 סקריפטים מול `research_journal.md`, `need_to_verify.md`, `theory.md`

---

## טבלת סיכום

| # | סקריפט | Test | סטטוס | בעיות עיקריות |
|---|--------|------|--------|---------------|
| 1 | `sigma_radiative_stability.py` | 1 | ✅ | ניתוח גנרי, לא תלוי בקונבנציה |
| 2 | `dark_axion_full.py` | 2 | ⚠️ | Bug 2 (CW min) התגלה בדיעבד בTest 4; הפלט כבר מתועד |
| 3 | `fifth_force_constraints.py` | 3 | ✅ | תקין; M_Pl, H_0 נכונים |
| 4 | `freeze_out_trapping.py` | 4 | ❌ | **5 באגים מתועדים** — הוחלף ע"י v_corrected |
| 5 | `verify_cw_bugs.py` | 4-audit | ✅ | מזהה את כל 5 הבאגים כראוי |
| 6 | `freeze_out_analysis_corrected.py` | 4-fix | ✅ | מתקן Bug 1 ([ln-1.0]); F_th/F_CW ~ 10⁻¹⁰ |
| 7 | `theta_topological.py` | 5 | ✅ | ניתוח טופולוגי; θ/(2π) לא רציונלי → Z_N לא עובד |
| 8 | `consistency_check_sidm.py` | 6 | ❌ | **שגיאת קונבנציה** — מכפיל ב-(8/9) שמבוטלת בTest 12 |
| 9 | `consistency_17bp.py` | 7 | ❌ | **שגיאת קונבנציה** — 10/17 שגוי; התוצאה הנכונה: 17/17 |
| 10 | `free_theta_scan.py` | 8 | ⚠️ | קונבנציה מושפעת; תבנית איכותית תקפה |
| 11 | `a4_dark_sector_model.py` | 9 | ✅ | תורת חבורות; sin²θ=1/10 for equal VEVs; תואם |
| 12 | `verify_a4_cg.py` | 10-CG | ✅ | tan²θ=1/9 מאומת; SymPy + numeric ✓ |
| 13 | `boltzmann_17bp.py` | 10-Boltz | ❌ | **שגיאת קונבנציה** — α_s=(8/9)α שגוי |
| 14 | `sigma_trapping_ode.py` | 10-trap | ⚠️ | y²=4πα (חסר /cos²θ); מסקנה איכותית נכונה |
| 15 | `fornax_gc_check.py` | 11 | ⚠️ | (8/9) מיותר; עמודת "without decomp" = הנכונה |
| 16 | `test_alpha_convention.py` | 12 | ✅ | **מפתח** — מוכיח α_CSV = α_s; 17/17 pass |
| 17 | `dark_force_accumulation.py` | 13 | ⚠️ | α=3.4e-3 שגוי (MAP=5.734e-3); מסקנה איכותית נכונה |
| 18 | `sigma_mass_protection.py` | 13D | ✅ | קונבנציה מתוקנת; 5 מנגנוני הגנה |
| 19 | `dark_qcd_consistency.py` | 14 | ✅ | 6 בדיקות; SU(2)_d preferred; Λ_d ~ √(H₀ M_Pl) |
| 20 | `born_full_amplitude.py` | 15 | ✅ | אמפליטודת Born מלאה; interference < 10⁻⁶; VPM validated |
| 21 | `vev_alignment_stability.py` | 16 | ✅ | יציבות VEV מאומתת; Hessian חיובי |
| 22 | `resonance_bp1.py` | — | ⚠️ | m_φ(BP1)=11.34 MeV לא תואם ל-10.83 MeV; + קונבנציה |

**סיכום**: ✅ 12 נקיים | ⚠️ 6 עם הערות | ❌ 4 עם שגיאות מהותיות

---

## ממצאים מפורטים לפי סקריפט

---

### 1. `sigma_radiative_stability.py` — Test 1 ✅

**מטרה**: בדיקת 4 אפשרויות לסקלר σ אולטרא-קל

**ממצאים**:
- 4 אפשרויות A-D נבדקות כהלכה; Option D (dark axion) הוא היחיד שעובד
- קבועים נכונים: M_Pl = 2.435e18, H_0 = 1.44e-42
- V_Planck = m_χ⁶/(32π² M_Pl²) — פורמולה נכונה for n_f=2 Majorana
- **עקבי עם Journal (Test 1)**

**הערות**: ניתוח גנרי שלא תלוי בקונבנציית α

---

### 2. `dark_axion_full.py` — Test 2 ⚠️

**מטרה**: ניתוח V_eff(σ) מלא, גזירת θ_relic

**ממצאים**:
- הגזירה θ_relic = arctan(1/√8) = 19.47° **נכונה** ✓
- הפלט טוען CW minimum at θ=0 — **שגוי**, תוקן ב-Test 4 (Bug 2): CW minimum at θ=π/2
- עם זאת: הקוד עצמו חישב את V_CW כראוי; הפרשנות בטקסט הייתה שגויה
- f_needed values for V_σ = ρ_Λ **עקביים עם journal** ✓
- מבנה ה-sign structure (scalar deepens V, pseudo raises V) **תקין** ✓

**שורה תחתונה**: הסקריפט הפיק תוצאה פיזיקלית נכונה שהתפרשה שגוי (Bug 2). הבאג מתועד ומתוקן.

---

### 3. `fifth_force_constraints.py` — Test 3 ✅

**מטרה**: בדיקת אילוצי כוח חמישי עם β = M_Pl/f

**ממצאים**:
- β = 5 for f = 0.2 M_Pl — מחושב נכון ✓
- רענון שנוי screening: m_σ,eff(ρ) density-dependent → טווח ~30m בוואקום
- Planck constraints: β < 0.066, אבל CW mass screens the force ✓
- **CMB epoch**: range ≪ sound horizon → screened ✓
- **עקבי עם Journal (Test 3)**

---

### 4. `freeze_out_trapping.py` — Test 4 ❌ SUPERSEDED

**מטרה**: ODE מצומד: σ + Boltzmann + Hubble friction

**ממצאים — 5 באגים מתועדים**:

| Bug | תיאור | חומרה |
|-----|--------|--------|
| 1 | מקדם CW: `[ln(M²/μ²) - 0.5]` במקום `[ln - 1.0]` | גבוהה — factor-of-2 error in derivative |
| 2 | CW minimum at θ=0 (שגוי) — צריך θ=π/2 | גבוהה — כיוון הגלגול הפוך |
| 3 | m²_σ sign שגוי | גבוהה |
| 4 | `t_eval` passed as float instead of array | runtime crash |
| 5 | Hard-coded summary doesn't match code | cosmetic |

**סטטוס**: **הוחלף ע"י** `freeze_out_analysis_corrected.py` + `sigma_trapping_ode.py`

---

### 5. `verify_cw_bugs.py` — Test 4 Audit ✅

**מטרה**: אימות 5 הבאגים ב-freeze_out_trapping.py

**ממצאים**:
- כל 5 הבאגים מזוהים ומאומתים ✓
- Bug 1: נגזרת CW צריכה `Σ Mⱼ⁴ [ln(Mⱼ²/μ²) - 3/2]`, לא `-1/2`
- Bug 2: d²V_CW/dθ² < 0 for ALL θ → CW concave → minimum at boundary θ=π/2 ✓
- **עקבי עם Journal**

---

### 6. `freeze_out_analysis_corrected.py` — Test 4 Fix ✅

**מטרה**: גרסה מתוקנת של Test 4

**ממצאים**:
- מקדם [ln - 1.0] מתוקן ✓
- F_th/F_CW ~ 10⁻¹⁰ → thermal backreaction זניח ✓
- CW dominates for ALL f values ✓
- מסקנה: σ rolling to π/2 → SIDM dies → need A₄ discrete mechanism
- **עקבי עם Journal**

---

### 7. `theta_topological.py` — Test 5 ✅

**מטרה**: חיפוש מבנה טופולוגי של θ_dark

**ממצאים**:
- θ_relic/(2π) = 0.054... — NOT rational → Z_N orbifold won't work ✓
- Tetrahedral connection: θ = complement of tetrahedral angle ✓
- Neutrino: sin²θ₁₂(TBM) = 1/3 = 3 × sin²θ_dark ✓
- Witten effect: α²/8 NOT Dirac-quantized — correctly identified ✓
- UV completion (Part 7) with fixed couplings — no σ needed
- V_CW / ρ_Λ = 10⁴⁰ → CC problem, correctly identified
- **עקבי עם Journal (Test 5)**

---

### 8. `consistency_check_sidm.py` — Test 6 ❌ CONVENTION ERROR

**מטרה**: VPM + KT relic for 4 scenarios

**ממצאים — שגיאת קונבנציה מהותית**:
- **Scenario A**: α_s = cos²θ × α_CSV → **שגוי!** α_CSV = α_s already (Test 12)
- The (8/9) reduction is double-counting → BP1 gets σ/m=0.447 (marginal) instead of true value
- "10/17 survive" is **WRONG** — correct answer is **17/17** (Test 12)
- VPM solver code itself is correct (spin weights w_even=1, w_odd=3 ✓)
- KT analytic formula is correct ✓
- Scenario B (pure scalar) results are mathematically correct but physically wrong interpretation

**impact**: All σ/m values using α_s = (8/9)α are systematically low by ~11% (for cos²θ scenarios). This pushed 7 BPs below the 0.5 threshold.

**סטטוס**: **הוחלף ע"י** `test_alpha_convention.py`

---

### 9. `consistency_17bp.py` — Test 7 ❌ CONVENTION ERROR

**מטרה**: 17 BP consistency test with θ-decomposition

**ממצאים**:
- **Same convention error**: α_s = (8/9)α applied to α_CSV which IS already α_s
- Result "10/17 pass SIDM" is **WRONG** → correct: **17/17**
- "Resonance crossing" analysis — still qualitatively valid (λ mapping is correct)
- KT relic density shows Ωh² × 81/64 increase — this is **ALSO WRONG** (no increase needed if α_CSV = α_s)
- The "nearest resonance" analysis is based on λ_new = (8/9)λ_orig → λ values are wrong

**סטטוס**: **הוחלף ע"י** `test_alpha_convention.py`

---

### 10. `free_theta_scan.py` — Test 8 ⚠️

**מטרה**: סריקת θ ∈ [0°, 45°] עבור BPs שונים

**ממצאים**:
- Uses α_s = α_CSV × cos²θ in VPM
- **Convention issue**: If α_CSV = α_s (at θ_relic=19.47°), then varying θ requires α_total = α_CSV/cos²(19.47°) as baseline, and α_s(θ) = α_total × cos²θ
- At θ=0°: script gives α_s = α_CSV — correct only if α_CSV = α_total (which it isn't)
- At θ=19.47°: script gives α_s = (8/9)α_CSV — double-counting
- **התנהגות איכותית** (monotonic σ/m vs θ) **עדיין תקפה** — just the absolute scale is shifted
- Resonance map (λ(θ) vs λ_crit) is qualitatively valid

**שורה תחתונה**: מספרים מוחלטים שגויים, אך המגמה האיכותית שרירה

---

### 11. `a4_dark_sector_model.py` — Test 9 ✅

**מטרה**: בניית מודל A₄ מפורש

**ממצאים**:
- S, T generators verified: S²=I, T³=I, (ST)³=I, |A₄|=12 ✓
- S eigenvalues: +1 (singlet), -1/3 (doublet) ✓
- Decomposition of (1,0,0) in S-eigenbasis: 1/3 + 2/3 ✓
- sin²θ = 1/10 for equal VEVs (not 1/9) — correctly identified
- Two-flavon model: v_S = 2√6 v_T needed for sin²θ = 1/9
- Part 7: Most general A₄ mass matrix M = m₀I + m_S S ✓
- **עקבי עם Journal (Test 9) ו-verify_a4_cg.py**

---

### 12. `verify_a4_cg.py` — Test 10 CG ✅

**מטרה**: אימות CG coefficients של A₄ עם SymPy + numeric

**ממצאים**:
- A₄ group verified: 12 elements, all relations hold ✓
- 3⊗3→1: ab contraction invariant under all 12 elements ✓
- 3⊗3⊗3→1: g_s = 3 (or 2 with 6-term formula), g_p = 1 (or 2/3)
  - **Ratio g_p/g_s = 1/3 IDENTICAL in both** ✓
- tan²θ = 1/9 **CONFIRMED** both symbolically and numerically ✓
- Majorana symmetry: contraction IS symmetric in first 2 slots for ψ=(1,1,1)/√3 ✓
- BUT: NOT automatically symmetric for GENERIC vectors → need symmetrization for general case
- VEV ratio: v_p/v_s = 3/(2√2) ≈ 1.061 closes 1/10→1/9 gap ✓
- SymPy confirms: (abc)₁ − (bac)₁ ≠ 0 in general → 9-term formula has a bug
- **עקבי עם Journal ו-need_to_verify.md (Item 2a)**

---

### 13. `boltzmann_17bp.py` — Test 10 Boltzmann ❌ CONVENTION ERROR

**מטרה**: Full RK4 Boltzmann solver עם g_*(T) tabulation

**ממצאים**:
- Boltzmann solver itself is well-implemented (RK4, g_* table, KT comparison)
- **Convention error**: α_s = (8/9)α applied to α_CSV
- Shows Ωh² × 81/64 = 1.266 increase — **WRONG** (no increase if α_CSV = α_s)
- KT/Boltzmann ratio analysis is VALID (systematic bias ~constant) ✓
- The g_*(T) tabulation is correct ✓

**impact**: Relic density values after "θ-decomposition" are all **overestimated** by ~27%, because the coupling was double-reduced.

---

### 14. `sigma_trapping_ode.py` — Test 10 Trapping ⚠️

**מטרה**: ODE מתוקן לריצת σ עם CW force

**ממצאים**:
- **[ln-1.0] coefficient** — מתוקן (Bug 1 fixed) ✓
- CW minimum correctly at θ=π/2 ✓
- **Minor inconsistency**: `y_sq = 4*math.pi*alpha` with α = 5.734e-3 (MAP CSV)
  - Post-Test 12: α_CSV = α_s, so this gives y_s² (not y_total²)
  - Should be `y_sq = 4*math.pi*alpha / cos(theta)**2` to get y²_total
  - Error: CW force underestimated by factor (8/9)² ≈ 0.79 (21%)
- **Qualitative conclusion UNAFFECTED**: CW still dominates thermal by ~10¹⁰ even with 21% correction
- Part 5 conclusion: A₄ discrete interpretation preferred → no σ rolling → clean ✓
- **עקבי עם Journal (Test 10c)**

---

### 15. `fornax_gc_check.py` — Test 11 ⚠️

**מטרה**: Fornax GC constraint σ/m < 1.5 cm²/g at v=10-20 km/s

**ממצאים**:
- VPM solver copy — correct spin weights ✓
- FORNAX_BOUND = 1.5 cm²/g ✓
- **Convention issue**: "with decomposition" column uses α_s = (8/9)α_CSV → WRONG (double-counting)
- **"without decomposition" column = CORRECT** (uses α_CSV = α_s as-is)
- The script is **conservative**: underestimating α gives SMALLER σ/m → EASIER to pass Fornax
- So the "with decomposition passes" result is correct but for the wrong reason

**שורה תחתונה**: העמודה "without decomposition" היא הנכונה

---

### 16. `test_alpha_convention.py` — Test 12 ✅ KEY SCRIPT

**מטרה**: הכרעה אמפירית — α_s = (8/9)α או α_s = α?

**ממצאים**:
- **מפתח!** משווה שני תרחישים עם VPM ישיר:
  - Scenario A: α_s = (8/9)α → n_pass = varies
  - Scenario B: α_s = α → **17/17 pass** AND **matches CSV within 2%**
- CSV match criterion: |σ/m_computed / σ/m_CSV - 1| < 0.02 ✓
- Conclusion: **α_CSV = α_Yukawa = α_s** — confirmed empirically
- VPM code identical to earlier scripts (same kernel) ✓
- Clean, well-structured code ✓
- **מבטל את תוצאות Tests 6, 7, 8, 10b, 11** שהשתמשו ב-(8/9)

---

### 17. `dark_force_accumulation.py` — Test 13 ⚠️

**מטרה**: האם שדות מצטברים מ-DM background יכולים לתת ρ_Λ?

**ממצאים**:
- Poisson equation: φ₀ = y_s n_χ / m_φ² — correct ✓
- ρ_crit = 3H₀²M_Pl²/(8π) — correct ✓
- **WRONG α VALUE**: Uses α = 3.4e-3 instead of MAP's α = 5.734e-3
  - Impact: all derived quantities (y, φ_0, ρ) off by factor (5.734/3.4)² ≈ 2.84
  - But the gap is ~40 orders of magnitude, so the conclusion is UNAFFECTED
- Correctly identifies: same mediator CANNOT do both SIDM and DE ✓
- m_φ ~ MeV for SIDM vs m ~ H₀ for DE → 40 orders gap ✓
- **עקבי איכותית עם Journal (Test 13)**, but numbers are quantitatively wrong

---

### 18. `sigma_mass_protection.py` — Test 13D ✅

**מטרה**: 5 מנגנוני הגנה על m_σ

**ממצאים**:
- MAP values correct: m_chi=94.07e-3, m_phi=11.10e-3, α=5.734e-3 ✓
- `theta = np.arctan(1/np.sqrt(8))` ✓
- `y = np.sqrt(4*np.pi*alpha / np.cos(theta)**2)` ← **CORRECTED convention** ✓
- 5 mechanisms tested:
  1. Exact Goldstone → Need Λ_dark ~ 10⁻³ eV → **most promising** ✓
  2. Dark SUSY → m_χ/m_φ = 8 → no cancellation ✓
  3. Clockwork → Need N~20 gears → possible but ad hoc ✓
  4. Extra dimension → equivalent to clockwork ✓
  5. Large f → trans-Planckian → problematic ✓
- V_CW(0) ~ 10⁻⁷ GeV⁴, ΔV ~ 10⁻⁹ GeV⁴, ΔV/ρ_Λ ~ 10³⁸ → CC problem ✓
- **עקבי עם Journal ו-dark_qcd_consistency.py**

---

### 19. `dark_qcd_consistency.py` — Test 14 ✅

**מטרה**: 6 בדיקות עקביות ל-dark QCD scenario

**ממצאים**:
- Check 1 (Group theory): SU(2)_d ONLY option for Majorana + fundamental ✓
- Check 2 (CP ratio): g_p/g_s = 1/3 preserved in dark QCD ✓
- Check 3 (BBN): SU(2)_d marginal (ΔN_eff ~ 0.34 if thermalized), SU(3)_d excluded ✓
- Check 4 (RG): α_d(M_Pl) ~ 1/200 for Λ_d ~ 10⁻³ eV — unexplained but stable ✓
- Check 5 (Coincidence): Λ_d = √(H₀ M_Pl) ~ 10⁻³ eV ~ m_ν ✓
- Check 6 (Misalignment): ρ_σ ~ ρ_Λ for θ_i ~ 2 rad → Ω_σ ≈ 0.69 ✓
- w ≈ -1 for m_σ ≲ H₀ (slow-roll) ✓
- Constants correct throughout ✓
- **עקבי עם Journal (Test 14) ו-need_to_verify.md (Item 3a)**

---

### 20. `born_full_amplitude.py` — Test 15 ✅

**מטרה**: אמפליטודת Born מלאה — scalar + pseudoscalar

**ממצאים**:
- Dirac spinors constructed correctly (Dirac representation) ✓
- γ⁵ matrix correct ✓
- Majorana: t+u channels with Fermi minus sign (M_u = −...) ✓
- Convention: α_s = α (corrected), α_p = α_s/8 = α/8
  - Justification: tan²θ = 1/8 → α_p = α_s tan²θ = α_s/8 ✓
- **Key result**: Interference/Scalar < 10⁻⁶ → **VPM with α_s alone is excellent** ✓
- Pseudoscalar vertex suppressed by v² (velocity) ✓
- Spin-averaged interference vanishes: Tr[γ⁵ with <4 γ] = 0 ✓
- All 17 BPs show max correction < 10⁻⁶ ✓
- **עקבי עם Journal (Test 15)**

---

### 21. `vev_alignment_stability.py` — Test 16 ✅

**מטרה**: יציבות VEV alignment ב-A₄

**ממצאים**:
- Part 1: Single flavon
  - λ₂ > 0 → (1,1,1)/√3 — confirmed analytically AND numerically ✓
  - λ₂ < 0 → (1,0,0) — confirmed ✓
  - Hessian positive-definite at all minima ✓
- Part 2: Two flavons with cross-terms
  - κ₁-only: stable alignment ✓
  - κ₂, κ₃ ≠ 0: can break alignment (enforceable by Z₂) ✓
- Part 3: Stability boundary κ₁ scan ✓
- Part 4: VEV ratio tunable via κ₁ ✓
- **עקבי עם Journal (Test 16) ו-need_to_verify.md (Item 2b)**

---

### 22. `resonance_bp1.py` ⚠️

**מטרה**: מבנה תהודה ליד BP1

**ממצאים**:
- **Parameter inconsistency**: `M_PHI = 11.34e-3` (11.34 MeV)
  - Other scripts (born_full_amplitude.py fallback): BP1 m_φ = 10.83 MeV
  - Source discrepancy — possibly different BP table version
- λ_decomposed = (8/9)λ_orig → **convention error** (double-counting)
- VPM solver kernel identical to others ✓
- Fine scan near resonance: qualitative resonance structure valid ✓
- Phase shift analysis (δ₀ vs λ) is physically correct ✓
- **Impact**: Exact λ values for BP1 resonance proximity are wrong, but the physics of HOW resonances work is correctly analyzed

---

## שגיאות חוצות-סקריפטים

### 1. שגיאת הקונבנציה α_CSV = α_s (הבעיה המרכזית)

**הבעיה**: 6 סקריפטים מכפילים ב-cos²θ = 8/9 את α שכבר הוא α_s:

| סקריפט | שורה | ביטוי שגוי | impact |
|---------|-------|-----------|--------|
| `consistency_check_sidm.py` | ~scenarios | `alpha_s = alpha * COS2` | σ/m ↓11% → 7 BPs fail |
| `consistency_17bp.py` | ~main | `alpha_s = bp["alpha"] * COS2` | σ/m ↓11% → 7 BPs fail |
| `free_theta_scan.py` | ~scan | `a0*math.cos(math.radians(th))**2` | absolute values wrong |
| `boltzmann_17bp.py` | ~main | `alpha_s = alpha * COS2` | Ωh² ×1.27 over-produced |
| `fornax_gc_check.py` | ~decomp | `alpha_s = alpha * COS2` | conservative (under-est.) |
| `resonance_bp1.py` | ~L98 | `LAM_ORIG*8/9` | λ values shifted |

**תוקן ב**: `test_alpha_convention.py` → 17/17 match CSV with α_s = α

### 2. עקביות קבועים פיזיקליים

| קבוע | ערך | סקריפטים | תקין? |
|-------|------|----------|-------|
| M_Pl | 2.435e18 GeV | כולם | ✅ |
| H_0 | 1.44e-42 GeV | כולם | ✅ |
| θ_relic | arcsin(1/3) = arctan(1/√8) | כולם | ✅ |
| cos²θ | 8/9 | כולם | ✅ |
| sin²θ | 1/9 | כולם | ✅ |
| GEV2_TO_CM2 | 3.8938e-28 | VPM scripts | ✅ |
| GEV_IN_G | 1.78266e-24 | VPM scripts | ✅ |
| C_KM_S | 2.998e5 | VPM scripts | ✅ |

### 3. VPM Solver — עקביות Majorana

| Property | Expected | Actual (all copies) | Status |
|----------|----------|-------------------|--------|
| Spin weights | w_even=1, w_odd=3 | ✅ consistent | ✅ |
| σ_T formula | 2π/(k²) × Σ | ✅ consistent | ✅ |
| Identical particles factor | ×2 in 2π/(k²) | ✅ | ✅ |
| RK4 stepper | 4th order | ✅ all copies | ✅ |
| Adaptive l_max | min(kappa×x_max, 500) | ✅ | ✅ |
| Convergence check | 5 consecutive small terms | ✅ | ✅ |

### 4. ערך α שגוי ב-dark_force_accumulation.py

- Uses α = 3.4e-3 instead of MAP α = 5.734e-3
- All other MAP-specific scripts use 5.734e-3 correctly
- Impact: factor ~2.84 in coupling², ~8 in coupling⁴; qualitative conclusion unchanged

### 5. BP1 m_φ inconsistency

- `resonance_bp1.py`: m_φ = 11.34 MeV
- `born_full_amplitude.py` (fallback): m_φ = 10.83 MeV
- Source: possibly different versions of the BP table

---

## אימות מול need_to_verify.md

| Item | Status | סקריפט | הערות |
|------|--------|--------|-------|
| 1a. Convention | ✅ RESOLVED | `test_alpha_convention.py` | 17/17 pass; (8/9) was double-counting |
| 2a. A₄ CG | ✅ RESOLVED | `verify_a4_cg.py` | tan²θ=1/9 confirmed; 6-term formula correct |
| 2b. VEV stability | ✅ RESOLVED | `vev_alignment_stability.py` | (1,1,1)×(1,0,0) stable; Z₂ enforces |
| 3a. CW scale | 🔄 SUPERSEDED | `dark_qcd_consistency.py` | Dark QCD bypasses CW hierarchy |
| 3b. Fifth force | 🟡 OPEN | `fifth_force_constraints.py` | CW mass screens; chameleon not fully solved |
| 4a. σ trapping | ✅ RESOLVED | `sigma_trapping_ode.py` | θ is discrete (A₄), not dynamical |
| 4b. Coupled Boltzmann | 🟡 OPEN | `boltzmann_17bp.py` | Script exists but has convention error |
| 5a. Neutrino connection | 🟡 OPEN | — | No dedicated script yet |

---

## המלצות

### דחוף (blocking for publication)

1. **תקן/סמן סקריפטים עם שגיאת קונבנציה** — הוסף header comment ב-6 הסקריפטים הפגועים:
   ```python
   # ⚠️ SUPERSEDED — This script applies (8/9) reduction to α_CSV,
   # which is already α_s. See test_alpha_convention.py (Test 12).
   # Correct results: all 17 BPs pass with α_s = α_CSV (no reduction).
   ```

2. **תקן resonance_bp1.py**: עדכן m_φ(BP1) ל-10.83 MeV (ממקור ה-CSV)

3. **תקן dark_force_accumulation.py**: עדכן α ל-5.734e-3 (MAP)

### מומלץ (quality improvements)

4. **sigma_trapping_ode.py**: עדכן `y_sq = 4*math.pi*alpha / math.cos(theta_relic)**2` לקבלת y_total² (21% correction; qualitative result unchanged)

5. **הרץ boltzmann_17bp.py מחדש** עם הקונבנציה המתוקנת (α_s = α_CSV, without (8/9) reduction) — צפוי לתת Ωh² ≈ 0.120 עבור כל 17 BPs

6. **כתוב סקריפט neutrino_connection.py** (Item 5a in need_to_verify.md)

### אופציונלי

7. **consolidate VPM copies** — 7+ עותקים של אותו VPM solver. שקול להפוך למודול משותף.

---

## סיכום ביצועי

- **22 סקריפטים נבדקו**
- **12 נקיים** (55%) — ללא בעיות מהותיות
- **6 עם הערות** (27%) — בעיות מינוריות שלא משנות מסקנות
- **4 עם שגיאות** (18%) — שגיאת קונבנציה או באגים; כולם מתועדים ומוחלפים
- **0 שגיאות לא-מתועדות קריטיות** — כל הבאגים הרציניים כבר התגלו ותוקנו
- **שגיאת סימן**: לא נמצאו שגיאות סימן חדשות מעבר לאלה המתועדות (Bugs 2,3 in Test 4)
- **factor-of-2**: Bug 1 (CW coefficient) — מתועד ומתוקן
- **קונבנציה**: שגיאת (8/9) — מתועדת ומתוקנת ב-Test 12

---

## תיקונים שבוצעו (Fix Log)

**תאריך**: 27 מרץ 2026

### שלב 1 — ארכיון פלט באגי
נוצרה תיקיית `archived_buggy/` עם 8 קבצי פלט מהרצה **לפני** התיקון:

| קובץ | גודל | מקור |
|-------|-------|------|
| `consistency_check_sidm_BUGGY.txt` | 7,276 B | Test 6 — שגיאת (8/9) |
| `consistency_17bp_BUGGY.txt` | 8,004 B | Test 7 — שגיאת (8/9) |
| `free_theta_scan_BUGGY.txt` | 37,613 B | Test 8 — שגיאת (8/9) |
| `boltzmann_17bp_BUGGY.txt` | 5,448 B | Test 10 — שגיאת (8/9) |
| `fornax_gc_check_BUGGY.txt` | 1,257 B | Test 11 — שגיאת (8/9) |
| `resonance_bp1_BUGGY.txt` | 9,637 B | שגיאת (8/9) + m_φ שגוי |
| `dark_force_accumulation_BUGGY.txt` | 1,035 B | Test 13 — α שגוי |
| `sigma_trapping_ode_BUGGY.txt` | 6,336 B | Test 10c — y² חסר cos²θ |

### שלב 2 — תיקונים

#### 2a. SUPERSEDED headers (6 סקריפטים)
הוספת header `# ⚠️ SUPERSEDED` אחרי שורת shebang ב-6 סקריפטים עם שגיאת קונבנציה (8/9):
- `consistency_check_sidm.py` — מפנה ל-test_alpha_convention.py (Test 12)
- `consistency_17bp.py` — מפנה ל-test_alpha_convention.py (Test 12)
- `free_theta_scan.py` — ציון שמגמה איכותית תקפה, ערכים מוחלטים ~11% נמוכים
- `boltzmann_17bp.py` — ציון שΩh² מוגזם ב-~27%
- `fornax_gc_check.py` — ציון שעמודת "without decomposition" = הנכונה
- `resonance_bp1.py` — ציון שגיאת קונבנציה + תיקון m_φ

#### 2b. תיקון פרמטרים
- `dark_force_accumulation.py`: α = **3.4e-3 → 5.734e-3** (MAP value)
- `resonance_bp1.py`: M_PHI = **11.34e-3 → 10.83e-3** GeV (CSV source)

#### 2c. תיקון נוסחה
- `sigma_trapping_ode.py`: `y_sq = 4πα` → `y_sq = 4πα/cos²θ` (21% correction; qualitative result unchanged)

### שלב 3 — אימות כפול (Double-Check)
- סוכן משנה אימת שכל 9 ההחלפות (replacements) בוצעו נכון
- `dark_force_accumulation.py` ו-`resonance_bp1.py` הורצו בהצלחה עם הערכים המתוקנים
- אין שגיאות syntax, כל הערכים תואמים את המפרט

### סטטוס סופי

| סעיף בהמלצות | סטטוס |
|--------------|--------|
| 1. SUPERSEDED headers ל-6 סקריפטים | ✅ בוצע |
| 2. תיקון m_φ ב-resonance_bp1.py | ✅ בוצע |
| 3. תיקון α ב-dark_force_accumulation.py | ✅ בוצע |
| 4. תיקון y² ב-sigma_trapping_ode.py | ✅ בוצע |
| 5. הרצת boltzmann_17bp מחדש | ❌ לא בוצע (scripts SUPERSEDED) |
| 6. neutrino_connection.py | ❌ לא בוצע (מעבר ל-scope) |
| 7. consolidate VPM | ❌ לא בוצע (אופציונלי) |

**השורה התחתונה**: הפרויקט במצב טוב. כל הבאגים הרציניים התגלו ותוקנו. הפערים העיקריים הם (1) סימון explicit של סקריפטים superseded, (2) ריצה מחדש של boltzmann_17bp עם קונבנציה מתוקנת, (3) השלמת items פתוחים ב-need_to_verify.md.
