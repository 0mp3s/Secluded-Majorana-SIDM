# Need to Verify — Dark EM Duality + A₄ Framework

## Summary

This document lists all claims that require independent verification before this framework can be considered publishable. Each item specifies **what** to check, **where** (which project/directory), **how** (method), and **priority**.

---

## 1. SIDM Consistency Under θ-Decomposition

### 1a. VPM σ/m with α → α_s — Convention Correction

| | |
|---|---|
| **Claim** | ~~10/17 BPs survive σ/m ≥ 0.5 with cos²θ = 8/9 reduction~~ → **17/17 pass** |
| **Current status** | ✅ VERIFIED — see `test_alpha_convention.py` |
| **Result** | **Convention error found** (26 Mar 2026, cross-agent discussion): α_CSV in `sweep_17bp_results.csv` is already α_Yukawa = α_s, NOT α_total. The VPM solver (`core/v22_raw_scan.py` line 16) uses α directly as the Yukawa coupling. Therefore the (8/9) reduction was double-counting. Empirically confirmed: with α_s = α_CSV (no reduction), **17/17 BPs pass σ/m(30) ≥ 0.5**. The earlier 10/17 count was wrong. A₄ is SIDM-degenerate: it does not change any SIDM observable. |
| **Priority** | ✅ RESOLVED (corrected) |

### 1b. σ/m(30) Threshold Sensitivity

| | |
|---|---|
| **Claim** | ~~All 17/17 pass with 0.45 cut; 10/17 with 0.5 cut~~ → **17/17 pass with 0.5 cut** (convention corrected) |
| **What to verify** | Literature scan — what is the actual observational lower bound on σ/m? |
| **Where** | Literature review (not code) |
| **How** | Check Elbert+2015, KTY16, Kaplinghat+2016 for exact dwarf σ/m constraints |
| **Priority** | 🟡 MEDIUM — affects interpretation, not calculation |

---

## 2. A₄ Clebsch-Gordan Coefficients

### 2a. 3⊗3⊗3 → 1 Contraction

| | |
|---|---|
| **Claim** | (ψ̄ψξ_s)₁ = 3, (ψ̄ψξ_p)₁ = 1, giving tan²θ = 1/9 |
| **Current status** | ✅ VERIFIED — see `verify_a4_cg.py` |
| **Result** | 9-term formula fails S-invariance for general vectors (bug in paper). Correct 6-term contraction gives g_s=2, g_p=2/3. BUT: ratio g_p/g_s = 1/3 is IDENTICAL → tan²θ = 1/9 CONFIRMED. SymPy: tan²θ = 1/9 exact. Majorana symmetry OK. sin²θ = 1/10 for equal VEVs; 6% VEV correction (v_p/v_s=1.061) gives exact 1/9. |
| **Priority** | ✅ RESOLVED |

### 2b. VEV Alignment Stability

| | |
|---|---|
| **Claim** | ⟨ξ_s⟩ ∝ (1,1,1) and ⟨ξ_p⟩ ∝ (1,0,0) are stable minima of A₄-invariant potential |
| **Current status** | ✅ VERIFIED — see `vev_alignment_stability.py` |
| **Result** | Single flavon: λ₂ > 0 → (1,1,1), λ₂ < 0 → (1,0,0), all Hessians positive. Two flavons: alignment stable for κ₁-only cross-coupling (21/21 points). κ₂,κ₃ ≠ 0 breaks alignment — enforceable by Z₂ symmetry (ξ_p → -ξ_p). |
| **Priority** | ✅ RESOLVED |

### 2c. sin²θ = 1/10 vs 1/9 — Is the Gap Physical?

| | |
|---|---|
| **Claim** | A₄ with equal VEVs gives sin²θ = 1/10; we need 1/9; gap closed by v_p/v_s = 1.061 |
| **What to verify** | Can higher-order operators in the flavon potential naturally give v_p/v_s ≈ 1.06? |
| **Where** | `dark-energy-T-breaking/` — analytic calculation |
| **How** | 1. Add dim-6 operators to the two-flavon potential <br> 2. Check if radiative corrections shift VEV ratio <br> 3. Compute the fine-tuning measure (how sensitive is v_p/v_s to UV parameters?) |
| **Priority** | 🟡 MEDIUM — 6% deviation could be natural or could require tuning |

---

## 3. Dark Axion as Dark Energy

### 3a. Coleman-Weinberg Potential Scale

| | |
|---|---|
| **Claim** | V_CW for dark axion σ can contribute to ρ_Λ with appropriate f |
| **Current status** | ⚠️ CW alone insufficient — BUT dark QCD scenario bypasses this entirely |
| **Update (26 Mar 2026)** | Dark QCD consistency check (Test 14) shows σ = dark pion with $m_\sigma = \Lambda_d^2/f \sim H_0$ via misalignment mechanism. CW scale no longer the DE source — σ's frozen misalignment energy gives $\Omega_\sigma = 0.69$ for $\theta_i \sim 2$ rad. See `dark_qcd_consistency.py`. |
| **What to verify** | (a) Full dark QCD lattice estimates for Λ_d; (b) σ potential beyond leading order; (c) CW corrections to dark pion mass |
| **Priority** | 🟡 MEDIUM — dark QCD resolves the CW hierarchy, but new checks needed for confinement sector |

### 3b. Fifth Force Constraints

| | |
|---|---|
| **Claim** | For f ≈ 0.2 M_Pl, fifth force bounds are satisfied via chameleon/symmetron screening |
| **Current status** | Estimated in `fifth_force_constraints.py`; screening assumed, not proven |
| **What to verify** | Solve the chameleon field equation in realistic density profiles |
| **Where** | `dark-energy-T-breaking/` — new script |
| **How** | 1. Implement chameleon screening for our specific σ potential <br> 2. Check Eöt-Wash, MICROSCOPE, Cassini bounds <br> 3. Verify that β < 10⁻² in solar system |
| **Priority** | 🟡 MEDIUM — needed for paper but generic screening usually works at f ~ M_Pl |

---

## 4. Freeze-Out Dynamics

### 4a. σ Trapping at θ_relic During Freeze-Out

| | |
|---|---|
| **Claim** | The dark axion σ gets trapped at θ_relic = arcsin(1/3) during freeze-out |
| **Current status** | ✅ VERIFIED (negative result) — see `sigma_trapping_ode.py` |
| **Result** | CW dominates thermal by 10¹⁰. σ rolls to π/2 for ALL initial conditions and ALL f values. σ trapping as dynamical mechanism DOES NOT WORK. **Resolution**: θ is NOT a dynamical field — it's a GROUP THEORY CONSTANT from A₄ CG decomposition. The ratio g_p/g_s = 1/3 is protected by discrete A₄ symmetry. No σ trapping needed → no trapping problem. |
| **Priority** | ✅ RESOLVED (reframed: θ is not dynamical) |

### 4b. Coupled Boltzmann + σ Evolution

| | |
|---|---|
| **Claim** | The σ-dependent coupling doesn't destabilize the standard freeze-out picture |
| **What to verify** | Full numerical Boltzmann with time-dependent α_s(σ(T)) and α_p(σ(T)) |
| **Where** | `Secluded-Majorana-SIDM/` — modify Boltzmann solver to accept α(T) |
| **How** | 1. Feed σ(T) trajectory from freeze_out_trapping into Boltzmann solver <br> 2. Compare Ωh² with fixed-coupling result <br> 3. Estimate correction to relic density |
| **Priority** | 🟡 MEDIUM — probably small effect (σ changes slowly vs H) |

---

## 5. Neutrino Connection

### 5a. A₄ Model Consistency with Neutrino Data

| | |
|---|---|
| **Claim** | Same A₄ group gives sin²θ₁₂ = 1/3 for neutrinos and sinθ = 1/3 for dark sector |
| **What to verify** | Build complete A₄ model with BOTH sectors; check no conflicts in quantum numbers |
| **Where** | `dark-energy-T-breaking/` — new theory document |
| **How** | 1. List all fields with A₄ × U(1)_D × SM quantum numbers <br> 2. Check anomaly cancellation <br> 3. Verify no unwanted couplings between dark and visible sectors <br> 4. Compare with Altarelli-Feruglio SUSY SU(5)×A₄ (arXiv:0802.0090) |
| **Priority** | 🟡 MEDIUM — important for paper motivation but doesn't affect dark sector phenomenology |

### 5b. TBM Corrections (θ₁₃ ≠ 0)

| | |
|---|---|
| **Claim** | The 1/3 connection persists despite TBM being approximate (θ₁₃ ≈ 8.5°) |
| **What to verify** | Does A₄ breaking that gives θ₁₃ ≠ 0 also shift θ_dark from arcsin(1/3)? |
| **Where** | `dark-energy-T-breaking/` — analytic |
| **How** | Estimate: if A₄ → A₄+corrections, what is δθ_dark/δθ₁₃? |
| **Priority** | 🟢 LOW — second-order effect, but important for completeness |

---

## 6. Reproducibility of SIDM Results

### 6a. VPM Solver Validation

| | |
|---|---|
| **Claim** | VPM solver copied from SIDM project gives correct σ_T |
| **Current status** | Matches SIDM benchmarks (Test 6 in journal) |
| **What to verify** | Compare with independent VPM implementation or micrOMEGAs |
| **Where** | `Secluded-Majorana-SIDM/core/v22_raw_scan.py` (original) |
| **How** | 1. Pick 3 BPs, compute σ_T with published Tulin/Yu code or analytic limits <br> 2. Check Born limit (λ≪1) and classical limit (λ≫1) analytically <br> 3. Verify Majorana spin weights (w_even=1, w_odd=3) |
| **Priority** | 🟢 LOW — already validated in SIDM project, but good to document |

---

## 7. Dark Sector Production Mechanism

### 7a. Test 21 — FIMP Freeze-In Production

| | |
|---|---|
| **Claim** | Dark sector produced via freeze-in (FIMP) with coupling λ_hs ~ 5×10⁻⁴, giving T_D = 200 MeV as effective initial temperature |
| **What to verify** | Compute Ω_χ h² via freeze-in integral: $\Omega_{FIMP} \propto \int dT \, \frac{g_{SM} T}{H} \Gamma_{prod}(T)$; check if result is consistent with Ω_DM h² = 0.120 |
| **Where** | `dark-energy-T-breaking/` — new script `test21_fimp_production.py` |
| **How** | 1. Integrate production rate Γ(T) from T_RH down to T_D <br> 2. Compare to observed relic density <br> 3. If Ω_FIMP > Ω_DM → overproduction, need smaller coupling <br> 4. If Ω_FIMP << Ω_DM → T_D is NOT the decoupling temperature, something else sets the abundance |
| **Priority** | 🔴 HIGH — T_D = 200 MeV remains an assumption until this is verified |

### 7b. Test 22 — ΔN_eff Including φ → 2σ

| | |
|---|---|
| **Claim** | ~~ΔN_eff = 0.153 + 0.027 ≈ 0.180~~ **תוקן**: ΔN_eff ≈ 0 ב-secluded model (χ,φ הם NR ב-BBN). אם σ מתרמלז: 0.214. |
| **Result** | ✅ **DONE** (28 Mar 2026). הנוסחה הישנה (massless limit) **שגויה**: ספרה d.o.f. ב-$T_D$ אבל ΔN_eff נמדד ב-BBN. Check 4 תוקן עם Boltzmann suppression. |
| **Where** | `test22_neff_phi_decay.py` + `test_physics_checks.py` (Check 4) |
| **Priority** | ✅ **RESOLVED** |

---

## Verification Priority Summary

| Priority | Item | Where | Blocking? |
|---|---|---|---|
| ✅ DONE | Full Boltzmann relic + convention correction | `test_alpha_convention.py` | Resolved — α_CSV = α_Yukawa; (8/9) was wrong; 17/17 pass |
| ✅ DONE | A₄ CG cross-check | `verify_a4_cg.py` | Resolved — tan²θ=1/9 confirmed |
| ✅ DONE | σ trapping during freeze-out | `sigma_trapping_ode.py` | Resolved — θ is discrete, not dynamical |
| ✅ DONE | **PI-2: θ_A₄ instanton stability** | `SIDM-LAGANJIAN_INTEGRAL/test_PI2_instanton_stability.py` | Resolved — $S_E\sim10^{121}$, bubble wall = Hubble radius; ABSOLUTELY STABLE |
| � MED | CW potential scale | `dark-energy-T-breaking/` | Dark QCD bypasses CW hierarchy; new checks needed for confinement sector |
| 🟡 MEDIUM | σ/m threshold literature | Literature review | No |
| ✅ DONE | VEV alignment stability | `vev_alignment_stability.py` | Resolved — (1,1,1)×(1,0,0) stable with κ₁-only cross-coupling |
| 🟡 MEDIUM | sin²θ gap (1/10 vs 1/9) | `dark-energy-T-breaking/` | No |
| 🟡 MEDIUM | Fifth force screening | `dark-energy-T-breaking/` | No |
| 🟡 MEDIUM | Coupled Boltzmann + σ(T) | `Secluded-Majorana-SIDM/` | No |
| 🟡 MEDIUM | A₄ + neutrino consistency | `dark-energy-T-breaking/` | No |
| � HIGH | **Test 21: FIMP production** | `dark-energy-T-breaking/` | **Yes — T_D=200 MeV assumption** |
| ✅ DONE | **Test 22: ΔN_eff with φ→2σ** | `dark-energy-T-breaking/` | **Resolved — ΔN_eff ≈ 0 (secluded), 0.214 (if σ thermalizes)** |
| �🟢 LOW | TBM corrections to θ_dark | `dark-energy-T-breaking/` | No |
| 🟢 LOW | VPM solver validation | `Secluded-Majorana-SIDM/` | No |
