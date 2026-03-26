# opusB — Mixed Coupling Verification Files

**Agent B (Claude Opus 4.6), 23 March 2026**

## Files

| File | Condition | Type | Status |
|------|-----------|------|--------|
| `amplitude_derivation.md` | 1: s-wave ≠ 0 | Analytic (MD) | ✅ Done |
| `nr_potential_derivation.md` | 4: VPM unchanged | Analytic (MD) | ✅ Done |
| `nr_potential_check.py` | 4: VPM unchanged | Numerical (Python) | ✅ Ready to run |
| `validation_opusA_condition1.md` | 1: validate A | Review (MD) | ✅ Done |
| `validation_opusA_condition2.md` | 2: validate A | Review (MD) | ✅ Done |
| `validation_opusA_condition3.md` | 3: validate A | Review (MD) | ✅ Done |

## Quick Run

```bash
cd mixed_coupling/opusB
python nr_potential_check.py
```

## Summary of Conditions

1. **Amplitude** — $a_0 = 2\pi\alpha_s\alpha_p/m_\chi^2 \neq 0$ ✅ (corrected from B's original factor-4 error)
2. **2D scan** — A's band analysis validated: PASS. Band 0.61 decades, $\alpha_s/\alpha_p \in [13, 212]$.
3. **Coupled Boltzmann / Cannibal** — A's scan validated; CRITICAL bug found in ⟨σv²⟩ formula (×500 overestimate). Qualitative conclusion survives.
4. **NR potential** — Analytic: corrections $\sim (m_\phi/m_\chi)^2 \sim 10^{-7}$. Numerical check ready.
