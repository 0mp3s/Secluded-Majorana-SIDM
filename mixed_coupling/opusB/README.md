# opusB — Mixed Coupling Verification Files

**Agent B (Claude Opus 4.6), 23 March 2026**

## Files

| File | Condition | Type | Status |
|------|-----------|------|--------|
| `amplitude_derivation.md` | 1: s-wave ≠ 0 | Analytic (MD) | ✅ Done |
| `nr_potential_derivation.md` | 4: VPM unchanged | Analytic (MD) | ✅ Done |
| `nr_potential_check.py` | 4: VPM unchanged | Numerical (Python) | ✅ Ready to run |

## Quick Run

```bash
cd mixed_coupling/opusB
python nr_potential_check.py
```

## Summary of Conditions

1. **Amplitude** — $a_0 = \pi\alpha_s\alpha_p/(2m_\chi^2) \neq 0$ ✅
2. **4D scan** — Pending (needs coupled Boltzmann first)
3. **Coupled Boltzmann** — Pending
4. **NR potential** — Analytic: corrections $\sim (m_\phi/m_\chi)^2 \sim 10^{-7}$. Numerical check ready.
