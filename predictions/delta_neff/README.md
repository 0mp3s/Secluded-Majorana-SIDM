# Extra Radiation ΔN_eff Prediction

## Physics

The light scalar mediator $\phi$ forms a separate thermal bath that decouples
from the Standard Model at $T_{\rm dec} \sim m_\chi$.  After decoupling, SM
entropy transfers reheat the SM relative to the dark sector:

$$\frac{T_\phi}{T_{\rm SM}} = \left(\frac{g_{*S}(T_{\rm SM})}{g_{*S}(T_{\rm dec})}\right)^{1/3}$$

A massless $\phi$ would contribute

$$\Delta N_{\rm eff}^{\rm massless} = \frac{4}{7}\left(\frac{T_\phi}{T_\nu}\right)^4 \approx 0.036$$

However, in our model $m_\phi \sim 10\text{–}14\;\text{MeV}$, so at the BBN epoch
($T_{\rm SM} \sim 1\;\text{MeV}$, $T_\phi \sim 0.5\;\text{MeV}$) the mediator
is non-relativistic with $m_\phi / T_\phi \sim 20$.  The energy density is
Boltzmann-suppressed:

$$\rho_\phi \propto e^{-m_\phi / T_\phi} \sim e^{-20} \approx 10^{-9}$$

**Result:** $\Delta N_{\rm eff} \lesssim 10^{-8}$ at BBN and CMB epochs —
automatically consistent with all cosmological bounds.

## Data Source

- **Planck 2018:** $N_{\rm eff} = 2.99 \pm 0.17$ (arXiv:1807.06209)
- **SM prediction:** $N_{\rm eff} = 3.044$
- **CMB-S4 forecast:** $\sigma(N_{\rm eff}) = 0.03$ (arXiv:1610.02743)

## Expected Outcome

| Parameter | Value |
|-----------|-------|
| $T_\phi / T_\nu$ | ~0.50 |
| $\Delta N_{\rm eff}$ (massless limit) | 0.036 |
| $\Delta N_{\rm eff}$ at BBN | $\lesssim 10^{-8}$ |
| $\Delta N_{\rm eff}$ at CMB | $\sim 0$ |
| Planck 2018 status | ✓ SAFE |
| CMB-S4 status | ✓ SAFE |

The Boltzmann suppression from $m_\phi \gg T_\phi$ makes this prediction
trivially consistent.  This is an important **null prediction**: it guarantees
that the light mediator does not spoil BBN or CMB constraints.

## Usage

```bash
python predict_neff.py
```

Output → `output/delta_neff.png`

## References

- Planck Collaboration (2018) — arXiv:1807.06209
- CMB-S4 Science Case (2016) — arXiv:1610.02743
- Kolb & Turner — *The Early Universe*, Ch. 3
