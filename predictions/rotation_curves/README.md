# Rotation-Curve Diversity Prediction

## Physics

Self-interacting dark matter (SIDM) thermalizes the inner halo within a
radius $r_1$ where the scattering optical depth reaches unity:

$$\rho_{\rm NFW}(r_1)\,\frac{\sigma}{m}\,\sigma_v\,t_{\rm age} \simeq 1$$

Inside $r_1$ the halo develops an isothermal core with approximately constant
density $\rho_0 \approx \rho_{\rm NFW}(r_1)$.  This **lowers** the circular
velocity at small radii relative to NFW.

The "diversity problem" (Oman+2015): at fixed $V_{\rm max}$, observed galaxies
show a wide range of $V(2\,\text{kpc})$ that CDM/NFW cannot reproduce.
Velocity-dependent SIDM can produce diversity because $\sigma/m$ varies with
the characteristic halo velocity.

## Data

| File | Content | Source |
|------|---------|--------|
| `sparc_subset.csv` | 12 galaxies with $V_{\rm max}$, $V(2\,\text{kpc})$, core radii | Oh+2015, de Blok+2008, Adams+2014 |

Galaxies span three categories:
- **Low surface brightness** (DDO 154, DDO 87, IC 2574, NGC 2366, UGC 128, UGC 5750) — large cores, low $V(2\,\text{kpc})$
- **Spirals** (NGC 925, NGC 2403, NGC 2976, NGC 3198, UGC 2259) — intermediate
- **High surface brightness** (NGC 7814) — small core, high $V(2\,\text{kpc})$

## Expected Outcome

| Benchmark | λ | Expected behaviour |
|-----------|---|-------------------|
| BP1 | ~4 | Modest cores; good match to observed diversity |
| BP16 | ~4 | Similar to BP1 |
| MAP | ~333 | Very large σ/m at low v → large cores → may over-predict diversity for low-mass halos |

## Usage

```bash
python predict_core_sizes.py              # default config.json
python predict_core_sizes.py --config alt_config.json
```

Output → `output/rotation_curve_diversity.png`

## References

- Oman+2015 — "The unexpected diversity of dwarf galaxy rotation curves" (arXiv:1504.01437)
- Kaplinghat+2016 — "Dark Matter Halos as Particle Colliders" (arXiv:1508.03339)
- Kamada+2017 — "Self-Interacting Dark Matter Can Explain Diverse Galactic Rotation Curves" (arXiv:1611.02716)
- Lelli+2016 — "SPARC: Mass Models for 175 Disk Galaxies" (arXiv:1606.09251)
