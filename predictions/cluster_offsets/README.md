# Cluster-Merger Offset Prediction

## Physics

During a galaxy-cluster merger the dark-matter halo of each sub-cluster
passes through the other.  Self-interactions create a drag force

$$a_{\rm drag} \sim -\frac{\sigma}{m}\,\rho_{\rm DM}\,v^2$$

that decelerates the DM relative to the collisionless galaxies, producing
a measurable offset between the DM (lensing) and galaxy (light) centroids.

Harvey+2015 stacked 72 merging clusters and derived

$$\frac{\sigma}{m} < 0.47\;\text{cm}^2\,\text{g}^{-1} \quad (95\%\;\text{CL})$$

at typical infall velocities $v \sim 1000$–$3000\;\text{km/s}$.

## Data

| File | Content | Source |
|------|---------|--------|
| `clusters_data.csv` | 7 individual systems + Harvey+15 combined limit | Markevitch+2004, Harvey+2015, Massey+2015 |

## Expected Outcome

Our velocity-dependent cross section has $\sigma/m \propto v^{-4}$ in the
Born regime (high $v$).  At cluster velocities the predicted $\sigma/m$ is
orders of magnitude below the Harvey+15 upper bound:

| Benchmark | σ/m at 3000 km/s | Harvey+15 bound | Status |
|-----------|-------------------|-----------------|--------|
| BP1 (λ~4) | ~10⁻⁴ cm²/g | 0.47 cm²/g | ✓ |
| MAP (λ~49) | ~10⁻³ cm²/g | 0.47 cm²/g | ✓ |

This is a **generic success** of Yukawa-mediated models: the cross section
is strongly velocity-suppressed, naturally satisfying cluster bounds while
giving large $\sigma/m$ at dwarf-galaxy velocities.

## Usage

```bash
python predict_offsets.py
```

Output → `output/cluster_sigma_m.png`

## References

- Harvey+2015 — "The non-gravitational interactions of dark matter in colliding galaxy clusters" (arXiv:1503.07675)
- Markevitch+2004 — "Direct Constraints on the Dark Matter Self-Interaction Cross Section from the Merging Galaxy Cluster 1E 0657-56"
- Massey+2015 — "The behaviour of dark matter associated with four bright cluster galaxies in Abell 3827"
- Kahlhoefer+2014 — "Colliding clusters and dark matter self-interactions" (arXiv:1308.3419)
