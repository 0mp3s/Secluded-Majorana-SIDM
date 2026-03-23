# TODO — Next Steps

משימות פתוחות לאחר השלמת ה-consistency checks (predictions/).

| # | משימה | סוג | Impact | מאמץ | סטטוס |
|---|--------|------|--------|------|--------|
| 1 | **SPARC + baryons rotation curve fit** | ואלידציה | ⭐⭐⭐⭐⭐ | בינוני | 🔴 פתוח |
| 2 | **Fornax core profile** (Jeans solver) | תחזית | ⭐⭐⭐⭐ | בינוני | 🔴 פתוח |
| 3 | **Fermi-LAT dSph bounds** (φφ → dark, secondary γ) | constraint | ⭐⭐⭐ | בינוני | 🔴 פתוח |

## פירוט

### 1. SPARC + Baryons Fit
- Fit עקומות סיבוב מלאות: V²_tot = V²_bar + V²_SIDM
- נתוני SPARC (Lelli+2016): photometry 3.6μm → V_bar(r)
- פרמטר חופשי: Υ_* (mass-to-light ratio)
- ניבוי: Υ_* צריך לצאת 0.2–0.8 M_⊙/L_⊙ (פיזיקלי)
- אם עובד → פותר Diversity Problem (Oman+2015)

### 2. Fornax Core Profile
- לפתור Jeans equation עם פוטנציאל SIDM
- לנבא r_core [pc] ופרופיל ρ(r) מדויק
- בדיק מול GAIA DR3 + stellar kinematics
- תחזית אמיתית (blind prediction)

### 3. Fermi-LAT dSph Bounds
- χχ → φφ (secluded annihilation)
- φ → dark sector products (no SM tree-level)
- Secondary γ from bremsstrahlung / inverse Compton
- בדיקת upper limits מול Fermi-LAT dSph stacking

## הושלם

- [x] Gravothermal collapse timescales (predictions/gravothermal/)
- [x] SPARC rotation curve diversity — DM only (predictions/rotation_curves/)
- [x] Cluster merger offsets (predictions/cluster_offsets/)
- [x] ΔN_eff from light mediator (predictions/delta_neff/)
- [x] MCMC posterior sampling (stats_mcmc/)
- [x] χ² fit to 13 astrophysical systems (observations/)
- [x] Relic density — 17 viable BPs (relic_density/)
- [x] Born validation, literature cross-checks (cross_checks/)
- [x] Sommerfeld enhancement (cross_checks/)
- [x] Velocity-averaged cross sections (cross_checks/)
- [x] BBN / mediator cosmology (cosmology/)
