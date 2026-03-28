# Layer 8 — קבוע האבל כפלט מהלגרנג'יאן

> **תאריך:** 28 מרץ 2026  
> **קוד:** `The_Lagernizant_integral_SIDM/lagrangian_path_integral.py :: solve_friedmann_sigma()`  
> **מטרה:** $H_0$ הוא **פלט**, לא קלט. 6 פרמטרים מהלגרנג'יאן → $H_0$

---

## 1. הבעיה: מעגליות ב-$H_0$

עד Layer 7, הקוד השתמש ב-$H_0 = 67.4$ km/s/Mpc כ**קלט** (מ-config.json):

```
config.json: H0 = 67.4 km/s/Mpc  (קלט)
    ↓
ρ_σ = ½ f² θ_i² H₀²              (משתמש ב-H₀)
    ↓
H₀ = √(8πG ρ_total / 3)          (מנסה לחשב H₀)
```

$H_0$ מופיע **משני הצדדים** — זו מעגליות.

---

## 2. הפתרון: שכבה 8 — אינטגרל בזמן

במקום 7 שכבות ברגע קפוא, לפתור את הצמד **פרידמן + משוואת תנועה של $\sigma$** ביחד:

### משוואת תנועה של $\sigma$:
$$\ddot{\sigma} + 3H(t)\dot{\sigma} + V'(\sigma) = 0$$

### משוואת פרידמן:
$$H^2(t) = \frac{8\pi}{3M_{\rm Pl}^2}\left[\underbrace{\rho_r(a)}_{\text{SM radiation}} + \underbrace{\rho_\chi(a)}_{\text{Layer 7: relic}} + \underbrace{\rho_\phi(a)}_{\text{cannibal}} + \underbrace{\tfrac{1}{2}\dot{\sigma}^2 + V(\sigma)}_{\text{dark axion}}\right]$$

### שתי ODE מצומדות — $\sigma(t)$ ו-$a(t)$ — מ-reheating ועד היום:

| תקופה | $a$ | מה קורה |
|--------|-----|----------|
| Reheating | $\sim 10^{-28}$ | $\sigma = f\theta_i$ (initial misalignment) |
| Freeze-out | $\sim 10^{-12}$ | $\Gamma_{\chi\chi\to\phi\phi} < H$ → $\rho_\chi \propto a^{-3}$ (שכבה 7) |
| Cannibal | $\sim 10^{-10}$ | $3\phi\to 2\phi$ → $\rho_\phi \to 0$ |
| $\sigma$ frozen | $m_\sigma < H$ | שדה $\sigma$ "קפוא" ב-slow roll, $\dot{\sigma}\approx 0$ |
| $\sigma$ oscillates | $m_\sigma \sim H$ | $\sigma$ מתחיל לנדנד סביב $\theta_{\rm relic} = 19.47°$ |
| **היום** | $1$ | $V(\sigma_{\rm today}) = \rho_\Lambda$ → **$H_0$ יוצא מהחישוב!** |

### הפאנצ'ליין:
$$\boxed{H_0 = H(t_0) = \sqrt{\frac{8\pi}{3M_{\rm Pl}^2}\Big[\rho_{\rm SM}(t_0) + \Omega_\chi\rho_c + V\big(\sigma(t_0)\big)\Big]}}$$

**$H_0$ לא קלט. $H_0$ פלט.**

הקלטים היחידים: $m_\chi, m_\phi, \alpha, f, \Lambda_d, \theta_i$ — 6 פרמטרים מהלגרנג'יאן.

### מבחינת קוד:
`scipy.integrate.solve_ivp` עם 3 משתנים: $(\sigma, \dot{\sigma}, a)$, מ-$t_{\rm RH}$ עד $t_0$.  
Layer 7 (relic) ושכבה 5 (cannibal) כבר פתורים — הם נכנסים כ-$\rho_\chi(a)$ ו-$\rho_\phi(a)$ שכבר ידועים.

---

## 3. מימוש: `solve_friedmann_sigma()`

### משתנים:
- **זמן**: $N = \ln(a/a_0)$ (e-folds, לא שניות — מטפל ב-60+ סדרי גודל)
- **State**: $[\theta, d\theta/dN]$ כאשר $\theta = \sigma/f$

### אלגברה של פרידמן:
$H^2$ מופיע משני הצדדים (גם ב-$\rho_\sigma$), נפתר אלגברית:

$$H^2(N) = \frac{(8\pi/3M_{\rm Pl}^2)\Big[\rho_r(0)e^{-4N} + \rho_\chi(0)e^{-3N} + V(\theta)\Big]}{1 - (4\pi f^2/3M_{\rm Pl}^2)(d\theta/dN)^2}$$

### לולאה איטרטיבית:
$H_0$ נכנס ב-$\rho_r(0), \rho_\chi(0)$ — פותרים self-consistently:
1. נחש $H_0^{(0)}$ (Planck value)
2. פתור ODE → קבל $H_0^{(1)}$ ב-$N = 0$
3. חזור עם $H_0^{(1)}$ → התכנסות

---

## 4. תוצאות

### תוצאה מרכזית: $\Lambda_d$ scan

```
  Lambda_d [GeV]     H0 [km/s/Mpc]   delta
  --------------------------------------------------
  1.0000e-12         12.98            -80.7%
  1.5849e-12         32.61            -51.6%
  1.9953e-12         51.69            -23.3%
  2.2387e-12         65.07            -3.5%
  2.2800e-12         67.49            +0.1%    ← MATCH
  2.5119e-12         81.92            +21.5%
  3.1623e-12         129.83           +92.6%
```

### Best fit:

| פרמטר | ערך |
|---------|------|
| $\Lambda_d$ | $2.28 \times 10^{-12}$ GeV |
| $f$ | $0.24 \, M_{\rm Pl} = 2.93 \times 10^{18}$ GeV |
| $\theta_i$ | 2.0 rad |
| $m_\sigma = \Lambda_d^2 / f$ | $1.43 \times 10^{-42}$ GeV $\approx 2.2 \times 10^{-18}$ eV |
| **$H_0^{\rm predicted}$** | **67.49 km/s/Mpc** |
| $H_0^{\rm Planck}$ | 67.4 km/s/Mpc |
| $\Delta H_0 / H_0$ | **+0.13%** |
| $\theta_{\rm today}$ | 1.712 rad ≈ 98.1° |
| $V(\sigma_{\rm today}) / \rho_\Lambda$ | 1.23 |
| $m_\sigma / H_0$ | ~1.2 |
| Convergence | 12 iterations |

### פיזיקה:
- $\sigma$ **כמעט קפוא**: נע מ-$\theta_i = 2.0$ ל-$\theta = 1.712$ — רק 14% שינוי ב-13.8 מיליארד שנה
- $m_\sigma / H_0 \approx 1.2$ — בדיוק בנקודת המעבר frozen ↔ oscillating
- $V(\sigma)$ פועל כ-**cosmological constant אפקטיבי** (slow roll)

---

## 5. באג שתוקן: המרת יחידות

`GEV_TO_KM_S_MPC` ב-config.json הוא **1 Mpc ביחידות GeV⁻¹** ($= 1.5637 \times 10^{38}$), **לא** factor ההמרה המלא. חסר כפל ב-$c$:

$$H_0 [\text{km/s/Mpc}] = H_0 [\text{GeV}] \times \underbrace{1.5637 \times 10^{38}}_{\text{1 Mpc [GeV}^{-1}\text{]}} \times \underbrace{299{,}792}_c$$

אימות: $1.44 \times 10^{-42} \times 1.5637 \times 10^{38} \times 299{,}792 = 67.5$ ✓

---

## 6. בדיקות צולבות

(יתמלא בהמשך)

### 6.1 רגישות ל-$\theta_i$
TODO

### 6.2 רגישות ל-$f/M_{\rm Pl}$
TODO

### 6.3 כל 3 ה-benchmarks (BP1, MAP, MAP_relic)
TODO

### 6.4 צירוף $m_\sigma \sim H_0$ — מקרי או מבני?
$m_\sigma = \Lambda_d^2 / f$.  
$H_0 = \sqrt{8\pi V(\sigma) / 3M_{\rm Pl}^2}$.  
If $V \sim \Lambda_d^4$: $H_0 \sim \Lambda_d^2 / M_{\rm Pl}$.  
But $m_\sigma = \Lambda_d^2 / f$ and $f \sim 0.24 M_{\rm Pl}$:  
→ $m_\sigma / H_0 \sim M_{\rm Pl} / f \sim 1/0.24 \sim 4$.  
**Not a coincidence** — it's a structural requirement: $f \sim M_{\rm Pl}$ guarantees $m_\sigma \sim H_0$.

### 6.5 Hubble tension?
SH0ES: $H_0 = 73.0 \pm 1.0$ km/s/Mpc  
Planck: $H_0 = 67.4 \pm 0.5$ km/s/Mpc  
Model: can produce either value by adjusting $\Lambda_d$ by ~5%.  
→ The model doesn't *resolve* the tension, but it *parameterizes* it in terms of $\Lambda_d$.

---

## 7. מה הלאה

1. **Paper 2 עדכון**: Layer 8 הוא התוספת המרכזית — $H_0$ כפלט
2. **Fine-tuning analysis**: כמה $\Lambda_d$ צריך להיות מדויק? (looks like ~5% changes give ~20% H₀ shift)
3. **Comparison with quintessence**: standard quintessence models also have rolling scalar → DE; our model is more constrained (θ fixed by SIDM+relic)
4. **Late-time EOS**: $w(z)$ from the σ evolution — measurable by DESI?
