# 🔑 Hunt for H₀ — Finding the Car Keys

**מטרה**: למצוא את קבוע האבל $H_0$ כ**תוצאה** של הלגרנזיאן של SIDM — לא כקלט.

> "חייבים למצוא. זה כמו מפתחות של אוטו — הם איפשהו."

---

## נקודת המוצא — מה יש לנו

הלגרנזיאן:
$$\mathcal{L} = \frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2 + \bar\chi(i\!\not\!\partial - m_\chi)\chi - (g_s + g_p\gamma^5)\phi\bar\chi\chi$$

פרמטרים: $m_\chi, m_\phi, \alpha_D$ (+ $\theta_{A_4} = \arcsin(1/3)$ מ-A₄)

מה שכבר **נגזר** מזה:
- ✅ SIDM cross section $\sigma_T(v)$ — דרך VPM
- ✅ Relic density $\Omega_\chi h^2 = 0.120$ — דרך Boltzmann
- ✅ זווית CP $\theta = 19.47°$ — דרך A₄
- ✅ VPM = exact (fluctuation determinant) — PI-6
- ❌ $H_0$ — עדיין חסר

---

## מפת כל הכיוונים — 7 נתיבים

### נתיב 1️⃣: Dark QCD + Misalignment (PI-7, כבר נבדק)

**רעיון**: $\sigma$ כ-dark pion (pNGB), $V(\theta) = \Lambda_d^4(1-\cos\theta)$, misalignment → $\rho_\Lambda$

**מה עובד**:
- $\Omega_\sigma = 0.69$ ל-$\theta_i \sim 2$ rad, $f \sim 0.2 M_{Pl}$ ✅
- $m_\sigma \sim H_0$ מ-GMOR: $m_\sigma = \Lambda_d^2/f$ ✅
- $\Lambda_d \sim 2$ meV $= \sqrt{H_0 M_{Pl}}$ — סקאלת הנויטרינו ⚡
- A₄ מקבעת $\theta_{dark}$, tunneling יציב ($S_E \sim 10^{121}$) ✅

**מה לא עובד**:
- $\Lambda_d$ הוא **קלט** — לא נגזר מ-$\alpha_D, m_\chi, m_\phi$ ❌
- $\theta_i \sim 2$ rad פרמטר חופשי (PI-3: stochastic inflation לא נותן) ❌
- **סירקולריות**: $\Lambda_d = \sqrt{H_0 M_{Pl}}$ ← משתמש ב-$H_0$ כקלט ❌

**סטטוס**: ✅ עקבי, אבל **לא** prediction — consistency check

---

### נתיב 2️⃣: Transmutation ממדי — $\Lambda_d$ מ-$\alpha_d(\mu)$

**רעיון**: כמו ש-$\Lambda_{QCD}$ נובע מ-$\alpha_s(M_Z)$ דרך RG running:
$$\Lambda_d = \mu \cdot \exp\left(-\frac{2\pi}{b_0\,\alpha_d(\mu)}\right)$$

**איך זה נותן $H_0$**: אם $\alpha_d$ **קשור** ל-$\alpha_{SIDM}$ שלנו (למשל unification ב-$M_{GUT}$), אז:
$$\alpha_D(m_\phi) = 5.734 \times 10^{-3} \xrightarrow{\text{RG}} \alpha_d(M_{Pl}) \xrightarrow{\text{transmutation}} \Lambda_d \xrightarrow{\text{GMOR}} m_\sigma \xrightarrow{\text{Friedmann}} H_0$$

**מה חסר**:
- הקשר $\alpha_d \leftrightarrow \alpha_D$ (unification? portal? אותו gauge group?)
- $b_0$ (תלוי בתוכן SU(N)_d — Check 3 נתן SU(2)_d, $b_0 = 6.33$)
- Check 4 מצא: $\alpha_d(M_{Pl}) \sim 1/71$ נותן $\Lambda_d \sim 10^{-3}$ eV ← **אבל 1/71 הוא קלט**

**מה נדרש**: בדיקה — האם $\alpha_D = 5.734 \times 10^{-3}$ ו-$\alpha_d(M_{Pl}) = 1/71$ **עקביים** תחת RG של SU(2)_d? אם כן — $H_0$ נגזר.

**סטטוס**: 🔶 לא נבדק — **פוטנציאל גבוה ביותר**

---

### נתיב 3️⃣: אינטגרל זמן קוסמי — Boltzmann → Friedmann מלא

**רעיון** (של עומר): לא אינטגרל שדות (כבר עשינו). אינטגרל על **היסטוריה**:
$$S_{cosmo} = \int_0^{t_0} dt\, a^3(t)\, \mathcal{L}_{eff}(T(t))$$

שני אינטגרלים מקוננים:
```
אינטגרל פנימי (QFT):   ∫Dφ Dχ e^{iS}  →  V_eff(T), ⟨σv⟩, Y(T)
                              ↓
אינטגרל חיצוני (זמן):  ∫Da(t) e^{iS_cosmo}  →  H₀
```

**הנקודת אוכף של האינטגרל החיצוני = משוואת Friedmann**:
$$H^2(t) = \frac{8\pi G}{3}\Big[\rho_{rad}(a) + \rho_b(a) + \rho_\chi(a) + \rho_\Lambda\Big]$$

**מה זה נותן**: שרשרת **סגורה** הלגרנזיאן → $H_0$:
1. $\alpha_D, m_\chi \to Y_\infty \to \rho_{\chi,0}$ (Boltzmann — כבר יש)
2. SM: $\rho_b + \rho_r$ (ידוע)
3. $\rho_\Lambda$ — **חסר** (= הבעיה)
4. Friedmann: $H_0 = \sqrt{8\pi G \cdot \rho_{total}/3}$

**מה חסר**: $\rho_\Lambda$. בלי מקור ל-DE, $H_0$ יוצא 40-50 km/s/Mpc (DM+baryons+radiation בלבד).

**אפשרות רדיקלית**: אולי $\rho_\Lambda$ **מיותר** אם האינטגרל הקוסמי על $V_{eff}(T)$ מצטבר?

**סטטוס**: 🔵 הרעיון של עומר — צריך לבנות

---

### נתיב 4️⃣: CW Vacuum Energy — cancellation דינמי

**רעיון**: $V_{CW} \sim 10^{-7}$ GeV⁴ = $10^{40} \times \rho_\Lambda$. אבל $V_{CW}(T)$ **תלוי בטמפרטורה**.
אולי האינטגרל **על כל** ההיסטוריה (radiation → matter → DE) מייצר cancellation?

$$\rho_\Lambda^{eff} = \int_0^{t_0} dt\, \frac{dV_{CW}(T(t))}{dt}$$

**בעיה**: $V_{CW}$ הוא רגעי (static vacuum energy), לא מצטבר. כדי שזה יצטבר צריך מנגנון **דינמי** — למשל שדה $\sigma$ שמתגלגל ומצבר אנרגיה לאורך הזמן.

**קשר לנתיב 1**: אם $\sigma$ כ-QF field, האנרגיה הקינטית + פוטנציאלית שלו היום = $\rho_\Lambda$

**סטטוס**: 🔴 ספקולטיבי — צריך בדיקת היתכנות

---

### נתיב 5️⃣: Trace Anomaly — $\langle T^\mu_\mu \rangle$ מ-path integral

**רעיון**: ב-QCD, ה-trace anomaly נותן:
$$\langle T^\mu_\mu \rangle_{QCD} = \frac{\beta(g)}{2g} F_{\mu\nu}^a F^{a\mu\nu} + m_q\bar{q}q$$

זה הבסיס ל-95% ממסת הפרוטון (מאנרגיית binding של גלואונים, לא ממסת quark).

**אנלוגיה ל-dark QCD**: 
$$\langle T^\mu_\mu \rangle_{dark} = \frac{\beta_d(\alpha_d)}{2\alpha_d} G_{d,\mu\nu}^a G_d^{a\mu\nu}$$

אם ה-dark gluon condensate שורד עד היום → $\rho_{condensate} \sim \Lambda_d^4$

**מה מיוחד**: ה-trace anomaly **חייב** להופיע מהאינטגרל — הוא פיצ'ר של QFT (1-loop exact). ולכן $\Lambda_d^4$ הוא **תוצאה** של $\alpha_d$, לא פרמטר חופשי.

**חיבור לנתיב 2**: $\Lambda_d^4 = \Lambda_d^4(\alpha_d(M_{Pl}))$ — transmutation ממדי. ה-trace anomaly הוא ה**מנגנון** שבו $\alpha_d$ הופך ל-$\rho_\Lambda$.

**סטטוס**: 🔶 תיאורטי מבטיח — צריך חישוב

---

### נתיב 6️⃣: Friedmann ללא $\Lambda$ — $H_0$ מ-DM+baryons בלבד

**רעיון**: בדיקת בסיס. מה $H_0$ מתקבל אם **אין DE כלל**?

$$H_0^{(no\Lambda)} = \sqrt{\frac{8\pi G}{3}(\rho_\chi + \rho_b + \rho_r)} \approx 40 \text{ km/s/Mpc}$$

**למה זה חשוב**: כדי לדעת כמה $\rho_\Lambda$ חסר. הפער $67.4 - 40 = 27.4$ km/s/Mpc → צריך $\rho_\Lambda / \rho_{total} \approx 0.69$.

**סטטוס**: 🟢 חישוב פשוט — בסיס

---

### נתיב 7️⃣: $H_0$ מ-$\alpha_D$ ישירות (wild card)

**רעיון**: הקשר $\Lambda_d = \sqrt{H_0 M_{Pl}}$ + Friedmann $H_0^2 = 8\pi G \rho_\Lambda / (3\Omega_\Lambda)$ + $\rho_\Lambda = \Lambda_d^4 (1-\cos\theta_i)$ נותן:

$$H_0 = \frac{(1-\cos\theta_i)^{1/2}}{(3\Omega_\Lambda)^{1/2}} \cdot \frac{\Lambda_d^2}{M_{Pl}}$$

אם transmutation (נתיב 2) נותן $\Lambda_d = \Lambda_d(\alpha_d)$, ו-$\alpha_d = \alpha_D$ (unification), אז:

$$\boxed{H_0 = H_0(\alpha_D, \theta_i)}$$

עם $\alpha_D = 5.734 \times 10^{-3}$ ו-$\theta_i = 2$ rad — שני מספרים → $H_0$.

**הבעיה**: $\alpha_D \neq \alpha_d$ (SIDM coupling ≠ dark QCD coupling). אלא אם כן יש unification. 

**סטטוס**: 🔴 ספקולטיבי — צריך מודל unification

---

## טבלת סיכום — סדר עדיפויות

| # | נתיב | $H_0$ prediction? | פרמטרים חופשיים | מעשיות | עדיפות |
|---|---|---|---|---|---|
| **2** | **Transmutation $\Lambda_d(\alpha_d)$** | **כן — אם $\alpha_d = f(\alpha_D)$** | **0-1** | **💻 RG equation** | **⭐⭐⭐ ראשון** |
| **5** | **Trace anomaly** | **כן — same as 2** | **0** | **📐 analytic** | **⭐⭐⭐ ראשון** |
| **3** | **אינטגרל זמן (עומר)** | בהתניית $\rho_\Lambda$ | 1 ($\rho_\Lambda$) | 💻 numeric | ⭐⭐ שני |
| **6** | **Friedmann ללא Λ** | $H_0^{noDE}$ | 0 | 💻 simple | ⭐⭐ שני (baseline) |
| 1 | Dark QCD + misalignment | consistency only | 2 ($\Lambda_d, \theta_i$) | ✅ done (PI-7) | ⭐ כבר נעשה |
| 4 | CW cancellation | ספקולטיבי | ? | 🔴 hard | ⬜ |
| 7 | $\alpha_D$ ישיר | ספקולטיבי | 1 ($\theta_i$) | 🔴 needs unification | ⬜ |

---

## תוכנית עבודה

### שלב א — Baseline (חישוב פשוט)
- [ ] `baseline_friedmann.py`: $H_0$ מ-DM+baryons+radiation בלבד (ללא DE)
- [ ] `baseline_friedmann.py`: $H_0$ עם DE כסריקה על $\rho_\Lambda$ — מה הטווח?

### שלב ב — Transmutation (הכיוון המבטיח ביותר)
- [ ] `rg_transmutation.py`: RG running של $\alpha_d$ מ-$m_\phi$ עד $M_{Pl}$ ב-SU(2)_d
- [ ] בדוק: האם $\alpha_D(m_\phi) = 5.734 \times 10^{-3}$ → $\alpha_d(M_{Pl}) = 1/71$?
- [ ] אם כן: $\Lambda_d$ **נגזר** → $m_\sigma$ **נגזר** → $\rho_\Lambda$ **נגזר** → $H_0$ **נגזר**

### שלב ג — אינטגרל זמן (הרעיון של עומר)
- [ ] `cosmic_time_integral.py`: אינטגרציה של $\rho_{total}(t)$ מ-$t_{BBN}$ עד $t_0$
- [ ] בדוק cancellation דינמי ב-$V_{CW}(T(t))$
- [ ] פלט: $H_0(t)$ כפונקציה של הזמן הקוסמי

### שלב ד — Trace Anomaly
- [ ] `trace_anomaly.py`: $\langle T^\mu_\mu \rangle_{dark}$ מהלגרנזיאן
- [ ] חישוב condensate $\sim \Lambda_d^4$
- [ ] חיבור ל-$\rho_\Lambda$

---

## הקשר בין הנתיבים

```
                    Lagrangian
                   ╱          ╲
          derivative           path integral
         (EOM → VPM)          (Z → propagators)
              ↓                      ↓
         V(r) → σ_T/m         ⟨σv⟩ → Ωχh²
              ↓                      ↓
          SIDM ✅              ρ_DM ✅
                                     ↓
                              ρ_total = ρ_DM + ρ_b + ρ_r + ρ_Λ
                                                            ↑
                                                        ??? ←── THE HUNT
                                                            ↑
              path 2: αd → Λd⁴ ─────────────────────────────┘
              path 5: trace anomaly = Λd⁴ ──────────────────┘
              path 3: ∫dt V_eff(T(t)) ──────────────────────┘
                                     ↓
                                  Friedmann
                                     ↓
                              H₀ = √(8πGρ/3)
```

---

*Date: 28 Mar 2026*
*Location: `dark-energy-T-breaking/hunt_H0/`*
