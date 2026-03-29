# יומן מחקר — מה זה אומר
**תיקייה:** `what_does_it_mean/`  
**מטרה:** הבנה פיזיקלית מעמיקה של המודל — מה הלגרנזיאן *אומר* בפועל, שלב אחר שלב.

---

## 2026-03-27 — מסקנות מאינטגרל הדרך של הלגרנזיאן

### נקודת המוצא

$$\mathcal{L} = \frac{1}{2}(\partial_\mu\phi)^2 - \frac{1}{2}m_\phi^2\phi^2 + \bar{\chi}(i\!\not\!\partial - m_\chi)\chi - g_s\,\phi\bar{\chi}\chi$$

$$\alpha = \frac{g_s^2}{4\pi}, \qquad \lambda = \frac{\alpha m_\chi}{m_\phi}$$

---

### שלב 1 — חוקי פיינמן מאינטגרל הדרך

פונקציית היצירה:
$$Z[J,\eta,\bar\eta] = \int \mathcal{D}\phi\,\mathcal{D}\chi\,\mathcal{D}\bar\chi\;\exp\!\left(i\int d^4x\,\mathcal{L} + \text{sources}\right)$$

האינטגרציה הגאוסיאנית על האיבר הריבועי (שדות חופשיים) מניבה:

| אובייקט | ביטוי |
|---------|-------|
| פרופגטור סקלרי | $\Delta_F(q^2) = \dfrac{i}{q^2 - m_\phi^2 + i\varepsilon}$ |
| פרופגטור פרמיוני | $S_F(p) = \dfrac{i(\not\!p + m_\chi)}{p^2 - m_\chi^2 + i\varepsilon}$ |
| ורטקס אינטראקציה | $-ig_s$ (סקלר × ביליניאר מאיורנה) |

---

### שלב 2 — פוטנציאל יוקאווה

שני פרמיונים מחליפים סקלאר אחד (ערוץ $t$ + ערוץ $u$ עבור מאיורנה).  
בגבול NR: $t = -|\mathbf{q}|^2$, $\hat{V}(q) = -4\pi\alpha/(|\mathbf{q}|^2 + m_\phi^2)$.  
טרנספורמציית פורייה:

$$\boxed{V(r) = -\frac{\alpha}{r}\,e^{-m_\phi r}}$$

**מה זה אומר:** הורטקס $-ig_s$ יוצר **כוח אטרקטיבי** עם:
- טווח $r_0 = \hbar c / m_\phi$
- BP1: $r_0 \approx 15$ fm (~13× רדיוס פרוטון)
- MAP: $r_0 \approx 20$ fm (~17× רדיוס פרוטון)

---

### שלב 3 — קירוב בורן (תקף רק כאשר λ ≪ 1)

מקביעה ל-Born rule ופיתוח:

$$\sigma_T^{\rm Born,\,Maj} = \frac{16\pi\alpha^2}{m_\chi^2 v^4}\,f(\beta), \qquad \beta \equiv \frac{m_\chi v}{m_\phi}$$
$$f(\beta) = \ln(1+\beta^2) - \frac{\beta^2}{1+\beta^2}$$

**בדיקת תוקף** ב-$v = 30$ km/s:

| מודל | $\lambda$ | $\sigma_T^{\rm Born}$ [cm²/g] | $\sigma_T^{\rm VPM}$ [cm²/g] | Born/VPM |
|------|-----------|-------------------------------|------------------------------|----------|
| BP1 | 11.1 | 5.90e+01 | 6.70e-01 | **88×** — לא תקף |
| MAP | 33.3 | 2.44e+02 | 1.81e+00 | **135×** — לא תקף |

> **מסקנה:** Born מגזים פי 88–135 עבור הנקודות שלנו. חייבים VPM.

---

### שלב 4 — VPM (Variable Phase Method)

מאחר ש-$\lambda = 11$–$33 \gg 1$, חייבים לסכם גלים חלקיים מלאים:

$$\sigma_T = \frac{2\pi}{k^2}\sum_\ell w_\ell(2\ell+1)\sin^2\delta_\ell$$

- $w_\ell = 1$ עבור $\ell$ זוגי, $w_\ell = 3$ עבור $\ell$ אי-זוגי (סימטריית זהות מאיורנה)
- ה-prefactor $2\pi$ (ולא $4\pi$) נובע מהגורם $\frac{1}{2}$ של חלקיקים זהים
- הפאזה $\delta_\ell$ נגזרת מ-ODE של VPM עבור פוטנציאל יוקאווה

**תוצאה:** $\sigma_T(v)$ גדול ב-dwarf ($v \sim 30$ km/s), קטן ב-cluster ($v \sim 1000$–$5000$ km/s) → בדיוק מה שצפינו

---

### שלב 5 — אניהלציה χχ → φφ

ערוצי $t + u$ (חילוף פרמיון) עבור מאיורנה:

$$\langle\sigma v\rangle_{s\text{-wave}} = \frac{\pi\alpha^2}{4m_\chi^2} \qquad (m_\phi/m_\chi \to 0)$$

| מודל | $\langle\sigma v\rangle$ [cm³/s] | / Planck relic |
|------|----------------------------------|----------------|
| BP1 | 2.155e-26 | **0.72×** |
| MAP | 1.019e-26 | **0.34×** |

> BP1 קרוב לצפיפות הרליק של Planck ($3\times10^{-26}$ cm³/s).  
> MAP צריך עיבוי נוסף (Sommerfeld enhancement בזמן הקפאה).

---

### סיכום — מה הלגרנזיאן אומר פיזיקלית

| איבר בלגרנזיאן | משמעות פיזיקלית |
|----------------|-----------------|
| $\frac{1}{2}(\partial\phi)^2 - \frac{1}{2}m_\phi^2\phi^2$ | מדיאטור סקלרי עם מסה $m_\phi$ — קובע את טווח הכוח |
| $\bar\chi(i\!\not\!\partial - m_\chi)\chi$ | DM פרמיון מאיורנה — חלקיק = אנטי-חלקיק |
| $-g_s\phi\bar\chi\chi$ | הנקודה שממנה הכל נגזר: כוח יוקאווה + אניהלציה |

- **ורטקס אחד** → שני תהליכים: פיזור (SIDM) + אניהלציה (relic density)
- **הכוח אטרקטיבי** ← סימן מינוס ב-$V(r)$ ← שני ערוצי $t+u$ של מאיורנה מחזקים
- **המשטר רזוננטי** ($\lambda \gg 1$) הוא ה*יתרון* של המודל: מאפשר $\sigma/m$ גדול ב-dwarfs בלי לחרוג ב-clusters
- **Born לא מספיק** — ה-VPM הוא הכרחי, לא בחירה

---

**קבצים רלוונטיים לתיקייה זו:**
- `lagrangian_to_observables.py` — קוד שרץ את כל השרשרת
- `output/lagrangian_to_observables.png` — גרף שלושת השלבים

---

---

## 2026-03-27 — מה עוד ניתן לחקור דרך האינטגרל?

### הנקודה של עומר

כשגוזרים את $\mathcal{L}_{SIDM}$, הפעולה על קבועים היא אפס — הקבוע הקוסמולוגי **נאבד בנגזרת**.  
האינטגרל $\int\mathcal{D}\mathcal{L}$ **לא** מאבד קבועים — הוא עושה בדיוק ההפך: מצטבר עליהם.  
לכן ביצוע האינטגרל (path integral) על $\mathcal{L}_{SIDM}$ הוא הדרך הנכונה לשאול מה $\rho_\Lambda$ שייך לסקטור האפל.

### מפת השוואה — נגזרת vs. אינטגרל

| מה עשינו עם הנגזרת | מה מקביל לו דרך האינטגרל |
|---------------------|--------------------------|
| EOM → פוטנציאל יוקאווה → VPM → σ_T/m | אינטגרל fluctuation → תיקונים לσ_T מעבר ל-VPM |
| V_CW(θ) — אנרגיית ריק בנקודה קבועה | $Z(T)$ — אנרגיית ריק **כפונקציית טמפרטורה** |
| כוח חמישי β מ-$\partial\mathcal{L}/\partial\sigma$ | V_eff(σ, n_χ) מאינטגרל ב-medium |
| בדיקת לכידה ב-θ_relic (ODE) | P(θ_i) — הסתברות לקבלת θ_i מהינפלציה |
| יציבות ראדיאטיבית של σ (loop corrections) | RG זרימה ← אינטגרל Wilson |
| θ כפרמטר טופולוגי/A₄ | אינסטנטון ← tunneling בין θ=0 ל-θ=π/2 |

---

### Test PI-1: פוטנציאל תרמי $V_{eff}(\sigma, T)$ — האינטגרל כפונקציית טמפרטורה

**שאלה**: בטמפרטורת ה-freeze-out ($T_{fo} \sim m_\chi/20$), מה הפוטנציאל האפקטיבי שσ חש?  
**מה Test 4 עשה**: הוסיף ידנית "thermal backreaction" כ-$\sim n_\chi \cdot dM_{eff}/d\theta$ — אנליטי וגס.  
**מה האינטגרל נותן**: פוטנציאל תרמי מלא ב-Matsubara:

$$Z(T) = \int_{\text{periodic}} \mathcal{D}\phi\,\mathcal{D}\chi\; e^{-S_E[\phi,\chi]}$$

$$V_{eff}(\sigma, T) = V_{CW}(\sigma) + \frac{T^4}{2\pi^2}\int_0^\infty dp\; p^2\left[\ln\left(1 - e^{-E_\phi/T}\right) - \ln\left(1 + e^{-E_\chi/T}\right)\right]$$

**שאלת המחקר**:  
האם ב-$T = T_{fo}$ יש מינימום ב-$V_{eff}(\sigma, T)$ ב-$\theta_{relic}$?  
Test 4 מצא ש-CW מגלגל ל-$\pi/2$ — פוטנציאל תרמי יכול לשנות זאת בצורה דרמטית בסמוך ל-$T \sim m_\chi$.

**✗ נבדק — `test_PI1_thermal_veff.py` (commit `9833519`)**

| גודל | ערך ב-$T_{fo}$ | משמעות |
|------|---------------|--------|
| $\theta_{min}(T_{fo})$ | $1.57$ rad = $\pi/2$ | CW שולט, תרמי זניח |
| $\theta_{min}(T=0)$ | $1.57$ rad = $\pi/2$ | CW לא משתנה |
| $\theta_{min}(T=10m_\chi)$ | $0.01$ rad | תרמי דוחף ל-$\theta=0$ |
| $V_{QCD}/V_{CW}$ ב-$T_{fo}$ | $6.4\times10^{-116}$ | QCD מדוכא לחלוטין |
| $V_T/V_{CW}$ ב-$T_{fo}$ | $1.2\times10^{-4}$ | תרמי זניח (Boltzmann suppressed) |

**סיבה פיזיקלית**: $V_{CW}$ מינימום ב-$\theta=\pi/2$ כי $m_\chi(\pi/2) = m_\chi/3$ → אנרגיית לולאה נמוכה ביותר. הפוטנציאל התרמי $V_T$ מנסה לדחוף ל-$\theta=0$ בטמפ' גבוהות, אבל **Boltzmann-suppressed כבר ב-$T=T_{fo}$** כי $m_\chi/T_{fo} = x_{fo} \approx 20$.

**מסקנה**: **אין attractor תרמי ב-$\theta \sim 2$ rad.** $V_{eff}(\theta, T)$ תמיד ממזערת ב-$\theta=\pi/2$, לא ב-$\theta_i$.

**✗ PI-1 NEGATIVE RESULT — θ_i CONFIRMED AS FREE PARAMETER**

---

### Test PI-2: אינסטנטון — האם θ_A₄ יציב מול tunneling?

**שאלה**: גם אם A₄ קובע $\theta = \arcsin(1/3)$ כ"מינימום" UV, האם הוואקום הקוונטי יכול לנהור ממנו?  
**כיצד**: דרך ה-Euclidean path integral מחפשים פתרון אינסטנטון:

$$\frac{d^2\sigma}{d\tau^2} = \frac{dV_{CW}}{d\sigma}$$

עם תנאי גבול $\sigma(\tau \to \pm\infty) = \theta_{relic}$, $\sigma(\tau=0) = \theta_{\pi/2}$ (נקודת ה"עבר").

**שיעור הנהירה**:
$$\Gamma_{tunnel} \sim \exp\left(-S_E^{instanton}\right), \quad S_E \sim \frac{\Delta V \cdot r_0^4}{\hbar}$$

כאשר $\Delta V = V_{CW}(\pi/2) - V_{CW}(\theta_{relic})$ ו-$r_0 = 1/m_\sigma$.

**למה זה חשוב**: אם $\Gamma_{tunnel} \ll H_0^4$ — θ יציב לגיל היקום ✅.  
אם $\Gamma_{tunnel} \gg H_0^4$ — A₄ לא מחזיק את הזווית על פני 14 מיליארד שנה ❌.

**תוצאות (committed `8873be1`)**:

| שיטה | $\log_{10}(S_E)$ | יציב? |
|------|----------------|-------|
| Dimensional estimate | 120.9 | ✅ כן |
| Thin-wall (analytic) | 121.5 | ✅ כן |
| Thin-wall (numeric) | 121.5 | ✅ כן |

- $S_E \sim (f/m_\sigma)^2 \sim (2\times10^{60})^2 \sim 10^{121}$
- $S_{E,crit} \approx 0$ (סף יציבות = כמעט כלום)
- $\log_{10}(N_{bubbles}) \approx -1.4\times10^{121}$ — **אפס מוחלט** של אירועי tunneling

**מסקנה פיזיקלית**: דופן הבועה חייבת לפרוש את $r_0 \sim 1/m_\sigma \sim H_0^{-1}$ — **רדיוס הוביל**. אינסטנטון שכזה אסור קוסמולוגית, לא מדוכא אנרגטית. $\theta_{A_4} = \arcsin(1/3)$ קפוא לכל גיל היקום.

**✅ PI-2 PASS — ABSOLUTELY STABLE**

---

### Test PI-3: Stochastic Path Integral — P(θ_i) מהינפלציה

**שאלה**: מה הסבירות שσ/f ינחת קרוב ל-2 rad (הדרוש ל-H₀ = 67.4 ב-PI-7)?  
**המנגנון**: בינפלציה, כל שדה סקלרי עובר **diffusion קוונטי** ב-de Sitter:

$$\langle(\delta\sigma)^2\rangle = \frac{H_{inf}^2}{4\pi^2} \quad \text{per e-fold}$$

לאחר $N$ e-folds: $P(\theta_i) = \mathcal{N}\exp\left(-\frac{\theta_i^2}{2\sigma_{rms}^2}\right)$ כאשר $\sigma_{rms} = \frac{H_{inf}}{2\pi f}\sqrt{N}$.

**✗ נבדק — `test_PI3_stochastic_inflation.py` (commit `8fb75b9`)**

| גודל | ערך | משמעות |
|------|-----|--------|
| $\sigma_{rms}$ (H_inf=10¹³ GeV) | $4.2\times10^{-6}$ rad | כמעט אפס |
| $\theta_i / \sigma_{rms}$ | $4.75\times10^5$ — רחוק $5\times10^5$ סטיות תקן | כיוון עצום |
| $P(\theta_i=2)$ | $\sim 10^{-4.9\times10^{10}}$ | fine-tuning קיצוני |
| $H_{inf}$ לטבעיות 1σ | $4.75\times10^{18}$ GeV | 5 סדרי גודל **מעל** גבול Planck |
| Drift/Diffusion ratio | $2.3\times10^{-98}$ | פוטנציאל $V(\theta)$ זניח לחלוטין |

**מסקנה**: $\theta_i = 2$ rad **לא** נוצר טבעית על-ידי de Sitter diffusion בסקלות אינפלציה CMB-מותרות.

$\theta_i$ הוא **פרמטר חופשי** של המודל (זווית misalignment קלאסית).

**השלכה על PI-7**: H₀ = 67.4 הוא consistency check נתון $\theta_i$, לא תחזית ללא פרמטר חופשי — **אלא אם PI-1 יראה attractor תרמי שמקבע $\theta_i \sim 2$ rad**.

**שלוש אפשרויות**:
1. (A) $\theta_i$ כתנאי התחלה קלאסי — H₀ = consistency check
2. (B) attractor תרמי ב-$V_{eff}(\theta, T)$ — PI-1 יבדוק זאת
3. (C) anthropic/landscape: patches עם $\theta_i \sim 2$ נותנות $\rho_\Lambda$ הנכונה

**✗ PI-3 NEGATIVE RESULT — θ_i IS A FREE PARAMETER**

---

### Test PI-4: Wilsonian RG — מה קורה לα בזרימה ל-UV?

**שאלה**: האם $\alpha_{dark}$ ב-UV מגיע מ-UV fixed point, Landau pole, או נגמר בחלון ה-SIDM?  
**שיטה**: Wetterich equation (functional RG):

$$\partial_t \Gamma_k = \frac{1}{2}\text{Tr}\left[\left(\Gamma_k^{(2)} + R_k\right)^{-1} \partial_t R_k\right]$$

כאשר $t = \ln(\Lambda/k)$, $\Gamma_k$ הוא ה-1PI effective action עם regulator $R_k$ שחוסם מומנטומים נמוכים.

**מה זה אומר בפועל**: טרייל אינטגרציה על שכבות של מומנטום מ-$m_\phi$ עד $M_{Pl}$. בכל שכבה: איך $\alpha$ משתנה?

**שאלת המחקר**:  
- האם $\alpha_{dark}$ זורם לנקודה קבועה ב-UV? (כמו $\alpha_{QCD}$ — קטן ב-UV, גדל ב-IR)  
- או שיש Landau pole לפני $M_{Pl}$? (אז המודל לא UV-complete)  
- האם ה-A₄ symmetry מגן על ה-ratio $\alpha_p/\alpha_s$ תחת RG flow?

---

### Test PI-5: In-Medium Path Integral — V_eff(σ, n_χ) בגלקסיה

**שאלה**: בתוך גלקסייה ננסית ($n_\chi \gg n_{\chi,cosmological}$), האנרגיה האפקטיבית של σ משתנה?  
**שיטה**: אינטגרל עם צפיפות ב-medium:

$$V_{eff}(\sigma, n_\chi) = V_{CW}(\sigma) + n_\chi \cdot M_{eff}(\sigma)$$

כאשר $M_{eff}(\sigma) = \sqrt{m_\chi^2 + y^2\phi^2/4}$ × פונקציית θ.

**מה זה אומר**: מסת σ אפקטיבית תלויה בצפיפות המקומית:

$$m_{\sigma,eff}^2(n_\chi) = \left.\frac{d^2 V_{eff}}{d\sigma^2}\right|_{\sigma_0} = m_{\sigma,CW}^2 + n_\chi \frac{d^2 M_{eff}}{d\sigma^2}$$

**שאלת המחקר** (הנגזרת ממה שמצאנו ב-Test 3):  
Test 3 טען שהמסה CW $\sim 10^{-12}$ eV מגינה מהכוח החמישי — אבל זה היה הערכה. האינטגרל ב-medium נותן את הטיפול **המדויק** של ה-chameleon screening:

- בחלל ריק: $m_{\sigma,eff} = m_{\sigma,CW} \sim 10^{-12}$ eV → טווח 200 ק"מ
- בגלקסיה ננסית: $m_{\sigma,eff} \gg m_{\sigma,CW}$ → טווח קצר עוד יותר → σ מכווץ (screened)
- תיקוף Test 3 מהאינטגרל, ולא מהנגזרת בלבד

---

### Test PI-6: Fluctuation Determinant — תיקוני σ_T מעבר ל-VPM

**שאלה**: האם יש תיקונים ל-$\sigma_T$ שה-VPM מחמיץ?  
**שיטה**: אינטגרל הדרך סביב הפתרון הקלאסי (פוטנציאל יוקאווה):

$$Z_{scatter} = e^{iS_{cl}[\phi_{cl}]} \cdot \underbrace{\left[\det\left(-\partial^2 - m_\phi^2 - g_s^2 n_\chi(r)\right)\right]^{-1/2}}_{\text{fluctuation determinant}}$$

ה-determinant שווה לסכום על כל eigenvalues של אופרטור פיזור — שהוא **בדיוק** ה-VPM phase shift sum!

$$\ln\det = \text{Tr}\ln = \sum_\ell (2\ell+1) \ln e^{2i\delta_\ell} = 2i\sum_\ell (2\ell+1)\delta_\ell$$

**מסקנה אנליטית**: ה-VPM הוא **אינטגרל ה-Gaussian הקלאסי** של path integral הפיזור — לא אפרוקסימציה, אלא **חישוב מדויק** של ה-one-loop correction לאמפליטודה. לא צריך ללכת מעבר אלא אם כן $\lambda \gg \lambda_{crit}$ עם non-perturbative bound states.

---

### סיכום מפת הדרכים

| Test | שאלה | מה הנגזרת נתנה | מה האינטגרל מוסיף |
|------|------|----------------|-------------------|
| PI-1 | האם θ נלכד בַּ-freeze-out? | ✗ **FAIL (negative)** | $\theta_{min}(T_{fo})=\pi/2$; $V_T/V_{CW}\sim10^{-4}$; θ_i **פרמטר חופשי** |
| PI-2 | האם θ_A₄ יציב קוונטית? | ✅ נבדק — PASS | $S_E\sim10^{121}$, bubble wall = Hubble radius → **ABSOLUTELY STABLE** |
| PI-3 | מאיפה θ_i ~ 2 rad? | ✗ **FAIL (negative)** | $\sigma_{rms}=4\times10^{-6}$ rad; gap 5 סדרי גודל; θ_i **פרמטר חופשי** |
| PI-4 | האם α ב-UV תקין? | לא נבדק | RG flow עד $M_{Pl}$ |
| PI-5 | האם screening ב-גלקסיה מדויק? | Test 3: הערכה | V_eff(n_χ): exact chameleon |
| PI-6 | האם VPM exact? | Test 15: σ_p זניח | ✅ VPM = fluctuation determinant |

**סדר עדיפות מוצע**:
1. ~~**PI-2 (אינסטנטון)**~~ — ✅ **DONE** — $S_E\sim10^{121}$, θ_A₄ יציב לחלוטין. `8873be1`
2. ~~**PI-3 (stochastic inflation)**~~ — ✗ **DONE (negative)** — θ_i פרמטר חופשי; gap 5 סדרי גודל מ-CMB. `8fb75b9`
3. ~~**PI-1 (V_eff תרמי)**~~ — ✗ **DONE (negative)** — $\theta_{min}(T_{fo})=\pi/2$; אין attractor; θ_i פרמטר חופשי. `9833519`
4. **Test 21 (FIMP/relic density)** — ✗ **DONE** — `<σv>_bare = πα_D²/4m_χ² = 3.4×10⁻²⁰ cm³/s >> 3×10⁻²⁶` (Planck target) פקטור 10⁶; κ=5.3×10⁻⁴ → dark sector in thermal equilibrium, NOT FIMP; Sommerfeld-enhanced Boltzmann required for Ω h²=0.12. `efe9b8a`
5. **PI-5 (in-medium)** — מאמת Test 3 כמותית
6. **Test 22 (ΔN_eff + φ→2σ)** — 🔴 **NEXT** — BBN/CMB constraint on dark sector dof
