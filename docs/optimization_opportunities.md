# הזדמנויות ייעול — תהליך חישוב VPM Scan
**תהליך:** `core/v22_raw_scan.py` — Order 1, Raw Parameter Scan
**תאריך:** 25 מרץ 2026
**מצב נוכחי:** ריצה ~7+ שעות (l_max=82, grid 50×70×4×200)

---

## א. לוג קובץ + ניתוח חלקי בזמן ריצה

### בעיה
הסקריפט מדפיס progress רק לחלון נפרד (לא נגיש מה-terminal של VS Code).
אין אפשרות לעקוב אחרי ריצה ברקע בלי לפתוח את החלון ידנית.

### פתרון מוצע
כתוב שורת JSONL לקובץ `data/scan_progress.jsonl` לכל m_chi שמסתיים:
```json
{"ts": "14:32:01", "done": 32, "total": 50, "elapsed_s": 12340, "eta_s": 2300,
 "raw_so_far": 89234, "rep_so_far": 1201,
 "avg_lmax_used": 34.2, "max_lmax_used": 82, "early_exit_pct": 73.1,
 "ms_per_point": 0.22}
```

### שדות חדשים בשורת ה-print הקיימת
```
[32/50] m_chi=45.62 GeV  (12340s, ETA 2300s)  raw: 89,234  rep: 1,201
        l_max: avg=34.2 max=82 | early_exit: 73% | ms/point: 0.22
```

### תועלת
- ניתן לקרוא את `scan_progress.jsonl` בזמן ריצה: `Get-Content data/scan_progress.jsonl -Tail 1`
- ניתוח חלקי אפשרי לפני סיום: כמה נקודות ברות-קיימא יש כבר, l_max ממוצע בפועל
- ETA מדויק יותר בריצות עתידיות על בסיס הנתונים מהריצה הנוכחית

---

## ב. הערכת זמן ריצה לפני תחילה (Dry-run mode)

### בעיה
לא ניתן לדעת כמה תיקח ריצה לפני שמתחילים — עד עכשיו למדנו בדיעבד.

### פתרון מוצע — מצב `--estimate`
רץ **2% מה-grid** (5 m_chi × 7 m_phi × sample α) ומחשב extrapolation:

```
$ python core/v22_raw_scan.py --estimate

[Dry-run: 5/50 m_chi × 7/70 m_phi × 10/200 alpha]
Sampling 350 points (0.1% of grid)...

Results:
  avg l_max used:      34.2  (max observed: 82)
  early_exit_rate:     71%
  ms per point (fast): 0.19ms
  ms per point (slow): 2.1ms

Estimated runtime:
  Optimistic (all fast): ~1.8h
  Realistic  (71% fast): ~4.2h   ← best estimate
  Pessimistic (all slow): ~19.5h

Heavy workers expected: ~6/50 m_chi (12%)
Peak l_max expected:    ~82 at m_chi ≈ 150-300 GeV
```

---

## ג. מצבי דיוק (Accuracy Presets)

### בעיה
כרגע יש רק הגדרה אחת. לא ניתן לבחור trade-off בין מהירות לדיוק.

### פתרון מוצע — הוסף ל-`global_config.json`
```json
"accuracy_preset": "fast"
```

| Preset | l_max buffer | N_steps (max) | early_exit threshold | זמן צפוי | שגיאה ברזוננסים |
|---|---|---|---|---|---|
| `fast`        | +5  (במקום +20) | 4,000  | 1e-3 | ~45 דק'  | ~3-5%   |
| `standard`    | +20             | 8,000  | 1e-4 | ~3h      | ~0.5-1% |
| `publication` | +30             | 12,000 | 1e-5 | ~8h      | <0.1%   |

**הדפסה בתחילת ריצה:**
```
Accuracy preset: fast
  l_max buffer: +5    →  truncation error ~3-5% near resonances
  N_steps: 4000       →  RK4 discretization error <0.1%
  early_exit: 1e-3    →  partial-wave truncation <0.3%
  TOTAL expected error: ~3-5% (acceptable for parameter scan)
```

---

## ד. שיפורי דיוק ומהירות אלגוריתמיים

### ד.1 Early-exit אדפטיבי לפי proximity לרזוננס
**רעיון:** סף ה-early-exit תלוי בעוצמת האינטראקציה:
```python
# רחוק מרזוננס (sigma קטן)
if peak_contrib < 10:   threshold = 1e-3  # 2x מהיר מ-standard
# קרוב לרזוננס (sigma גדול)
elif peak_contrib > 100: threshold = 1e-5  # יותר מדויק
else:                    threshold = 1e-4  # standard
```
**תועלת:** שיפור ~30% בזמן ריצה ללא אובדן דיוק.

### ד.2 שמירת epsilon לכל נקודה בפלט CSV
**רעיון:** הוסף עמודה `lmax_epsilon` לכל שורה בפלט:
```python
epsilon = contrib_last / sigma_sum  # שגיאה יחסית אמיתית
```
**תועלת:** לכל נקודה בפלט יודעים בדיוק מה רמת הדיוק שלה. ניתן לסנן נקודות עם epsilon גבוה.

### ד.3 Caching לנקודות סימטריות
**רעיון:** אם שתי נקודות שונות (m_chi, m_phi, α) נותנות אותם (κ, λ) — חשב פעם אחת.
**מגבלה:** דורש lookup table — מורכב יחסית ליישום.

### ד.4 Early termination של m_chi שלם
**רעיון:** אם לאחר 10% מנקודות ה-m_phi שום נקודה לא עוברת את ה-SIDM cuts — דלג על שאר.
**תועלת:** ייתכן שחלק מה-50 m_chi values לא אמורים לתת תוצאות כלל.

---

## ה. מדדים להוסיף ל-`runs_log.csv` בסיום ריצה

| עמודה | תוכן |
|---|---|
| `avg_lmax` | ממוצע l_max בפועל לכל הנקודות |
| `max_lmax` | מקסימום l_max שנצפה |
| `early_exit_pct` | אחוז נקודות שסיימו early-exit |
| `heavy_workers` | כמה m_chi workers היו "כבדים" (>2x זמן ממוצע) |
| `ms_per_point_avg` | מילישניות ממוצע לנקודה |

**תועלת:** נוסחת ETA מדויקת לריצות עתידיות:
$$T_{\text{עתידי}} = N_{\text{grid}}^{\text{חדש}} \times \overline{t_{\text{point}}} \times \frac{\overline{l_{max}}^{\text{חדש}}}{\overline{l_{max}}^{\text{נוכחי}}} \div N_{\text{workers}}$$

---

## ז. סריקת רגישות ל-SIDM Cuts (חדש)

### רקע
ב-commit `0f73dec` שונה `sigma_m_1000_hi` מ-0.1 ל-0.47 cm²/g (Harvey+2015, 72 cluster mergers).
שינוי זה הכפיל/שילש את מספר הנקודות הכשירות — מה שגם חשף את bug ה-l_max.
הערך 0.47 הוא **95% CL upper limit** — לא ערך מוחלט. קיים אי-ודאות אובזרבציונלי אמיתי.

### ערך מחקרי
- **ניתוח רגישות** (sensitivity analysis): כיצד ה-viable parameter space משתנה ל-σ/m(cluster) שונה?
- הוכחה שהתוצאות (MAP, relic abundance, phenomenology) חזקות ביחס לאי-ודאות ב-cluster constraint.
- **ציטוט**: Robertson+2021, Sagunski+2021, Andrade+2022 נותנים ערכים 0.1–1.5 cm²/g (model-dependent, relaxation vs. scattering).

### ערכים לסריקה
```python
sigma_m_1000_hi_values = [0.1, 0.2, 0.3, 0.47, 0.7, 1.0]  # cm^2/g
```
כל ריצה = pipeline מלא (raw_scan → relic → MCMC) עם ערך שונה של הפרמטר.

### השפעה על זמן ריצה
- ריצה אחת: ~9 שעות (עם l_max נכון)
- 6 ערכים × 9 שעות = ~54 שעות (אם sequential). אפשר להריץ 2-3 במקביל.
- ריצה עם ערך נמוך (0.1): מהירה יותר — פחות נקודות עוברות filter

### מה לתעד לכל ריצה
הוסף ל-`runs_log.csv`:
- `sigma_m_1000_hi` כפרמטר מפתח
- `n_viable` (כמה נקודות עברו את ה-filter)
- MAP point per run (עשוי להשתנות!)

### עדיפות
**גבוהה** — זה הנתיב הישיר לסקציה "SIDM Constraint Robustness" בפרסום.

---

## ח. ריצות מקביל על שירות ענן

### בעיה
סריקת 6 ערכי `sigma_m_1000_hi` sequentially = ~54 שעות wall time.
המחשב תפוס ואי אפשר להריץ שום דבר אחר.

### פתרון: Cloud Spot Instances

#### Modal.com — הכי פשוט לPython
- API Python-native — מעלים קוד ישירות, ללא setup server
- משיק 6 jobs במקביל בפקודה אחת
- Serverless: משלמים רק זמן ריצה אמיתי
- Numba+Linux מהיר ~15% יותר מ-Windows
- **עלות אומדנת**: ~$0.15/h × 9h × 6 runs = **~$8–12 סה"כ**
- **Wall time**: ~5–6 שעות (במקום 54h)

#### AWS EC2 Spot Instances — הכי זול
| Instance | Cores | מחיר spot | עלות ×6 ריצות |
|---|---|---|---|
| c5.4xlarge | 16 | ~$0.07/h | ~$4 |
| c6i.8xlarge | 32 | ~$0.12/h | ~$6 |

- עם 32 cores: כל ריצה מתקצרת ל-~5h (יותר workers)
- צריך: git clone + pip install, S3 לפלטים

#### Vast.ai / RunPod — זול + ממשק פשוט
- CPU-only: $0.03–0.05/h לcore
- ממשק GUI, לא צריך AWS account
- **עלות**: ~$3–5 סה"כ ל-6 ריצות

### תנאי מוקדם לענן
הסקריפטים צריכים לקבל `sigma_m_1000_hi` כ-CLI argument (`--sigma-cluster 0.47`)
במקום לקרוא ישירות מ-`global_config.json` — כדי שכל instance ירוץ עם ערך שונה.
זה שינוי קטן (~5 שורות) שגם מועיל לריצות מקומיות.

### עדיפות
**גבוהה** — מאפשר להשלים את כל ניתוח הרגישות ב-5 שעות במקום שבוע.
עלות: ~$5–12 חד-פעמי.

---

## ו. סיכום עדיפויות

| שיפור | השפעה על מהירות | השפעה על דיוק | קושי יישום |
|---|---|---|---|
| לוג JSONL לקובץ | — | — | נמוך ⭐ |
| הדפסת avg_lmax + early_exit% | — | — | נמוך ⭐ |
| Accuracy presets | ×3-5 מהיר יותר | מתכוונן | נמוך ⭐⭐ |
| --estimate dry-run | — | — | בינוני ⭐⭐ |
| Early-exit אדפטיבי | ×1.3 מהיר יותר | +שיפור | בינוני ⭐⭐ |
| epsilon בפלט CSV | — | ידוע מדויק | נמוך ⭐ |
| Early termination m_chi | ×? | — | בינוני ⭐⭐ |
| Richardson extrapolation | ×0.5 איטי יותר | ×10 יותר מדויק | גבוה ⭐⭐⭐ |
| CLI arg לפרמטרי cuts | תנאי מוקדם לענן | — | נמוך ⭐ |
| ריצות ענן מקביל (Modal/AWS) | ×10 wall time | — | בינוני ⭐⭐ |

**המלצה לריצה הבאה:** יישם קודם את הלוג + epsilon + presets — נמוך בקושי, גבוה בתועלת.

---

## ט. TODO — GPU FP32 (אחרי אימות fast)

### מטרה
להאיץ את סקן הרגישות (6 ערכי sigma_m_1000_hi) מ-~3.5h ל-~30 דקות.

### תנאי מוקדם — אימות דיוק FP32 vs FP64
לפני שמשתמשים ב-FP32 לתוצאות, חובה לוודא שהשגיאה קטנה מסף:
```
python core/validate_fp32_accuracy.py
```
בודק ~50 נקודות representative (כולל ליד רזוננסים בהן השגיאה הכי גדולה):
- `sigma_T_vpm_fp32(mc, mp, alpha, v)` vs `sigma_T_vpm(mc, mp, alpha, v)` (FP64)
- קריטריון: שגיאה יחסית < 1% בכל הנקודות
- אם עוברים → GPU FP32 בטוח לסקן רגישות
- אם נכשלים → FP32 לא בטוח ליד רזוננסים, להשתמש רק ל-coarse scan

### סיכון עיקרי
ליד רזוננס, $\delta_l \approx \pi/2$ ו-$\sin^2(\delta_l)$ רגיש מאוד לשגיאת עיגול.
FP32 מצבר ~$10^{-7}$ per RK4 step × 4000-12000 steps → עד $10^{-3}$ בפאזה.
זה עלול ל-shift את מיקום הרזוננס ולפגוע בספירת viable points.

### קובץ ליצור: `core/v22_raw_scan_gpu_fp32.py`
- `@cuda.jit` kernel — RK4 + partial wave sum, FP32
- כל thread = נקודה אחת (mc, mp, alpha, v)
- 2.8M threads ל-grid מלא
- thread divergence בגלל `l_max` שונה per thread — לנהל עם `__syncwarp`

### עדיפות
**בינונית** — רלוונטי רק לסקן רגישות (6 ריצות × ~35 min = 3.5h כבר סביר).
לא לפני שהרצת fast הראשונה עוברת ומאומתת.

**לאחר השלמת Order 1:** הוסף CLI args לסקריפטים → הפעל pipeline מלא על `sigma_m_1000_hi = [0.1, 0.2, 0.3, 0.47, 0.7, 1.0]` דרך Modal או AWS (~$10, ~5h wall time).
