"""One-shot VPM check: MAP parameters, no cooking, exact tool from Test 6."""
import sys
sys.path.insert(0, '.')
from consistency_check_sidm import sigma_T_vpm

# MAP parameters — GeV (VPM expects GeV)
m_chi  = 94.07e-3   # GeV  (= 94.07 MeV)
m_phi  = 11.10e-3   # GeV  (= 11.10 MeV)
alpha  = 5.734e-3   # dimensionless

print("VPM solver — MAP (94.07 MeV, 11.10 MeV, alpha=5.734e-3)")
print("Tool: sigma_T_vpm() from consistency_check_sidm.py (Test 6, verbatim)")
print("Parameters: UNCHANGED. No cooking.")
print()
print(f"  {'v':>10}  {'sigma/m':>16}  target                   result")
rows = [
    (30,    0.5,  50.0,  "0.5-50   dwarfs"),
    (200,   0.0,  35.0,  "<35      MW-sat"),
    (1000,  0.1,   1.0,  "0.1-1.0  clusters"),
    (2000,  0.0,   3.0,  "<3.0     MACS"),
    (3000,  0.0,   1.25, "<1.25    Bullet Cluster"),
]
passed = 0
for v, lo, hi, tgt in rows:
    sm = sigma_T_vpm(m_chi, m_phi, alpha, float(v))
    ok = "PASS" if lo <= sm <= hi else "FAIL"
    if ok == "PASS":
        passed += 1
    print(f"  {v:>10}  {sm:>16.4f}  {tgt:25s}  [{ok}]")

print()
print(f"  Result: {passed}/{len(rows)} constraints passed")
print()
if passed == len(rows):
    print("  CONCLUSION: MAP is SIDM-consistent. No MCMC rerun needed.")
elif passed >= 3:
    print("  CONCLUSION: MAP partially consistent. Check failures carefully.")
else:
    print("  CONCLUSION: MAP is NOT SIDM-consistent. Bullet Cluster constraint needed in MCMC.")
