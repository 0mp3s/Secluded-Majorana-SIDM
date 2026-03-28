import csv, sys, os

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
f_path = os.path.join(_ROOT, 'The_derivative_of_Lagernizan_SIDM', 'data', 'archive',
                      'T10_dirac_mimicry_2026_03_27_r033.csv')

with open(f_path) as f:
    rows = list(csv.DictReader(f))

bps = {}
for r in rows:
    bp = r['bp']
    if bp not in bps:
        bps[bp] = {'chi2': float(r['chi2_per_dof']), 'devs': []}
    bps[bp]['devs'].append((float(r['v_km_s']), abs(float(r['deviation_pct']))))

print("BP     chi2/dof   max_dev%   at_v(km/s)   status")
print("-"*55)
for bp in bps:
    mv, md = max(bps[bp]['devs'], key=lambda x: x[1])
    chi2 = bps[bp]['chi2']
    tag = 'DISTINGUISHABLE' if chi2 > 0.01 else 'DEGENERATE'
    print(f"{bp:6s} {chi2:.5f}    {md:6.2f}%     {mv:6.0f}       {tag}")

print()
print("Deviation at each v for each BP:")
print(f"{'v(km/s)':>12}", end="")
for bp in bps:
    print(f"  {bp:>8}", end="")
print()
all_v = sorted(set(float(r['v_km_s']) for r in rows))
for v in all_v:
    print(f"{v:12.1f}", end="")
    for bp in bps:
        match = [r for r in rows if r['bp']==bp and abs(float(r['v_km_s'])-v)<1]
        if match:
            print(f"  {abs(float(match[0]['deviation_pct'])):8.2f}%", end="")
        else:
            print(f"  {'N/A':>8}", end="")
    print()
