#!/usr/bin/env python3
"""
V8 - v22_raw_scan.py
====================
SIDM parameter scan storing ALL viable alpha values (raw output),
plus a compressed representative set (1 per cell).

Changes from V7/v18:
  - Stores EVERY (m_chi, m_phi, alpha) that passes SIDM cuts
  - Outputs two CSVs:
      all_viable_raw_v8.csv          - every viable sample
      all_viable_representative_v8.csv - 1 per (m_phi, resonance) cell
  - Reports both counts separately
  - Everything else identical: same VPM solver, same grid, same cuts

VPM ODE (attractive Yukawa V = -alpha e^{-m_phi r}/r):
  d(delta_l)/dx = +(lam e^{-x} / (kappa*x)) [j_hat_l(kappa*x) cos(delta) - n_hat_l(kappa*x) sin(delta)]^2
  lam = alpha * m_chi / m_phi > 0
"""
import sys, os, math, time
import numpy as np
from numba import jit, prange
from scipy.special import spherical_jn, spherical_yn
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from output_manager import timestamped_path
from run_logger import RunLogger
from global_config import GC

if sys.stdout.encoding != 'utf-8':
    sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', buffering=1)
    sys.stderr = open(sys.stderr.fileno(), mode='w', encoding='utf-8', buffering=1)
os.environ['PYTHONIOENCODING'] = 'utf-8'

# ==============================================================
#  Constants (sourced from global_config.json)
# ==============================================================
_PC = GC.physical_constants()
_CC = GC.cosmological_constants()

GEV2_TO_CM2 = _PC["GEV2_to_cm2"]
GEV_IN_G    = _PC["GeV_in_g"]
C_KM_S      = _PC["c_km_s"]
C_CM_S      = _PC["c_cm_s"]
M_PL        = round(_PC["m_pl_GeV"], -17)   # 1.22e19 — 2-sig-fig precision (validation section only)
M_E         = _PC["m_e_GeV"]
GEV_TO_SEC  = _PC["hbar_GeV_s"]
HBAR        = _PC["hbar_GeV_s"]
ALPHA_EM    = _PC["alpha_em"]
G_STAR      = _CC["g_star_s_20_90_GeV"]


# ==============================================================
#  Spherical Bessel functions (Numba)
# ==============================================================

@jit(nopython=True, cache=True)
def sph_jn_numba(l, z):
    """Spherical Bessel function j_l(z) via upward recurrence (Numba JIT)."""
    if z < 1e-30:
        return 1.0 if l == 0 else 0.0
    j0 = math.sin(z) / z
    if l == 0:
        return j0
    j1 = math.sin(z) / (z * z) - math.cos(z) / z
    if l == 1:
        return j1
    j_prev, j_curr = j0, j1
    for n in range(1, l):
        j_next = (2 * n + 1) / z * j_curr - j_prev
        j_prev = j_curr
        j_curr = j_next
        if abs(j_curr) < 1e-300:
            return 0.0
    return j_curr


@jit(nopython=True, cache=True)
def sph_yn_numba(l, z):
    """Spherical Bessel function y_l(z) via upward recurrence (Numba JIT)."""
    if z < 1e-30:
        return -1e300
    y0 = -math.cos(z) / z
    if l == 0:
        return y0
    y1 = -math.cos(z) / (z * z) - math.sin(z) / z
    if l == 1:
        return y1
    y_prev, y_curr = y0, y1
    for n in range(1, l):
        y_next = (2 * n + 1) / z * y_curr - y_prev
        y_prev = y_curr
        y_curr = y_next
        if abs(y_curr) > 1e200:
            return y_curr
    return y_curr


# ==============================================================
#  VPM phase shift solver (identical to V7/v18)
# ==============================================================

@jit(nopython=True, cache=True)
def _vpm_rhs(l, kappa, lam, x, delta):
    """RHS of the VPM ODE: d(delta_l)/dx for Yukawa potential."""
    if x < 1e-20:
        return 0.0
    z = kappa * x
    if z < 1e-20:
        return 0.0
    jl = sph_jn_numba(l, z)
    nl = sph_yn_numba(l, z)
    j_hat = z * jl
    n_hat = -z * nl
    cd = math.cos(delta)
    sd = math.sin(delta)
    bracket = j_hat * cd - n_hat * sd
    if not math.isfinite(bracket):
        return 0.0
    pot = lam * math.exp(-x) / (kappa * x)
    val = pot * bracket * bracket
    return val if math.isfinite(val) else 0.0


@jit(nopython=True, cache=True)
def vpm_phase_shift(l, kappa, lam, x_max=50.0, N_steps=4000):
    """Compute partial-wave phase shift delta_l via RK4 integration of the VPM ODE."""
    if lam < 1e-30 or kappa < 1e-30:
        return 0.0
    x_min = max(1e-5, 0.05 / (kappa + 0.01))
    if l > 0:
        x_barrier = l / kappa
        if x_barrier > x_min:
            x_min = x_barrier
    h = (x_max - x_min) / N_steps
    delta = 0.0
    for i in range(N_steps):
        x = x_min + i * h
        k1 = _vpm_rhs(l, kappa, lam, x, delta)
        k2 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k1)
        k3 = _vpm_rhs(l, kappa, lam, x + 0.5*h, delta + 0.5*h*k2)
        k4 = _vpm_rhs(l, kappa, lam, x + h, delta + h*k3)
        delta += h * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return delta


# ==============================================================
#  sigma_T for identical Majorana fermions
# ==============================================================

@jit(nopython=True, cache=True)
def sigma_T_vpm(m_chi, m_phi, alpha, v_km_s):
    """Cross section sigma_T/m [cm^2/g] for identical Majorana fermions.

    For identical Majorana fermions sigma_elastic = sigma_transfer exactly
    (each spin channel has definite parity so <cos theta> = 0).
    Sums partial waves with weights w=1 (even l) and w=3 (odd l).
    Args: m_chi [GeV], m_phi [GeV], alpha (coupling), v_km_s [km/s].
    """
    v = v_km_s / C_KM_S
    mu = m_chi / 2.0
    k = mu * v
    kappa = k / m_phi
    lam = alpha * m_chi / m_phi
    if kappa < 1e-15:
        return 0.0

    if kappa < 5:
        x_max, N_steps = 50.0, 4000
    elif kappa < 50:
        x_max, N_steps = 80.0, 8000
    else:
        x_max, N_steps = 100.0, 12000

    l_max = min(max(3, int(kappa) + 3), 80)

    sigma_sum = 0.0
    for l in range(l_max + 1):
        delta = vpm_phase_shift(l, kappa, lam, x_max, N_steps)
        weight = 1.0 if l % 2 == 0 else 3.0
        contrib = weight * (2*l + 1) * math.sin(delta)**2
        sigma_sum += contrib
        if l > int(kappa) + 1 and sigma_sum > 0:
            if contrib / sigma_sum < 1e-3:
                break

    sigma_GeV2 = 2.0 * math.pi * sigma_sum / (k * k)
    sigma_cm2 = sigma_GeV2 * GEV2_TO_CM2
    return sigma_cm2 / (m_chi * GEV_IN_G)


# ==============================================================
#  Scan grid  — configurable via --config or vpm_scan/config.json
# ==============================================================
from config_loader import load_config as _load_config
_SCAN_CFG = _load_config(__file__).get("grid", {})

N_CHI   = _SCAN_CFG.get("n_chi", 50)
N_PHI   = _SCAN_CFG.get("n_phi", 70)
N_RES   = _SCAN_CFG.get("n_res", 4)
N_ALPHA = _SCAN_CFG.get("n_alpha", 200)

_m_chi_range = _SCAN_CFG.get("m_chi_range_GeV", [0.1, 100.0])
_m_phi_range = _SCAN_CFG.get("m_phi_range_GeV", [0.1e-3, 200e-3])

M_CHI_VALS = np.logspace(np.log10(_m_chi_range[0]), np.log10(_m_chi_range[1]), N_CHI)
M_PHI_VALS = np.logspace(np.log10(_m_phi_range[0]), np.log10(_m_phi_range[1]), N_PHI)
LAM_CRITS  = np.array(_SCAN_CFG.get("lambda_crits", [1.68, 6.45, 14.7, 26.0]))


# ==============================================================
#  RAW scan — stores ALL viable alphas per (m_phi, resonance)
# ==============================================================

@jit(nopython=True, cache=True)
def scan_one_mchi_raw(mc, m_phi_arr, lam_crits, n_phi, n_res, n_alpha,
                      sd_lo, sd_hi, sc_hi):
    """Scan all (m_phi, resonance, alpha) for one m_chi.
    Returns arrays of ALL viable samples (not compressed).

    SIDM viability cuts (from global_config.json, sourced from literature):
        sd_lo  — sigma/m(30 km/s) lower bound  [KTY16]
        sd_hi  — sigma/m(30 km/s) upper bound  [KTY16]
        sc_hi  — sigma/m(1000 km/s) upper bound [Harvey+2015, conservative]

    Uses a pre-allocated buffer. Returns count of valid entries.
    Max possible hits per m_chi = n_phi * n_res * n_alpha.
    """
    max_buf = n_phi * n_res * n_alpha
    buf_mphi  = np.zeros(max_buf)
    buf_alpha = np.zeros(max_buf)
    buf_sd    = np.zeros(max_buf)
    buf_sc    = np.zeros(max_buf)
    buf_ires  = np.zeros(max_buf, dtype=np.int32)
    count = 0

    # Also track representative (1 per cell, closest to sigma=5)
    rep_max = n_phi * n_res
    rep_mphi  = np.zeros(rep_max)
    rep_alpha = np.zeros(rep_max)
    rep_sd    = np.zeros(rep_max)
    rep_sc    = np.zeros(rep_max)
    rep_valid = np.zeros(rep_max, dtype=np.int32)
    rep_dist  = np.full(rep_max, 1e10)

    for ip in range(n_phi):
        mp = m_phi_arr[ip]
        for ir in range(n_res):
            lam_c = lam_crits[ir]
            alpha_c = lam_c * mp / mc
            if alpha_c > 0.5 or alpha_c < 1e-7:
                continue
            a_lo = alpha_c * 0.7
            a_hi = alpha_c * 1.3
            rep_idx = ip * n_res + ir

            for ia in range(n_alpha):
                frac = ia / (n_alpha - 1.0)
                alpha = a_lo * (a_hi / a_lo) ** frac

                sd = sigma_T_vpm(mc, mp, alpha, 30.0)
                if sd < sd_lo or sd > sd_hi:
                    continue
                sc = sigma_T_vpm(mc, mp, alpha, 1000.0)
                if sc >= sc_hi:
                    continue

                # Store raw hit
                if count < max_buf:
                    buf_mphi[count]  = mp
                    buf_alpha[count] = alpha
                    buf_sd[count]    = sd
                    buf_sc[count]    = sc
                    buf_ires[count]  = ir
                    count += 1

                # Update representative (closest to sigma=5)
                dist = abs(sd - 5.0)
                if dist < rep_dist[rep_idx]:
                    rep_dist[rep_idx]  = dist
                    rep_mphi[rep_idx]  = mp
                    rep_alpha[rep_idx] = alpha
                    rep_sd[rep_idx]    = sd
                    rep_sc[rep_idx]    = sc
                    rep_valid[rep_idx] = 1

    return (buf_mphi[:count], buf_alpha[:count], buf_sd[:count],
            buf_sc[:count], buf_ires[:count], count,
            rep_mphi, rep_alpha, rep_sd, rep_sc, rep_valid)


# ==============================================================
#  Scipy reference (for validation)
# ==============================================================

def _scipy_ref(l, kappa, lam, x_max=50.0):
    from scipy.integrate import solve_ivp
    x_min = max(1e-6, 0.05 / (kappa + 0.01))
    if l > 0:
        x_barrier = l / kappa
        if x_barrier > x_min:
            x_min = x_barrier
    def rhs(x, y):
        z = kappa * x
        if z < 1e-30: return [0.0]
        jl = spherical_jn(l, z)
        nl = spherical_yn(l, z)
        j_hat, n_hat = z * jl, -z * nl
        pot = lam * math.exp(-x) / (kappa * x)
        cd, sd = math.cos(y[0]), math.sin(y[0])
        b = j_hat * cd - n_hat * sd
        if not math.isfinite(b): return [0.0]
        v = pot * b * b
        return [v] if math.isfinite(v) else [0.0]
    sol = solve_ivp(rhs, [x_min, x_max], [0.0], method='RK45',
                    rtol=1e-10, atol=1e-13, max_step=0.5)
    return sol.y[0, -1] if sol.success else float('nan')


# ==============================================================
#  Validation
# ==============================================================

def section_0_validation():
    print("=" * 75)
    print("SECTION 0: VPM Solver Validation")
    print("=" * 75)
    print()

    # 0a: Free particle
    print("  [0a] Free particle (lam=0)")
    all_zero = True
    for l in [0, 1, 2, 5, 10, 20]:
        d = vpm_phase_shift(l, 2.0, 0.0)
        ok = abs(d) < 1e-10
        all_zero = all_zero and ok
        tag = "PASS" if ok else "FAIL"
        print(f"    l={l:3d}: delta = {d:.6e}  [{tag}]")
    print(f"  [0a] {'PASS' if all_zero else 'FAIL'}")
    print()

    # 0b: sigma_T_vpm accuracy vs scipy
    print("  [0b] sigma_T_vpm accuracy vs scipy ground truth")
    mc_t, mp_t, a_t = 75.0, 2.0e-3, 5e-5
    mu_t = mc_t / 2.0
    v_t = 30.0 / C_KM_S
    k_t = mu_t * v_t
    kap_t = k_t / mp_t
    lam_t = a_t * mc_t / mp_t
    l_ref = int(kap_t) + 6
    print(f"    kappa={kap_t:.3f}, lam={lam_t:.3f}")
    s_ref = 0.0
    for l in range(l_ref + 1):
        d = _scipy_ref(l, kap_t, lam_t)
        w = 1.0 if l % 2 == 0 else 3.0
        c = w * (2*l+1) * math.sin(d)**2
        s_ref += c
    sig_ref = 2*math.pi * s_ref / k_t**2 * GEV2_TO_CM2 / (mc_t * GEV_IN_G)
    sig_numba = sigma_T_vpm(mc_t, mp_t, a_t, 30.0)
    err_pct = abs(sig_numba - sig_ref) / sig_ref * 100
    print(f"    scipy sigma/m = {sig_ref:.6f}")
    print(f"    Numba sigma/m = {sig_numba:.6f}")
    print(f"    agreement = {err_pct:.2f}%")
    converged = err_pct < 3.0
    print(f"  [0b] {'PASS' if converged else 'FAIL'} (<3% tolerance)")
    print()

    overall = all_zero and converged
    print(f"  OVERALL: {'ALL PASSED' if overall else 'SOME FAILED'}")
    return overall


# ==============================================================
#  Main Scan
# ==============================================================

def section_1_scan():
    print()
    print("=" * 75)
    print("SECTION 1: Raw Parameter Scan (all viable alphas)")
    print("=" * 75)
    total = N_CHI * N_PHI * N_RES * N_ALPHA
    n_workers = min(14, multiprocessing.cpu_count())
    _cuts = GC.sidm_cuts()
    print(f"  Grid: {N_CHI} x {N_PHI} x {N_RES} x {N_ALPHA} = {total:,}")
    print(f"  SIDM (raw scan): sigma/m(30) in [{_cuts['sigma_m_30_lo']},{_cuts['sigma_m_30_hi']}], sigma/m(1000) < {_cuts['sigma_m_1000_hi']}  [from global_config.json]")
    print(f"  Parallel workers: {n_workers}\n")

    all_raw = []
    all_rep = []
    t0 = time.time()
    done = 0

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {}
        for ic in range(N_CHI):
            mc = float(M_CHI_VALS[ic])
            f = pool.submit(_scan_worker, ic, mc)
            futures[f] = (ic, mc)

        for f in as_completed(futures):
            ic, mc = futures[f]
            raw_list, rep_list = f.result()
            all_raw.extend(raw_list)
            all_rep.extend(rep_list)
            done += 1
            elapsed = time.time() - t0
            eta = elapsed / done * (N_CHI - done) if done > 0 else 0
            hs = f"  +{len(raw_list)} raw, +{len(rep_list)} rep" if raw_list else ""
            print(f"  [{done:2d}/{N_CHI}] m_chi={mc:8.3f} GeV  "
                  f"({elapsed:6.1f}s, ETA {eta:5.0f}s)  "
                  f"raw: {len(all_raw)}, rep: {len(all_rep)}{hs}", flush=True)

    # Sort by m_chi for consistent output
    all_raw.sort(key=lambda p: (p['m_chi'], p['m_phi']))
    all_rep.sort(key=lambda p: (p['m_chi'], p['m_phi']))

    print(f"\n  Done: {time.time()-t0:.1f}s")
    print(f"  Raw viable samples: {len(all_raw)}")
    print(f"  Representative cells: {len(all_rep)}")
    return all_raw, all_rep


def _scan_worker(ic, mc):
    """Worker function for parallel scan. Called per m_chi value."""
    # Numba cache=True means JIT is only compiled once on disk,
    # each subprocess loads the cached version automatically.
    _cuts = GC.sidm_cuts()
    (raw_mphi, raw_alpha, raw_sd, raw_sc, raw_ires, raw_count,
     rep_mphi, rep_alpha, rep_sd, rep_sc, rep_valid) = \
        scan_one_mchi_raw(mc, M_PHI_VALS, LAM_CRITS,
                          N_PHI, N_RES, N_ALPHA,
                          _cuts["sigma_m_30_lo"],
                          _cuts["sigma_m_30_hi"],
                          _cuts["sigma_m_1000_hi"])

    raw_list = []
    for j in range(raw_count):
        raw_list.append({
            'm_chi': mc, 'm_phi': float(raw_mphi[j]),
            'alpha': float(raw_alpha[j]),
            'sigma_30': float(raw_sd[j]),
            'sigma_1000': float(raw_sc[j]),
            'resonance': int(raw_ires[j])
        })

    rep_list = []
    for j in range(len(rep_valid)):
        if rep_valid[j] > 0:
            rep_list.append({
                'm_chi': mc, 'm_phi': float(rep_mphi[j]),
                'alpha': float(rep_alpha[j]),
                'sigma_30': float(rep_sd[j]),
                'sigma_1000': float(rep_sc[j])
            })

    return raw_list, rep_list


# ==============================================================
#  Analysis & Output
# ==============================================================

def section_2_analysis(raw_points, rep_points):
    print()
    print("=" * 75)
    print("SECTION 2: Results")
    print("=" * 75)
    print()

    if not raw_points:
        print("  No viable points found.")
        return None

    # --- Raw analysis ---
    n_raw = len(raw_points)
    mc_raw = np.array([p['m_chi'] for p in raw_points])
    mp_raw = np.array([p['m_phi'] for p in raw_points]) * 1000  # MeV
    al_raw = np.array([p['alpha'] for p in raw_points])
    sd_raw = np.array([p['sigma_30'] for p in raw_points])
    sc_raw = np.array([p['sigma_1000'] for p in raw_points])
    lam_raw = al_raw * mc_raw / (mp_raw / 1000)
    ratio_raw = sd_raw / np.maximum(sc_raw, 1e-20)

    print(f"  === RAW VIABLE SAMPLES: {n_raw} ===")
    print(f"    m_chi : {mc_raw.min():.3f} - {mc_raw.max():.3f} GeV")
    print(f"    m_phi : {mp_raw.min():.3f} - {mp_raw.max():.3f} MeV")
    print(f"    alpha : {al_raw.min():.2e} - {al_raw.max():.2e}")
    print(f"    lam   : {lam_raw.min():.2f} - {lam_raw.max():.2f}")
    print(f"    sigma/m(30): {sd_raw.min():.2f} - {sd_raw.max():.2f}")
    print(f"    sigma/m(1k): {sc_raw.min():.6f} - {sc_raw.max():.6f}")
    print(f"    ratio: {ratio_raw.min():.0f} - {ratio_raw.max():.0f}x")
    print()

    bbn_raw = mp_raw > 2 * M_E * 1000
    print(f"  BBN filter (m_phi > 1.022 MeV): {int(bbn_raw.sum())}/{n_raw}")
    print()

    # --- Representative analysis ---
    n_rep = len(rep_points)
    mc_rep = np.array([p['m_chi'] for p in rep_points])
    mp_rep = np.array([p['m_phi'] for p in rep_points]) * 1000
    bbn_rep = mp_rep > 2 * M_E * 1000

    print(f"  === REPRESENTATIVE SET (1 per cell): {n_rep} ===")
    print(f"  BBN filter: {int(bbn_rep.sum())}/{n_rep}")
    print()

    # --- Per-mass distribution ---
    print(f"  {'m_chi':>10}  {'raw':>6}  {'raw_BBN':>7}  {'rep':>5}  {'rep_BBN':>7}")
    unique = np.unique(np.round(mc_raw, 2))
    for m in unique:
        mask_raw = np.abs(mc_raw - m) < 0.01
        mask_rep = np.abs(mc_rep - m) < 0.01
        n_r = int(mask_raw.sum())
        n_rb = int((mask_raw & bbn_raw).sum())
        n_p = int(mask_rep.sum())
        n_pb = int((mask_rep & bbn_rep).sum())
        print(f"  {m:10.3f}  {n_r:6d}  {n_rb:7d}  {n_p:5d}  {n_pb:7d}")
    print()

    # --- Compression ratio ---
    print(f"  Compression ratio: {n_raw} raw / {n_rep} representative = {n_raw/max(n_rep,1):.1f}x")
    print()

    # --- Best benchmark (BBN-safe representative, closest to sigma=5) ---
    dist_rep = np.array([abs(p['sigma_30'] - 5.0) for p in rep_points])
    dist_rep[~bbn_rep] = 1e10
    bi = int(np.argmin(dist_rep))
    bp = rep_points[bi]
    lb = bp['alpha'] * bp['m_chi'] / bp['m_phi']
    print(f"  Best benchmark (representative, BBN-safe):")
    print(f"    m_chi={bp['m_chi']:.6f} GeV, m_phi={bp['m_phi']*1000:.6f} MeV")
    print(f"    alpha={bp['alpha']:.6e}, lam={lb:.4f}")
    print(f"    sigma/m(30)={bp['sigma_30']:.4f}, sigma/m(1k)={bp['sigma_1000']:.6f}")
    print(f"    ratio={bp['sigma_30']/bp['sigma_1000']:.0f}x")

    # --- Save CSVs ---
    csv_raw = timestamped_path("all_viable_raw_v8")
    with open(csv_raw, 'w') as f:
        f.write("m_chi_GeV,m_phi_MeV,alpha,sigma_m_30,sigma_m_1000,resonance_idx\n")
        for p in raw_points:
            f.write(f"{p['m_chi']:.6f},{p['m_phi']*1e3:.10e},{p['alpha']:.10e},"
                    f"{p['sigma_30']:.6f},{p['sigma_1000']:.8f},{p['resonance']}\n")
    print(f"\n  Saved raw: {csv_raw}")

    csv_rep = timestamped_path("all_viable_representative_v8")
    with open(csv_rep, 'w') as f:
        f.write("m_chi_GeV,m_phi_MeV,alpha,sigma_m_30,sigma_m_1000\n")
        for p in rep_points:
            f.write(f"{p['m_chi']:.6f},{p['m_phi']*1e3:.10e},{p['alpha']:.10e},"
                    f"{p['sigma_30']:.6f},{p['sigma_1000']:.8f}\n")
    print(f"  Saved representative: {csv_rep}")

    return bp, csv_raw, csv_rep


# ==============================================================
#  Main
# ==============================================================

def main():
    print("=" * 75)
    print("V8 - v22_raw_scan.py")
    print("SIDM Parameter Scan: Raw + Representative Output")
    print(f"  CPU cores available: {multiprocessing.cpu_count()}")
    print("=" * 75)
    t_start = time.time()

    _cuts = GC.sidm_cuts()
    _rl_params = {
        "n_chi": N_CHI, "n_phi": N_PHI, "n_res": N_RES, "n_alpha": N_ALPHA,
        "sigma_m_30_lo": _cuts["sigma_m_30_lo"],
        "sigma_m_30_hi": _cuts["sigma_m_30_hi"],
        "sigma_m_1000_hi": _cuts["sigma_m_1000_hi"],
    }

    with RunLogger(
        script="core/v22_raw_scan.py",
        stage="0 - VPM Scan",
        params=_rl_params,
        data_source="grid scan from scratch",
    ) as rl:
        # JIT warmup — compiles Numba cache before spawning workers
        print("\n  Warming up Numba JIT...", flush=True)
        _ = sigma_T_vpm(10.0, 1e-3, 1e-4, 30.0)
        print("  JIT ready.\n")

        section_0_validation()
        raw, rep = section_1_scan()
        bp, csv_raw, csv_rep = section_2_analysis(raw, rep)

        rl.set_n_viable(len(raw))
        rl.add_output(str(csv_raw))
        rl.add_output(str(csv_rep))
        rl.set_notes(f"raw={len(raw)}, rep={len(rep)}, BBN-safe={sum(1 for p in rep if p['m_phi']*1e3 > 1.022)}")

    elapsed = time.time() - t_start
    print(f"\n  Total runtime: {elapsed:.1f}s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
