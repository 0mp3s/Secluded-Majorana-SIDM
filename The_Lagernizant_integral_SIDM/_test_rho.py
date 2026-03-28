import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from lagrangian_path_integral import solve_friedmann_sigma, M_PL_GEV

params = {
    "Lambda_d": 2.28e-12,
    "f_gev": 0.24 * M_PL_GEV,
    "theta_i": 2.0,
    "omega_chi_h2": 0.1200,
    "m_chi_gev": 30.0,
    "alpha": 0.03,
    "n_points": 2000,
}

for ti in [1.5, 2.0, 2.5, 3.0]:
    params["theta_i"] = ti
    r = solve_friedmann_sigma(**params)
    h0 = r["H0_km_s_Mpc"]
    conv = r["converged"]
    vf = r["V_sigma_final_GeV4"]
    print(f"  ti={ti}  H0={h0:.2f}  converged={conv}  V_final={vf:.3e}")
