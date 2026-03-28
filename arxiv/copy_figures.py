"""
arxiv/copy_figures.py
=====================
Copy all figures referenced in main.tex from their module output
directories into arxiv/figures/.  Python equivalent of copy_figures.ps1,
compatible with the pipeline runner.
"""
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEST = Path(__file__).resolve().parent / "figures"

FIGURES = {
    "v31_island_of_viability.png":    "relic_density/output/v31_island_of_viability.png",
    "v31_bp1_velocity_profile.png":   "relic_density/output/v31_bp1_velocity_profile.png",
    "v33_observational_comparison.png": "observations/output/v33_observational_comparison.png",
    "v34_chi2_fit.png":               "observations/output/v34_chi2_fit.png",
    "v38_corner.png":                 "stats_mcmc/output/v38_corner.png",
    "v38_lambda_posterior.png":       "stats_mcmc/output/v38_lambda_posterior.png",
    "v38_autocorr_diagnostic.png":    "stats_mcmc/output/v38_autocorr_diagnostic.png",
    "grid_convergence.png":           "vpm_scan/output/grid_convergence.png",
    "v37_velocity_averaged.png":      "cross_checks/output/v37_velocity_averaged.png",
}


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    ok = 0
    for name, rel_path in FIGURES.items():
        src = REPO_ROOT / rel_path
        dst = DEST / name
        if src.exists():
            shutil.copy2(src, dst)
            print(f"[OK]      {name}")
            ok += 1
        else:
            print(f"[MISSING] {rel_path}")

    total = len(list(DEST.glob("*")))
    print(f"\nDone. {ok}/{len(FIGURES)} copied, {total} figures in arxiv/figures/")


if __name__ == "__main__":
    main()
