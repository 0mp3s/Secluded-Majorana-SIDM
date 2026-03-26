# Copy all figures referenced in main.tex into arxiv/figures/
# Run from the repo root:  .\arxiv\copy_figures.ps1

$root = Split-Path $PSScriptRoot -Parent
$dest = Join-Path $PSScriptRoot "figures"

$figures = @{
    "v31_island_of_viability.png"    = "relic_density\output\v31_island_of_viability.png"
    "v31_bp1_velocity_profile.png"   = "relic_density\output\v31_bp1_velocity_profile.png"
    "v33_observational_comparison.png"= "observations\output\v33_observational_comparison.png"
    "v34_chi2_fit.png"               = "observations\output\v34_chi2_fit.png"
    "v38_corner.png"                 = "stats_mcmc\output\v38_corner.png"
    "v38_lambda_posterior.png"       = "stats_mcmc\output\v38_lambda_posterior.png"
    "v38_autocorr_diagnostic.png"    = "stats_mcmc\output\v38_autocorr_diagnostic.png"
    "grid_convergence.png"           = "vpm_scan\output\grid_convergence.png"
    "v37_velocity_averaged.png"      = "cross_checks\output\v37_velocity_averaged.png"
}

foreach ($fig in $figures.GetEnumerator()) {
    $src = Join-Path $root $fig.Value
    $dst = Join-Path $dest $fig.Key
    if (Test-Path $src) {
        Copy-Item $src $dst -Force
        Write-Host "[OK] $($fig.Key)"
    } else {
        Write-Host "[MISSING] $($fig.Value)" -ForegroundColor Red
    }
}

Write-Host "`nDone. $((Get-ChildItem $dest -File).Count) figures in arxiv/figures/"
