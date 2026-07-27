# O0 — Unified analytical framework (§2.5)

Cross-cutting statistical framework that makes MANUSCRIPT §2.5 true. Script: `scripts/obj23/153_nested_framework.py`.
Plot level, sample_coverage ≥ 0.90. Model = OLS with **NEON domain fixed effects** + **site-clustered robust SEs**.

> **Why not a mixed model?** A site random intercept was specified but its variance collapsed to ~0 (singular)
> once domain fixed effects + canopy predictors were included — many NEON domains contain a single site, so the
> domain-fixed + site-random design is near-singular. OLS + domain + site-clustered SE gives identical fixed-effect
> inference with well-defined likelihoods (finite AIC/LRT). We report **R² beyond domain** (R²_full − R²_domain).
>
> **Response transform:** Hill numbers are log-transformed (satisfies homoscedasticity/normality: Breusch–Pagan
> p>0.4, Shapiro p>0.2); LCBD indices remain skewed under transformation and use the heteroscedasticity-consistent
> site-clustered SEs. Nested conclusions identical on raw and log scales.
> **VIF:** iterative to <5 (max 4.0). Height & gap metrics dropped (collinear with retained density/heterogeneity).

## Files
- `vif_screening.csv` / `vif_dropped.csv` — VIF<5 retained set (max VIF 4.0) & pruning log. Retained: structure {Rugosity, Vert_CV, VCI, LAI}, spectral {EVI}, dynamics {7 trends}, productivity {MODIS, PML}.
- `nested_models.csv` — M0→M1→M2→M3 R²/AIC/BIC/llf + LRT rows (ΔR², ΔAIC, χ², df, p) for M1→M2 (dynamics) and M2→M3 (function).
- `variance_partition.csv` — unique fractions (structure/spectral/dynamics/productivity) + shared, beyond domain.
- `coeffs_m3_clustered.csv` — full-model (M3) standardized betas with site-clustered SEs and p.
- `morans_I.csv` — Moran's I on M3 residuals (kNN k=8, 999 perms) + block uniques after adding a spatial polynomial trend.
- `interaction_fdr.csv` — O4 interaction p-values with Benjamini–Hochberg q (DHI terms excluded).

## Headline results (log-Hill; VIF<5; PRIMARY test = cluster-robust joint Wald, `cluster_wald.csv`)
> LRT/AIC (`nested_models.csv`) are secondary and more liberal; cluster-robust Wald accounts for 19-site clustering.
- **Dynamics beyond state:** Hill q1 p=2e-4, Hill q2 p=1e-4, nestedness p=0.012 — significant; **turnover p=0.61 (n.s.)**.
- **Function/GPP beyond state+dynamics:** **turnover p=4e-4** (strong); Hill q1 p=0.06, Hill q2 p=0.07, nestedness p=0.07 — marginal/n.s.
- Coefficients: fully-standardized in `coeffs_fullstd.csv` (VCI strongest for alpha; turnover-GPP product-inconsistent PML +0.42 / MODIS −0.25).
- **Variance partition (beyond domain):** alpha unique = structure(0.106) > dynamics(0.069) > spectral(0.011) > productivity(0.008); **spectral unique ≈ 0 for beta** (turnover 0.005, nestedness 0.001) — spectral tracks alpha only; dynamics contributes strongly to beta (nestedness 0.071).
- **Moran's I:** M3 residuals show significant positive spatial autocorrelation (I=0.05–0.15, p≤0.019), but block unique contributions are essentially unchanged after adding a spatial trend → inferences robust.
- **Interaction FDR:** 20/68 interactions survive Benjamini–Hochberg (mostly stand-age × spectral/structure and disturbance × Deep_Gap/SAVI/EVI).

## Note vs original proposal
Realized OLS+domain+cluster (not mixed) and R²-beyond-domain (not Nakagawa marginal R²) for the reasons above; GPP entered M3 as plot-level MODIS+PML (both VIF<5), with tower GPP retained as the site-level species–energy validation (O3 figure).
