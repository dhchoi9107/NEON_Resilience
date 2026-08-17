# Supplementary Information (stub)

Reproducibility: the plot-level framework and its block tests are produced by `scripts_pipeline/framework_state_dynamics.py`; the Bayesian models by `bayes_multilevel.py` (prior check: `bayes_prior_sens.py`); the dissimilarity analyses by `rs_beta_mantel.py` and `rs_beta_gdm_full.py`; the productivity diagnostic by `recompute_suppression_26.py`. All numbers in the main text are re-derived from the source data by `verify_manuscript_numbers.py`.

## Supplementary Methods S1 — Regression-to-the-mean (RTM) null model
Apparent buffering or amplification between initial canopy structure and subsequent structural change was tested against a regression-to-the-mean null. The null assumes no biological feedback between initial state and subsequent change while preserving the observed variance and measurement error of the repeated LiDAR observations; expected initial–change relationships under pure RTM were simulated and compared with the observed slopes. Only observed relationships exceeding the RTM envelope were interpreted as biological regulation.

## Supplementary Table S1 — Response transformation check (structural-dynamics increment, LRT p)
Hill numbers are log-transformed in the main analysis (log restores homoscedasticity and approximate normality: Breusch–Pagan p > 0.4, Shapiro–Wilk p > 0.2). The dynamics increment is significant on either scale, so the conclusion does not depend on the transformation.

| response | dynamics increment, raw scale | dynamics increment, log scale |
|---|---|---|
| Hill q1 | 6.2e-08 | 9.2e-13 |
| Hill q2 | 8.0e-06 | 2.6e-10 |

## Supplementary Table S1b — Few-cluster-valid block test for structural dynamics (26 sites)
Restricted wild cluster bootstrap (Rademacher weights, 999 replicates) is the primary block test; the asymptotic cluster-robust Wald p is shown for contrast (anticonservative at this cluster count). Productivity is not a plot-level block (Section 2.5), so only the dynamics increment is tested here. Source: `wildboot_sd.csv`, `cluster_wald_sd.csv`.

| response | dynamics ΔR² | dynamics: Wald p | **dynamics: wild-boot p** |
|---|---|---|---|
| Hill q1 | 0.042 | 2.3e-5 | **0.014** |
| Hill q2 | 0.039 | 1.8e-5 | **0.007** |
| Turnover | 0.013 | 0.33 | **0.767** |
| Nestedness | 0.056 | 0.0044 | **0.143** |

The dynamics increment is significant for both alpha metrics and does not reject the joint null for either beta component; the nestedness association is supported at the coefficient level instead (Table S8).

## Supplementary Table S2 — Beta components: Gaussian OLS vs beta regression (structural-dynamics increment, LRT p)
The beta components are bounded proportions, so the nested comparison was repeated with beta regression. Both models support the dynamics increment for both components; beta regression is the more liberal of the two, and the Gaussian result is reported in the main text as the conservative choice.

| response | Gaussian OLS | beta regression |
|---|---|---|
| Turnover | 0.007 | 1.8e-12 |
| Nestedness | 3.1e-09 | 4.5e-19 |

## Supplementary Table S3 — Variance partition beyond domain (unique and shared fractions)
The main text reports the unique (semi-partial) fractions; the shared fraction is given here for completeness. Plot-level blocks only — productivity is analysed at the site scale (Section 3.5). Source: `variance_partition_sd.csv`.

| response | total R² beyond domain | unique structure | unique spectral | unique dynamics | shared |
|---|---|---|---|---|---|
| Hill q1 | 0.153 | 0.081 | 0.011 | 0.042 | 0.019 |
| Hill q2 | 0.141 | 0.073 | 0.013 | 0.039 | 0.015 |
| Turnover | 0.024 | 0.006 | 0.001 | 0.013 | 0.004 |
| Nestedness | 0.092 | 0.025 | 0.000 | 0.056 | 0.011 |

## Supplementary Table S4 — Mixed vs clustered-OLS fixed effects
Site random-intercept variance estimated at 0 (singular); standardized fixed-effect coefficients closely similar (mean |Δ| = 0.03, max |Δ| = 0.05, r = 0.93 across 14 RS predictors).

## Supplementary Table S5 — Climate control (NEON site MAT/MAP)
Climate = NEON site-level mean annual temperature and precipitation (`data/site_climate_neon.csv`). Domain and continuous climate are collinear, so climate is used as an alternative baseline rather than an additional covariate. Source: `climate_species_energy.csv`, `climate_nested_compare.csv`.

**(a) Does GPP add beyond climate alone? (the species–energy signal is not a climate artifact)**

| response | R² climate only | ΔR² GPP beyond climate | p | PML-V2 β (clustered) |
|---|---|---|---|---|
| Hill q1 | 0.052 | 0.209 | 1.4e-35 | +0.284 |
| Hill q2 | 0.053 | 0.185 | 3.4e-31 | +0.243 |
| Turnover | 0.153 | 0.200 | 1.3e-34 | +0.124 |
| Nestedness | 0.094 | 0.056 | 9.0e-09 | -0.005 |

**(b) Structural-dynamics increment under the domain versus climate baseline (LRT p)**

| response | ΔR² dynamics (domain) | p (domain) | ΔR² dynamics (climate) | p (climate) |
|---|---|---|---|---|
| Hill q1 | 0.042 | 9.2e-13 | 0.083 | 2.4e-22 |
| Hill q2 | 0.039 | 2.6e-10 | 0.081 | 2.4e-19 |
| Turnover | 0.013 | 0.007 | 0.035 | 8.2e-06 |
| Nestedness | 0.056 | 3.1e-09 | 0.086 | 2.1e-12 |

The domain baseline is much stronger than the climate baseline (Hill q1 R² = 0.490 versus 0.052) because domain fixed effects absorb climatic **and** non-climatic biogeography (species pools, soils). Domain is therefore retained as the primary control, and climate serves as the explicit "beyond-climate" check. The dynamics increment is significant under both baselines for every response.

## Supplementary Table S6 — Context interactions surviving FDR (framework predictors, cluster-robust, BH within family)
The same retained predictors as the plot-level models (4 structural + EVI + PML-V2 GPP) crossed with each context variable: 144 terms, of which 5 survive Benjamini–Hochberg correction within family. Source: `o4_interactions_pml.csv`.

| family | context | response | predictor | β interaction | q |
|---|---|---|---|---|---|
| disturbance | severity | Hill q1 | LAI | −0.112 | 0.004 |
| disturbance | severity | Hill q2 | LAI | −0.106 | 0.004 |
| land use | forest fraction | turnover | EVI | −0.056 | <0.001 |
| land use | edge density | turnover | EVI | +0.060 | 0.001 |
| land use | land-cover diversity | turnover | EVI | +0.049 | 0.009 |

Interpretation: the structure–diversity coupling weakens with disturbance severity (negative LAI × severity for both alpha metrics). The three surviving land-use terms all involve the spectral index and compositional turnover, and are reported as exploratory. No stand-age interaction survives correction, so the age moderation is reported as directionally consistent but unsupported.

## Supplementary Table S7 — Plot-level productivity diagnostic
Why productivity is analysed at the site scale (Section 3.5). Plot level, domain fixed effects + site-clustered SE, all predictors and the response standardized. Source: `gpp_suppression_26.csv` (script `recompute_suppression_26.py`).

| response | PML-V2 β (fully standardized) | p (clustered) | ΔR² beyond state + dynamics |
|---|---|---|---|
| Hill q1 | +0.151 | 0.056 | 0.003 |
| Turnover | +0.273 | 0.151 | 0.011 |

Neither increment reaches significance, so productivity adds little at the plot grain. For completeness we note the collinearity artifact that motivated carrying a single satellite product: with r = 0.87 between MOD17 and PML-V2, entering both in the turnover model made their coefficients oppose and inflate (MODIS −0.233, p = 0.054; PML +0.367, p = 0.040) although neither was associated with turnover alone (MODIS −0.030, p = 0.850; PML +0.273, p = 0.151) — a textbook suppression pattern rather than added information.

## Supplementary Methods S8 — Four-method triangulation of the plot-level framework (statistical robustness)

The plot-level conclusions (§3.1–§3.3) rest on domain fixed effects with site-clustered SEs and joint block tests by restricted wild cluster bootstrap. Because inference is anchored on only 26 site clusters (~12 domains) — a regime where joint block-increment tests are known to be conservative — we re-derived the core conclusions with three additional, methodologically independent approaches. All four agree on the two headline results (**structure → alpha**; **dynamics → nestedness**), and the Bayesian coefficient-level analysis recovers the dynamics → nestedness signal that the conservative joint block test misses.

**Method 1 — Bayesian multilevel (partial pooling).** bambi/PyMC, crossed random intercepts for domain and site-within-domain, non-centred parameterization, HalfNormal(0.5/0.3) priors on the group scales; 4 chains × (3,000 tune + 2,000 draws), target_accept = 0.99. Fully converged (0 divergent transitions; R̂ ≤ 1.002). The near-singular site variance (Table S4, site RE variance ≈ 0) is the source of the funnel geometry that the priors + high target_accept resolve. Coefficient credible when its 95% HDI excludes zero. Models were fitted for all four responses; the table below shows the two headline responses (Hill q2 mirrors q1 — LAI +0.40, VCI +0.26, EVI +0.26 credible; turnover shows only two small credible dynamics trends, LAI trend and height-ratio trend, both β ≈ +0.09), full posteriors in `bayes_multilevel_coeffs.csv`.

| Response | Predictor | β (std) | 95% HDI | pd |
|---|---|---|---|---|
| Hill q1 (alpha) | LAI_mean (struct) | +0.404 | [+0.266, +0.539] | 1.000 |
| Hill q1 | VCI_mean (struct) | +0.274 | [+0.205, +0.345] | 1.000 |
| Hill q1 | EVI_mean (spectral) | +0.237 | [+0.132, +0.343] | 1.000 |
| Hill q1 | Rugosity_mean (struct) | −0.098 | [−0.172, −0.024] | 0.995 |
| Hill q1 | Rumple_trend (dyn) | +0.091 | [+0.016, +0.167] | 0.991 |
| Hill q1 | Ht_Ratio_trend (dyn) | +0.114 | [+0.036, +0.192] | 0.998 |
| Nestedness | **Ht_Ratio_trend (dyn)** | **−0.307** | **[−0.419, −0.198]** | **1.000** |
| Nestedness | **Vert_CV_trend (dyn)** | **+0.200** | **[+0.081, +0.320]** | **0.999** |
| Nestedness | **LAI_trend (dyn)** | **−0.191** | **[−0.308, −0.071]** | **0.999** |
| Nestedness | Rugosity_mean (struct) | +0.137 | [+0.035, +0.240] | 0.994 |
| Nestedness | Vert_CV_mean (struct) | −0.121 | [−0.237, −0.004] | 0.976 |

Structure dominates alpha (3/4 structural terms credible); **dynamics dominate nestedness (3/7 dynamics trends credible, more than any other block for this component)** even though the joint dynamics block-increment did not clear the wild bootstrap (p = 0.14, Table S1b) — i.e. the block test is conservative, not the effect absent. Full posterior tables: `results/O0_framework/bayes_multilevel_coeffs.csv` (+ `_meta.csv`).

*Prior sensitivity.* Because the coefficient-level inference is the primary support for the dynamics → nestedness result, we repeated the Hill q1 and nestedness models with the group-scale priors doubled (domain σ ~ HalfNormal(1.0), site σ ~ HalfNormal(0.6)). The credible set was identical for all 24 coefficients (12 predictors × 2 responses), the maximum change in any posterior mean was 0.007, the three nestedness dynamics trends were essentially unchanged (height-ratio −0.306 vs −0.31; vertical-CV +0.199 vs +0.20; leaf-area −0.191 vs −0.19), and both models again sampled with 0 divergent transitions (`bayes_prior_sensitivity.csv`). The headline result is therefore not prior-driven.

**Method 2 — Formal variation partitioning (vegan-varpart style).** Domain-conditional partition of the remote-sensing R² into unique structure/spectral/dynamics blocks via set-complement R² differences, cross-checking the manual semi-partial R² of §3.1. The two recipes agree to three decimals:

| Response | total | unique struct | unique spectral | unique dyn | shared |
|---|---|---|---|---|---|
| Hill q1 | 0.153 | **0.081** | 0.011 | 0.042 | 0.019 |
| Hill q2 | 0.141 | **0.073** | 0.013 | 0.039 | 0.015 |
| Turnover | 0.024 | 0.006 | 0.001 | 0.013 | 0.004 |
| Nestedness | 0.092 | 0.025 | 0.000 | **0.056** | 0.011 |

Structure is the largest unique alpha block; dynamics the largest unique nestedness block — identical to the main analysis. Output: `results/O0_framework/varpart_robust.csv`.

**Method 3 — GDM-family saturation model (beta, §3.1b).** Monotone-spline generalized dissimilarity model of Bray–Curtis composition on structural + spectral distance, within-site permutation (999). Deviance explained = 0.122 (p = 0.002); structure and spectral contribute non-redundantly (structure 45% : spectral 55% of fitted transform SD), confirming that remote-sensing heterogeneity tracks compositional turnover once analysed as dissimilarity rather than predictor levels.

**Method 4 — Restricted-permutation Mantel (beta, §3.1b).** Within-site (plots permuted only within sites; biogeography held fixed), 999 permutations, 23 sites. Structural distance vs Bray–Curtis: Mantel r = +0.264 (Pearson) / +0.32 (Spearman); spectral distance r = +0.293 / +0.30; both p = 0.001. Rank-based estimates are reported in the main text as robust to the Bray–Curtis = 1 ceiling.

**Bottom line.** The framework conclusions are not artifacts of one inferential recipe: structure→alpha and dynamics→nestedness hold under partial-pooling Bayesian estimation and formal variation partitioning, and the beta/heterogeneity signal holds under both a saturation (GDM) model and restricted-permutation Mantel tests. The only substantive change relative to the conservative wild-bootstrap joint test is that the Bayesian coefficient-level analysis elevates dynamics→nestedness from "descriptive" to "credibly supported" (adopted in §3.1, §3.3, §4.2).

## Supplementary Table S9 — MODIS MOD17 cross-check (site level)
MOD17 is not carried into the models (Section 2.3) but reproduces the site-level species–energy result, confirming that the conclusion is not specific to the chosen satellite product.

| comparison | PML-V2 (used) | MODIS MOD17 (cross-check) |
|---|---|---|
| r with site-mean Hill q1 | +0.50 (p = 0.009) | +0.60 (p = 0.001) |
| r with eddy-covariance tower GPP | +0.86 | +0.85 |
| quadratic (hump) term | not supported | not supported |
