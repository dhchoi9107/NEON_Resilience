# Supplementary Information (stub)

Reproducibility: robustness checks = `scripts/obj23/155_methods_robustness.py`; framework = `153_nested_framework.py`.

## Supplementary Methods S1 — Regression-to-the-mean (RTM) null model
Apparent buffering/amplification between initial canopy structure and subsequent structural change was tested against an RTM null. The null assumes no biological feedback between initial state and subsequent change while preserving the observed variance and measurement error of the repeated LiDAR observations; expected initial–change relationships under pure RTM were simulated and compared with the observed slopes. Only observed relationships exceeding the RTM envelope were interpreted as biological regulation. [Full simulation code: `scripts/obj23/90–92`; expand narrative for submission.]

## Supplementary Table S1b — Few-cluster-valid block tests (wild cluster bootstrap, 19 sites)
Restricted wild cluster bootstrap (Rademacher, 999 reps), the primary block test; asymptotic cluster-robust Wald shown for contrast (anticonservative with 19 clusters). Script `159_fewcluster_robust.py`, `wildboot_blocks.csv`.

| response | dynamics: Wald p | **dynamics: wild-boot p** | function: Wald p | **function: wild-boot p** |
|---|---|---|---|---|
| Hill q1 | 2e-4 | **0.058** | 0.059 | **0.128** |
| Hill q2 | 1e-4 | **0.062** | 0.074 | **0.114** |
| LCBD turnover | 0.61 | **0.764** | 4e-4 | **0.003** |
| LCBD nestedness | 0.012 | **0.241** | 0.071 | **0.167** |

Only the productivity→turnover increment survives the wild bootstrap; dynamics increments are marginal (alpha) to non-significant (beta). The variance partition (Table S3) is descriptive and unaffected.

## Supplementary Table S7 — Productivity–turnover: product-specificity & suppression; simple slopes
GPP suppression for LCBD turnover (`gpp_suppression.csv`): r(MODIS, PML) = 0.75; VIF between products = 2.25.

| model | MODIS b | MODIS p | PML b | PML p |
|---|---|---|---|---|
| MODIS only | −0.011 | 0.70 | — | — |
| PML only | — | — | +0.067 | 0.036 |
| joint | −0.051 | — | +0.083 | — |

→ turnover signal is PML-specific; MODIS shows none; joint coefficients oppose (suppression). No directional interpretation.

Simple slopes (`simple_slopes.csv`), Hill q1, ±1 SD moderator:
- GPP × stand age: slope = +0.06 (young, p = 0.64), +0.23 (median, p = 0.08), +0.39 (old, p = 0.004).
- VCI × severity: slope = +0.48 (low, p < 0.001), +0.35 (median, p < 0.001), +0.22 (high, p = 0.003).

## Supplementary Table S1 — Raw vs log-transformed responses (nested LRT p)
Hill numbers log-transformed (satisfies homoscedasticity/normality: Breusch–Pagan p > 0.4, Shapiro p > 0.2); conclusions identical to raw scale.

| response | dynamics (M1→M2) raw | log | function (M2→M3) raw | log |
|---|---|---|---|---|
| Hill q1 | 3.5e-11 | 2.2e-13 | 0.019 | 0.009 |
| Hill q2 | 4.9e-09 | 3.4e-12 | 0.032 | 0.018 |

## Supplementary Table S2 — LCBD: Gaussian OLS vs beta regression (nested LRT p)
LCBD are bounded proportions; beta regression as robustness. Conclusions identical except the nestedness productivity increment (bold), reported conservatively from the Gaussian model in the main text.

| response | model | dynamics (M1→M2) | function (M2→M3) |
|---|---|---|---|
| turnover | OLS | 0.0011 | 8.3e-07 |
| turnover | beta-reg | 2.6e-07 | 3.3e-05 |
| nestedness | OLS | 9.3e-07 | **0.094 (n.s.)** |
| nestedness | beta-reg | 6.2e-12 | **0.028 (sig.)** |

## Supplementary Table S3 — Commonality analysis (shared variance beyond domain)
Main text reports unique (semi-partial R²); shared fractions here. Source: `results/O0_framework/variance_partition.csv`.

| response | unique structure | unique spectral | unique dynamics | unique productivity | shared |
|---|---|---|---|---|---|
| Hill q1 | 0.106 | 0.010 | 0.069 | 0.008 | 0.016 |
| Hill q2 | 0.105 | 0.013 | 0.070 | 0.008 | 0.013 |
| LCBD turnover | 0.020 | 0.005 | 0.017 | 0.029 | 0.011 |
| LCBD nestedness | 0.049 | 0.001 | 0.071 | 0.008 | 0.021 |

## Supplementary Table S6 — Context interactions (framework predictors, cluster-robust, BH per family)
Same retained predictors as the nested models (4 structure + EVI + 2 GPP) × context. 168 terms; 11 survive BH. Script `158_o4_consistent.py`; full table `results/O0_framework/o4_interactions_consistent.csv`.

| family | context | response | predictor | β interaction | q |
|---|---|---|---|---|---|
| stand age | stand_age | Hill q1 | MODIS GPP | +0.093 | 0.008 |
| stand age | stand_age | Hill q1 | PML GPP | +0.063 | 0.013 |
| stand age | stand_age | Hill q2 | MODIS GPP | +0.084 | 0.008 |
| stand age | stand_age | Hill q2 | PML GPP | +0.056 | 0.017 |
| disturbance | severity | Hill q1 | VCI | −0.078 | 0.021 |
| disturbance | severity | Hill q1 | LAI | −0.112 | 0.028 |
| disturbance | severity | Hill q2 | VCI | −0.067 | 0.022 |
| disturbance | severity | Hill q2 | LAI | −0.109 | 0.022 |
| disturbance | recency | Hill q1 | VCI | −0.175 | 0.021 |
| disturbance | recency | Hill q2 | VCI | −0.165 | 0.021 |
| land use | forest fraction | turnover | EVI | −0.064 | <0.01 |

Interpretation: productivity–diversity coupling strengthens with stand age (positive GPP×age); structure–diversity coupling weakens with disturbance severity/recency (negative); EVI–turnover association modified by forest cover. (Supersedes the earlier SAVI-inclusive `interaction_fdr.csv`, which was not predictor-consistent with the nested models.)

## Supplementary Table S5 — Climate control (NEON site MAT/MAP)
Climate = NEON site-level mean annual temperature and precipitation (`data/site_climate_neon.csv`; script `156_climate_baseline.py`). Domain and continuous climate are collinear, so climate is used as an alternative baseline. All climate VIF < 1.3.

**(a) Does GPP add beyond climate alone? (species–energy is not a climate artifact)**

| response | R² climate only | ΔR² GPP beyond climate | p | MODIS β (clustered) | p |
|---|---|---|---|---|---|
| Hill q1 | 0.071 | **0.180** | <0.001 | **+0.197** | 0.019 |
| Hill q2 | 0.055 | 0.159 | <0.001 | +0.165 | 0.039 |
| LCBD turnover | 0.339 | 0.140 | <0.001 | +0.026 | 0.141 |
| LCBD nestedness | 0.022 | 0.027 | 0.003 | −0.014 | 0.061 |

**(b) RS increments under domain vs climate baseline (LRT p)**

| response | dynamics (domain) | dynamics (climate) | function (domain) | function (climate) |
|---|---|---|---|---|
| Hill q1 | <1e-4 | <1e-4 | 0.009 | 1e-4 |
| Hill q2 | <1e-4 | <1e-4 | 0.018 | 4e-4 |
| LCBD turnover | 0.001 | **0.149 (n.s.)** | <1e-4 | <1e-4 |
| LCBD nestedness | <1e-4 | <1e-4 | 0.094 (n.s.) | 0.596 (n.s.) |

Domain baseline R² (e.g., Hill q1 = 0.40) ≫ climate baseline R² (0.07): domain fixed effects absorb climatic **and** non-climatic biogeography (species pools, soils), so domain is retained as the primary, stronger control and climate serves as the explicit "beyond-climate" check. Note: the LCBD-turnover dynamics increment is baseline-dependent (significant under domain, n.s. under climate); its productivity increment is robust to both.

## Supplementary Table S4 — Mixed vs clustered-OLS fixed effects
Site random-intercept variance estimated at 0 (singular); standardized fixed-effect coefficients closely similar (mean |Δ| = 0.03, max |Δ| = 0.05, r = 0.93 across 14 RS predictors).
