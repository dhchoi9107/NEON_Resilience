# 2. Materials and Methods — 26-site update (drop-in sections)

> Only the sections that change with the 26-site expansion are given here, ready to
> replace the corresponding §2 blocks in MANUSCRIPT.md. Key substantive changes:
> (1) 19→26 sites and the study now spans **interior-Alaska boreal** forest (D19), so
> "conterminous United States" is corrected; (2) TREE's LiDAR and spectral predictors
> are sourced from the shared **STEI AOP flight box**; (3) **ORNL** lacks the .002
> spectral product and is excluded from spectral models; (4) cluster counts 19→26 and
> few-cluster language updated. Tower-GPP site count tagged **[re-run 26]** pending recount.

---

### 2.1 Study sites

We analyzed 26 forested sites of the U.S. National Ecological Observatory Network (NEON), spanning temperate and boreal biomes across 12 NEON eco-climatic domains from the conterminous United States to interior Alaska. Sites range from evergreen-needleleaf forests of the Pacific Northwest and Rocky Mountains, through deciduous-broadleaf and mixed forests of the eastern and southeastern United States, to boreal forests of interior Alaska (NEON Domain 19: BONA, DEJU, HEAL), encompassing broad gradients in climate, productivity, stand age, and species composition. Analyses were conducted at the level of NEON base and distributed plots and aggregated to the site and domain levels as appropriate (below).

Two NEON sites share an Airborne Observation Platform (AOP) flight box with a neighboring site, so their airborne data are published under the neighbor's site code: TREE (Treehaven) is flown within the STEI (Steigerwaldt) flight box. For TREE we therefore derived all airborne predictors (LiDAR structure and dynamics, spectral indices) from the STEI acquisitions, sampling them at TREE plot coordinates. One site (ORNL) lacked the bidirectional-reflectance (.002) spectral product entirely and was therefore included in structural and productivity analyses but excluded from analyses involving spectral predictors (spectral models: 25 sites).

### 2.3 Remote-sensing predictors

**Static canopy structure.** Airborne LiDAR from the NEON AOP was summarized into canopy structural metrics spanning height (mean and 95th-percentile canopy height, maximum height), vertical organization (foliage height diversity, vertical canopy index, vertical SD and CV), heterogeneity (Gini coefficient, rumple, rugosity, height ratio), gap structure (deep-gap fraction), and leaf area index, computed from 10 m foliage-structural-diversity rasters and sampled within a 40 m window at each plot. For the static-state analysis, each metric was represented by its multi-year mean at each plot. For TREE, structural metrics were computed from the STEI flight-box acquisitions (seven years, 2016–2025) at TREE plot coordinates.

**Static spectral state.** Spectral vegetation indices were computed from a single, consistent source — NEON bidirectional-reflectance (BRDF-corrected) surface reflectance (the .002 product) — as individual indices (NDVI, EVI, ARVI, SAVI). We deliberately used individual vegetation indices rather than a combined spectral-diversity metric, and excluded PRI and fPAR because of cross-year inconsistency in the provisional products. For TREE, spectral indices were extracted from the STEI flight-box VI tiles (four years) at TREE plot coordinates. ORNL was excluded from spectral analyses because the .002 VI product is not published for that site.

**Structural dynamics.** Repeat LiDAR acquisitions were used to estimate the interannual trend (per-year slope) of each structural metric at each plot, over the available observation window. These trends index the *rate and direction* of canopy structural change (e.g., height and leaf-area accumulation, gap opening) independent of the mean state. For TREE, trends were estimated from the STEI flight-box time series.

**Ecosystem productivity.** We used three independent estimates of gross primary productivity (GPP), none derived from vegetation greenness indices. Two are satellite-constrained model products extracted at plots via Google Earth Engine: MODIS MOD17 GPP and the Penman–Monteith–Leuning PML-V2 GPP (CAS/IGSNRR/PML/V2_v017), each averaged over available years. The third is directly measured eddy-covariance GPP from NEON flux towers (AmeriFlux FLUXNET-1F, GPP_NT_VUT_REF), available as a site-level multi-year mean for 22 of the 26 sites (four sites lacked a processed FLUXNET-1F GPP record: LENO, ORNL, SOAP, TEAK); no imputation was performed for sites lacking tower measurements. The boreal Alaska sites (BONA, DEJU, HEAL) extend the tower-GPP gradient down to ~380–540 gC m⁻² yr⁻¹, roughly a quarter of the most productive temperate sites (~2350).

### 2.4 Contextual covariates — Stand age (revised clause)

… Effective stand age was obtained from the GAMI global forest-age product (Global Age Mapping Integration, v3.1; 100 m; ensemble mean of 20 members for 2020), sampled at all 26 sites' plot coordinates from the EO-Forest STAC Zarr archive (GFZ Potsdam); re-extracted values reproduced the earlier 19-site extraction well (r = 0.83). A Global 30 m forest-age layer was retained as a complementary recency proxy but saturates at ~40 yr and was not used for the succession-trend analysis, which requires true continuous age. Because stand age carried substantial within-site homogeneity (ICC ≈ 0.75), age effects were interpreted at the site level (effective n ≈ 26 sites), and all age analyses were re-run with 20-year age bins to confirm robustness to age-estimate uncertainty.

### 2.5 Statistical analysis (revised clauses)

*(Predictor standardization, VIF screening, domain fixed effects + site-clustered SE, singular-mixed-model check, residual diagnostics, Moran's I, and FDR correction are unchanged except for the cluster count and few-cluster language below.)*

Relationships were estimated with linear models that included NEON domain as a fixed effect, separating broad biogeographic differences from within-domain ecological relationships, with standard errors clustered by site (26 clusters) to account for the non-independence of plots within sites. …

**Nested framework and few-cluster inference.** Analyses followed a nested hierarchical framework (M0 domain; M1 +static structure +spectral; M2 +structural dynamics; M3 +GPP) fitted over a common sample (n = 644 alpha plots, 579 beta plots). Because inference rested on 26 site clusters — a regime in which asymptotic cluster-robust standard errors remain anticonservative — block-level significance was assessed primarily with restricted wild cluster bootstrap tests (Rademacher weights, 999 replicates) of the null that the added block's coefficients were jointly zero; asymptotic cluster-robust joint Wald (F) tests, likelihood-ratio tests, and AIC are reported as secondary summaries. Both sequential partial R² and blockwise semi-partial (unique) R² are reported and labelled explicitly.

*(§2.5 O1–O4 sub-paragraphs and the interaction-screen count are unchanged; the O4 screen comprised 168 interaction terms as before, of which twelve survived FDR in the 26-site analysis — reconcile the surviving-term count in Results §3.5.)*

---

## Summary of Methods edits vs. the 19-site manuscript

| Location | Change |
|---|---|
| §2.1 | "19 … conterminous United States" → "26 … conterminous US to interior Alaska (12 domains)"; added TREE flight-box + ORNL spectral-exclusion paragraph |
| §2.3 structure/spectral/dynamics | Added TREE-via-STEI provenance; ORNL spectral exclusion (spectral models = 25 sites) |
| §2.3 productivity | Tower-GPP "16 of 19 sites" → **22 of 26 sites** (AmeriFlux zips already on disk in E:\FLUX; +6 new towers BONA/DEJU/HEAL/DELA/UKFS/YELL; LENO/ORNL/SOAP/TEAK lack GPP). Saved `data/site_tower_gpp_26.csv` |
| §2.4 | "effective n ≈ 19" → "≈ 26 sites" |
| §2.5 | "19 site clusters" → "26"; n = 644 alpha / 579 beta; few-cluster wording retained |
| §2.5 O4 | surviving interaction count 11 → 12 (reconcile with Results) |

**Resolved for all 26 sites:** tower GPP (22/26, `site_tower_gpp_26.csv`); stand age (GAMI v3.1 STAC Zarr, `plot_stand_age_gami_26.csv`, r = 0.83 vs original); land cover (ESA WorldCover, `plot_landuse_het_26.csv`); disturbance (NBR + Hansen + MTBS → `plot_disturbance_robust_26.csv`, 225 disturbed) and BACI (`plot_baci_diversity_26.csv`). O4/context analyses re-run (140/158/159_26, 71_26).

**Only remaining:** reconcile the land-use *direct*-effect model spec (26-site magnitude smaller and n.s.). Add to §2.3/§2.4: disturbance and land-cover products re-extracted for all 26 sites; single-digit UTM zones (Alaska, e.g. 6N) handled in the coordinate parser.
