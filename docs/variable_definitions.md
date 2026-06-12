# Variable Definitions & Computation Methods

## 1. Structural Diversity Metrics (compute_fsd.py)

Source: NEON LiDAR point clouds (DP1.30003.001), 10m grid, Sentinel-2 aligned.

### Pre-processing
- Height normalization via ground-point (Class 2) TIN interpolation
- Classification filter: Classes 1-5 only (unclassified, ground, low/med/high veg)
- Height range: -1m to 100m; vegetation cutoff Z0 = 0.5m
- Per-cell 6-sigma outlier removal

### CHM Metrics
| Band | Name | Formula | Description |
|------|------|---------|-------------|
| 1 | rumple | `1 + sd(z_all) / mean(z_all)` | Surface roughness index (all z >= 0) |
| 2 | top_rugosity | `sd(z_veg)` | SD of vegetation heights (z >= 0.5m) |
| 3 | mean_max_canopy_ht | `mean(z_veg)` | Mean vegetation height |
| 4 | max_canopy_ht | `max(z)` | Maximum height per cell |
| 5 | deepgap_fraction | `count(z < 0.5) / count(all)` | Proportion of sub-canopy points |

### Structural Metrics
| Band | Name | Formula | Description |
|------|------|---------|-------------|
| 6 | meanH | `mean(z_veg)` | Mean vegetation height (same as band 3) |
| 7 | vert_sd | `sd(z_veg)` | Vertical standard deviation |
| 8 | vertCV | `sd(z_veg) / mean(z_veg)` | Coefficient of variation |
| 9 | mean_sd | `mean(sd_subcells)` | Mean of 3m sub-cell SDs within 10m cell |
| 10 | sd_sd | `sd(sd_subcells)` | SD of 3m sub-cell SDs |
| 11 | GC | Gini coefficient | `(2*sum(z_sorted * rank) / sum(z_sorted) - (n+1)) / (n-1)` |

### LAI/Structure (Beer-Lambert LAD model, k=1.0)
| Band | Name | Formula | Description |
|------|------|---------|-------------|
| 12 | GFP | Gap Fraction Profile | Mean of `(1 - n_bin / n_above)` per 1m bin above Z0 |
| 13 | VCI | Vertical Complexity Index | Shannon entropy of point distribution in 1m bins, normalized |
| 14 | q25 | 25th percentile | Height quantile |
| 15 | q50 | 50th percentile (median) | Height quantile |
| 16 | q75 | 75th percentile | Height quantile |
| 17 | q95 | 95th percentile | Height quantile |
| 18 | HeightRatio | `count(z < 10m) / total * 100` | Percentage of points below 10m |
| 19 | FHD | Foliage Height Diversity | Shannon entropy of normalized LAD: `-sum(p * log(p))` where `p = LAD_i / sum(LAD)` |
| 20 | LAI | Leaf Area Index | `sum(LAD)` from Z0 upward; `LAD = log(pulse_in / pulse_out) / k` |
| 21 | LAI_subcanopy | Subcanopy LAI | `sum(LAD)` from Z0 to 5.0m |

---

## 2. Taxonomic Diversity (compute_taxonomic_diversity.py)

Source: NEON vegetation structure (DP1.10098.001) - tree inventory in 40x40m plots.

### Data Preparation
- Join `vst_mappingandtagging` (species ID) with `vst_apparentindividual` (annual measurements)
- Filter: live stems only (`plantStatus` contains "Live")
- Remove unknown/unidentified taxa
- Group by `plotID x year x taxonID`

### Alpha Diversity (per plot, per year)
| Metric | Formula | Reference |
|--------|---------|-----------|
| Richness (S) | Count of unique taxonID | - |
| Shannon (H') | `-sum(p_i * ln(p_i))` where `p_i = n_i / N` | Shannon 1948 |
| Simpson (1-D) | `1 - sum(p_i^2)` | Simpson 1949 |
| Abundance (N) | Total stem count | - |

### Beta Diversity (pairwise between plots within site)
| Metric | Formula | Reference |
|--------|---------|-----------|
| Bray-Curtis | `sum(|x_i - y_i|) / sum(x_i + y_i)` | Bray & Curtis 1957 |
| Jaccard | `1 - |A intersect B| / |A union B|` (presence/absence) | Jaccard 1912 |
| Sorensen (beta_sor) | `(b + c) / (2a + b + c)` | Sorensen 1948 |
| Turnover (beta_sim) | `min(b, c) / (a + min(b, c))` | Baselga 2010 |
| Nestedness (beta_sne) | `beta_sor - beta_sim` | Baselga 2010 |

Where: `a` = shared species, `b` = species only in site 1, `c` = species only in site 2.

### Gamma Diversity (per site, per year)
- Total species richness across all plots at a site
- Additive partitioning: gamma = mean(alpha) + beta

---

## 3. Spectral Diversity (compute_spectral_diversity.py)

Source: NEON AOP vegetation index products (1m resolution, 1km tiles).

### Input Bands
| Product | Bands |
|---------|-------|
| DP3.30026.002 (VI) | NDVI, EVI, ARVI, PRI, SAVI |
| DP3.30012.002 | LAI |
| DP3.30014.002 | fPAR |

### 10m Aggregation
Each 1m tile (1000x1000) is aggregated to 10m (100x100):
- **Mean** and **SD** per 10m cell per band
- Output: 14-band GeoTIFF (7 bands x 2 stats)

### Site-level Spectral Diversity Metrics
| Metric | Formula | Description |
|--------|---------|-------------|
| Rao's Q | `mean(d_ij)` for sampled pixel pairs | Mean pairwise Euclidean distance in 7D spectral space |
| Spectral CV | `mean(sd_band / mean_band)` across bands | Mean coefficient of variation |
| Spectral Shannon | Mean `-sum(p * log(p))` of binned values per band | Spectral entropy |

### PCA-based Functional Spectral Diversity
1. Stack 7 bands per pixel → PCA (site-level)
2. Extract PC1-PC3
3. Per plot (40x40m = ~1600 pixels at 1m):

| Metric | Formula | Reference |
|--------|---------|-----------|
| Spectral FRic | Convex hull volume in PC1-PC3 space | Villeger et al. 2008 |
| Spectral FDiv | `mean(dist_to_centroid) / max(dist_to_centroid)` | Villeger et al. 2008 |
| Spectral FEve | MST edge length regularity: `1 - sum(|p_i - 1/n|) / (2*(1-1/n))` | Villeger et al. 2008 |

---

## 4. Functional Diversity (compute_functional_diversity.py)

Source: Combined structural (FSD) + spectral (VI PCA) traits per plot.

### Trait Space
| Category | Traits | Source |
|----------|--------|--------|
| Structural | mean_max_canopy_ht, FHD, LAI, vert_sd, GC | FSD 10m raster |
| Spectral | NDVI, EVI, PRI, LAI, fPAR (mean per 10m cell) | Spectral 10m raster |

- All traits z-score normalized before analysis
- Per-plot: 4x4 = 16 pixels (40m plot / 10m cell)

### Metrics
| Metric | Formula | Reference |
|--------|---------|-----------|
| FRic | Convex hull volume in trait space (first 3 dims) | Villeger et al. 2008 |
| FDiv | `mean(dist_to_centroid) / max(dist_to_centroid)` | Villeger et al. 2008 |
| FEve | MST edge regularity (see spectral FEve) | Villeger et al. 2008 |
| Rao's Q | Mean pairwise Euclidean distance | Botta-Dukat 2005 |
| Beta functional | Mean pairwise centroid distance between plots | Laliberte et al. 2020 |

---

## 5. Environmental Heterogeneity (compute_env_heterogeneity.py)

### 5a. Compositional Heterogeneity
Source: FSD rasters (mean_max_canopy_ht band).

| Metric | Formula | Description |
|--------|---------|-------------|
| chm_mean | `mean(CHM)` | Site-level mean canopy height |
| chm_sd | `sd(CHM)` | Site-level canopy height variability |
| chm_cv | `sd / mean` | Coefficient of variation |
| chm_range | `max - min` | Total height range |
| chm_iqr | `Q75 - Q25` | Interquartile range |

### 5b. Configurational Heterogeneity (Gradient Surface Metrics)
Source: FSD canopy height with moving window analysis.
Reference: Smith et al. 2021 (geodiv R package).

| Metric | Formula | Window | Description |
|--------|---------|--------|-------------|
| Sa | `mean(|z - z_local_mean|)` | 100m, 500m, 1km | Average roughness |
| Sq | `sqrt(mean((z - z_local_mean)^2))` | 100m, 500m, 1km | RMS roughness |
| Ssk | `mean(dev^3) / Sq^3` | 100m, 500m, 1km | Surface skewness |
| Sku | `mean(dev^4) / Sq^4` | 100m, 500m, 1km | Surface kurtosis |

### 5c. Fragmentation
Source: NLCD 2021 (30m land cover), downloaded via MRLC WCS.

| Metric | Formula | Description |
|--------|---------|-------------|
| nlcd_shdi | `-sum(p_i * ln(p_i))` across land cover classes | Shannon Landscape Diversity Index |
| nlcd_patch_density | `n_patches / area_km2` | Patch density |
| nlcd_n_classes | Count of distinct NLCD classes | Class richness |
| nlcd_forest_proportion | `forest_pixels / total_pixels` | Forest coverage (classes 41-43, 90) |

---

## 6. Productivity / Dynamic Habitat Indices (compute_productivity_dhi.py)

Source: Multi-temporal FSD rasters (mean_max_canopy_ht band across years).
Reference: Radeloff et al. 2019; Coops et al. 2008.

### Height Change Computation
- For consecutive FSD years: `delta_h = CHM(t2) - CHM(t1) / (year2 - year1)` (annual rate)
- Spatial alignment: intersection bounding box across all years per site

### DHI Metrics (per pixel, then aggregated)
| Metric | Formula | Interpretation |
|--------|---------|---------------|
| Cumulative DHI | `sum(max(delta_h, 0))` across years | Total growth (available energy) |
| Minimum DHI | `min(delta_h)` across years | Environmental stress indicator |
| Variation DHI | `sd(delta_h)` across years | Environmental stability |
| Height trend | Linear slope of `height ~ year` | Long-term growth rate (m/yr) |

---

## 7. Mixed-Effects Models (run_mixed_models.R)

### Model 1a: Alpha Diversity ~ RS Diversity
```
alpha_shannon ~ spectral_RaoQ * structural_CHM_CV + Sa_500m + (1|domain/site)
```

### Model 1b: Beta Diversity ~ RS Diversity
```
beta_bray ~ spectral_RaoQ * structural_CHM_CV + Sa_500m + (1|domain/site)
```

### Model 2: Productivity ~ Functional Diversity
```
cumulative_DHI ~ functional_FRic * CHM_CV + (1|domain)
```

All predictors z-score standardized. Random effects: site nested within NEON domain.
Diagnostics: residual plots, Q-Q, VIF check.
R-squared: marginal (fixed effects) and conditional (fixed + random) via MuMIn::r.squaredGLMM.
