# NEON Resilience

NEON AOP LiDAR & spectral data processing pipeline for forest structural diversity analysis.

## Scripts

| Script | Description |
|--------|-------------|
| `compute_fsd.py` | Compute 21-band structural diversity metrics from LiDAR point clouds (10m, Sentinel-2 aligned) |
| `neon_lidar_batch.py` | Batch download NEON discrete-return LiDAR point clouds (DP1.30003.001) |
| `neon_veg_indices_batch.py` | Batch download NEON vegetation indices & related spectral products (BRDF-corrected) |

## Structural Diversity Metrics (21 bands)

| Group | Metrics |
|-------|---------|
| CHM | rumple, top_rugosity, mean_max_canopy_ht, max_canopy_ht, deepgap_fraction |
| Structural | meanH, vert_sd, vertCV, mean_sd, sd_sd, GC (Gini coefficient) |
| LAI/Structure | GFP, VCI, FHD, LAI, LAI_subcanopy |
| Quantiles | q25, q50, q75, q95 |
| Additional | HeightRatio |

## Downloaded Vegetation Products (BRDF-corrected, .002)

| Product | NEON ID | Description |
|---------|---------|-------------|
| VI | DP3.30026.002 | NDVI, EVI, ARVI, PRI, NDLI, NDNI, SAVI |
| LAI | DP3.30012.002 | Leaf Area Index |
| fPAR | DP3.30014.002 | Fraction of Photosynthetically Active Radiation |
| CWI | DP3.30019.002 | Canopy Water Indices |

## NEON API Token

Token stored at: `E:/neon_lidar/_code/.neon_token`  
Obtain from: https://data.neonscience.org → My Account → API Token

## Data Locations

- LiDAR point clouds: `E:/neon_lidar/DP1.30003.001/`
- Structural diversity rasters: `E:/neon_lidar/structural_diversity/`
- Vegetation indices: `E:/neon_lidar/vegetation_indices/`
