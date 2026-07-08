"""
NEON_v2 STEP 2 — Build final pooled dataset.
Combines:
  - Diversity responses + LiDAR predictors (mean/sd/trend) : reused, BRDF-irrelevant
    (lidar_pooled_predictors.csv, derived from prior pooled diversity work)
  - VI predictors (mean/sd/trend) : NEWLY computed from NEON .002 BRDF VI ONLY
    (plot_vi_neon_brdf.csv, single source = DP3.30026.002 bidirectional, 8 BRDF years)

Spectral temporal features are over the 8 NEON-BRDF years only -> 100% consistent source.
Output: NEON_v2/data/FINAL_v2_pooled.csv
"""
import os, numpy as np, pandas as pd
D = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"

vi = pd.read_csv(os.path.join(D, "plot_vi_neon_brdf.csv"))
VIS = ["NDVI", "EVI", "ARVI", "SAVI"]
print(f"VI per-plot-year rows: {len(vi)} | years: {sorted(vi['year'].unique())} | plots: {vi['plotID'].nunique()}")

def slope(y, x):
    if len(x) < 3: return np.nan
    return np.polyfit(x, y, 1)[0]

recs = []
for pid, g in vi.groupby("plotID"):
    g = g.sort_values("year")
    rec = {"plotID": pid, "VI_nyears": g["year"].nunique(),
           "VI_years": ",".join(map(str, g["year"].tolist()))}
    for v in VIS:
        col = f"{v}_mean"
        vals = g[col].values; yrs = g["year"].values
        ok = np.isfinite(vals)
        rec[f"{v}_mean"] = float(np.nanmean(vals)) if ok.sum() else np.nan
        rec[f"{v}_sd"] = float(np.nanstd(vals)) if ok.sum() >= 2 else np.nan
        rec[f"{v}_trend"] = slope(vals[ok], yrs[ok]) if ok.sum() >= 3 else np.nan
    recs.append(rec)
vif = pd.DataFrame(recs)
print(f"VI plot-level features: {len(vif)} plots | nyears dist:\n{vif['VI_nyears'].value_counts().sort_index().to_string()}")

lidar = pd.read_csv(os.path.join(D, "lidar_pooled_predictors.csv"))
final = lidar.merge(vif, on="plotID", how="left")
final["spectral_source"] = "NEON_DP3.30026.002_bidirectional_BRDF_only"
out = os.path.join(D, "FINAL_v2_pooled.csv")
final.to_csv(out, index=False)
print(f"\nDONE -> {out}")
print(f"shape: {final.shape} | plots with VI: {final['NDVI_mean'].notna().sum()} | cov>=0.9: {(final['sample_coverage']>=0.9).sum()}")
