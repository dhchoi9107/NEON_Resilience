"""
Rebuild FINAL_v2_pooled_26.csv = lidar_pooled_predictors_26.csv + pooled VI features
(from plot_vi_neon_brdf.csv). Mirrors NEON_v2/scripts/02_build_dataset.py exactly,
but for the 26-site file. Needed after TREE spectral was recovered from the STEI
flight box (TREE VI merged into plot_vi_neon_brdf.csv).

Matches the EXISTING FINAL_v2_pooled_26 schema (reindex to its columns), so all
downstream _26 scripts (153-159) keep working unchanged. Backs up the old file first.
"""
import os, shutil, numpy as np, pandas as pd
D = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
OUT = os.path.join(D, "FINAL_v2_pooled_26.csv")

vi = pd.read_csv(os.path.join(D, "plot_vi_neon_brdf.csv"))
VIS = ["NDVI", "EVI", "ARVI", "SAVI"]
print(f"VI rows: {len(vi)} | sites: {vi['siteID'].nunique()} | plots: {vi['plotID'].nunique()}")

def slope(y, x):
    if len(x) < 3: return np.nan
    return np.polyfit(x, y, 1)[0]

recs = []
for pid, g in vi.groupby("plotID"):
    g = g.sort_values("year")
    rec = {"plotID": pid, "VI_nyears": g["year"].nunique()}
    for v in VIS:
        vals = g[f"{v}_mean"].values; yrs = g["year"].values
        ok = np.isfinite(vals)
        rec[f"{v}_mean"] = float(np.nanmean(vals)) if ok.sum() else np.nan
        rec[f"{v}_sd"] = float(np.nanstd(vals)) if ok.sum() >= 2 else np.nan
        rec[f"{v}_trend"] = slope(vals[ok], yrs[ok]) if ok.sum() >= 3 else np.nan
    recs.append(rec)
vif = pd.DataFrame(recs)

lidar = pd.read_csv(os.path.join(D, "lidar_pooled_predictors_26.csv"))
final = lidar.merge(vif, on="plotID", how="left")
final["spectral_source"] = "NEON_DP3.30026.002_bidirectional_BRDF_only"

# Preserve exact existing schema (column set + order) so downstream _26 scripts are unaffected.
old = pd.read_csv(OUT)
missing = [c for c in old.columns if c not in final.columns]
extra = [c for c in final.columns if c not in old.columns]
assert not missing, f"rebuild is missing columns present before: {missing}"
if extra:
    print(f"note: dropping columns not in prior schema: {extra}")
final = final.reindex(columns=old.columns)

shutil.copy2(OUT, OUT + ".bak")
print(f"backup -> {OUT}.bak")
final.to_csv(OUT, index=False)

print(f"\nDONE -> {OUT}")
print(f"shape: {final.shape} | sites: {final.siteID.nunique()}")
print(f"plots with spectral: {final['NDVI_mean'].notna().sum()} (was {old['NDVI_mean'].notna().sum()})")
for s in ["TREE", "ORNL"]:
    sub = final[final.siteID == s]
    print(f"  {s}: {sub['NDVI_mean'].notna().sum()}/{len(sub)} plots with spectral")
sites_missing = sorted(final.groupby("siteID")["NDVI_mean"].apply(lambda x: x.notna().sum() == 0).loc[lambda z: z].index.tolist())
print(f"sites still fully missing spectral: {sites_missing}")
