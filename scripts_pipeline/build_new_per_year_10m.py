"""
Build per_year rows for the 7 NEW sites from 10m FSD tiles, using the EXACT
original method (_ARCHIVE_v1/NEON_FINAL/scripts/build_per_year.py):
  - 40x40m window (HALF=20) centered on plot easting/northing
  - np.nanmean per SKEEP band; plot center must fall within the raster bounds
  - VI (NDVI/EVI/ARVI/SAVI) merged same-plot-year from plot_vi_neon_brdf.csv
  - domain from vst_perplotperyear domainID
Output (NEW file, never overwrites per_year_v2):
  scripts_pipeline/_pipeline_state/per_year_newsites_10m.csv
Columns match per_year_v2.csv exactly.
"""
import os, sys, glob, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
sys.path.insert(0, r"C:\Users\star1\Documents\GitHub\NEON_Resilience")
from site_config import VEG_STRUCT_DIR

NEW = ["DELA", "LENO", "UKFS", "YELL", "BONA", "DEJU", "HEAL"]
HALF = 20  # 40x40m window, identical to original build_per_year
FSD_DIR = "E:/neon_lidar/structural_diversity"

# 21-band FSD raster order (band index = list.index(name)+1), identical to original
FSD_BANDS = ['rumple','top_rugosity','mean_max_canopy_ht','max_canopy_ht','deepgap_fraction','meanH',
 'vert_sd','vertCV','mean_sd','sd_sd','GC','GFP','VCI','q25','q50','q75','q95','HeightRatio','FHD','LAI','LAI_subcanopy']
SKEEP = {'Canopy_Ht':'mean_max_canopy_ht','Max_Ht':'max_canopy_ht','Rumple':'rumple','Rugosity':'top_rugosity',
 'Deep_Gap':'deepgap_fraction','Vert_SD':'vert_sd','Vert_CV':'vertCV','Gini':'GC','VCI':'VCI','FHD':'FHD',
 'LAI':'LAI','Q95':'q95','Ht_Ratio':'HeightRatio'}
# final column order == per_year_v2.csv
COLS = ['plotID','siteID','year','Canopy_Ht','Max_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Vert_CV',
        'Gini','VCI','FHD','LAI','Q95','Ht_Ratio','NDVI','EVI','ARVI','SAVI','domain']

# --- plot coords + domain (one row per plot) ---
ppy = pd.read_csv(str(VEG_STRUCT_DIR) + "/vst_perplotperyear.csv", low_memory=False)
meta = (ppy[ppy.siteID.isin(NEW)][['siteID','plotID','easting','northing','domainID']]
        .dropna(subset=['easting','northing']).drop_duplicates('plotID'))
coords = meta.set_index('plotID')
print(f"new-site plots with coords: {len(coords)}", flush=True)

# --- index 10m FSD rasters by (site, year); tile name = {year}_{site}_{N}_FSD_10m.tif ---
fsd_idx = {}
for f in glob.glob(os.path.join(FSD_DIR, "*_FSD_10m.tif")):
    p = os.path.basename(f).split("_")
    if p[1] in NEW:
        fsd_idx.setdefault((p[1], int(p[0])), []).append(f)
years_by_site = {}
for (s, y) in fsd_idx:
    years_by_site.setdefault(s, set()).add(y)
for s in NEW:
    print(f"  {s}: FSD years {sorted(years_by_site.get(s, []))}", flush=True)

def fsd_year(site, year, e, n):
    out = {k: np.nan for k in SKEEP}
    for f in fsd_idx.get((site, year), []):
        try:
            with rasterio.open(f) as ds:
                if not (ds.bounds.left <= e <= ds.bounds.right and ds.bounds.bottom <= n <= ds.bounds.top):
                    continue
                w = from_bounds(e - HALF, n - HALF, e + HALF, n + HALF, ds.transform)
                for lab, bn in SKEEP.items():
                    arr = ds.read(FSD_BANDS.index(bn) + 1, window=w).astype(float)
                    if np.isfinite(arr).any():
                        out[lab] = float(np.nanmean(arr))
                return out
        except Exception:
            continue
    return out

# --- one row per (plot, FSD year) where the plot is covered ---
rows = []
for pid, r in coords.iterrows():
    site, e, n, dom = r.siteID, float(r.easting), float(r.northing), r.domainID
    for year in sorted(years_by_site.get(site, [])):
        vals = fsd_year(site, year, e, n)
        if all(pd.isna(v) for v in vals.values()):
            continue  # plot not covered by this year's tiles
        rec = {'plotID': pid, 'siteID': site, 'year': int(year), 'domain': dom}
        rec.update(vals)
        rows.append(rec)

fsd = pd.DataFrame(rows)
print(f"\nFSD plot-years aggregated: {len(fsd)}", flush=True)

# --- merge VI (same plot-year) ---
vi = pd.read_csv("NEON_v2/data/plot_vi_neon_brdf.csv")
vic = vi[["plotID","year","NDVI_mean","EVI_mean","ARVI_mean","SAVI_mean"]].rename(
    columns={"NDVI_mean":"NDVI","EVI_mean":"EVI","ARVI_mean":"ARVI","SAVI_mean":"SAVI"})
out = fsd.merge(vic, on=["plotID","year"], how="left")
out = out.reindex(columns=COLS)

OUT = "scripts_pipeline/_pipeline_state/per_year_newsites_10m.csv"
out.to_csv(OUT, index=False)
print(f"\n-> {OUT}", flush=True)
print(f"rows: {len(out)}, sites: {out.siteID.nunique()}, plots: {out.plotID.nunique()}", flush=True)
print(f"structure non-null: Canopy_Ht {out.Canopy_Ht.notna().sum()}, VCI {out.VCI.notna().sum()}, LAI {out.LAI.notna().sum()}", flush=True)
print(f"VI matched: NDVI {out.NDVI.notna().sum()}/{len(out)}", flush=True)
print(out.groupby('siteID').agg(plot_years=('plotID','size'), plots=('plotID','nunique'),
      yrs=('year','nunique'), medHt=('Canopy_Ht','median')).round(2).to_string(), flush=True)
