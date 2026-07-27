"""
Recover TREE structural predictors from the STEI flight-box 10m FSD mosaics.
TREE (Treehaven) is flown as part of STEI's AOP box, so its LiDAR is published under
STEI. We sample STEI's 10m FSD mosaics at TREE plot centers (40x40m window), exactly
mirroring build_new_per_year_10m.py, then pool (mean/sd/trend/nyears) exactly like
build_new_pooled_26.py, and fill TREE's currently-NaN structure columns in
lidar_pooled_predictors_26.csv. Also appends TREE rows to per_year_v2_26.csv.

Never blind-overwrites: backs up each edited file to <file>.bak2.
"""
import os, glob, sys, shutil, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
sys.path.insert(0, r"C:\Users\star1\Documents\GitHub\NEON_Resilience")
from site_config import VEG_STRUCT_DIR

SITE = "TREE"
BOX = "STEI"          # flight-box site code that owns TREE's airborne tiles
HALF = 20             # 40x40m window, identical to build_new_per_year_10m
FSD_DIR = "E:/neon_lidar/structural_diversity"
D = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"

FSD_BANDS = ['rumple','top_rugosity','mean_max_canopy_ht','max_canopy_ht','deepgap_fraction','meanH',
 'vert_sd','vertCV','mean_sd','sd_sd','GC','GFP','VCI','q25','q50','q75','q95','HeightRatio','FHD','LAI','LAI_subcanopy']
SKEEP = {'Canopy_Ht':'mean_max_canopy_ht','Max_Ht':'max_canopy_ht','Rumple':'rumple','Rugosity':'top_rugosity',
 'Deep_Gap':'deepgap_fraction','Vert_SD':'vert_sd','Vert_CV':'vertCV','Gini':'GC','VCI':'VCI','FHD':'FHD',
 'LAI':'LAI','Q95':'q95','Ht_Ratio':'HeightRatio'}
STRUCT = list(SKEEP.keys())

# per_year_v2 column order
PY_COLS = ['plotID','siteID','year','Canopy_Ht','Max_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Vert_CV',
           'Gini','VCI','FHD','LAI','Q95','Ht_Ratio','NDVI','EVI','ARVI','SAVI','domain']

# --- TREE plot coords + domain ---
ppy = pd.read_csv(str(VEG_STRUCT_DIR) + "/vst_perplotperyear.csv", low_memory=False)
meta = (ppy[ppy.siteID == SITE][['siteID','plotID','easting','northing','domainID']]
        .dropna(subset=['easting','northing']).drop_duplicates('plotID').set_index('plotID'))
print(f"{SITE} plots with coords: {len(meta)}", flush=True)

# --- index STEI 10m mosaics by year: {year}_{BOX}_{N}_FSD_10m.tif ---
fsd_idx = {}
for f in glob.glob(os.path.join(FSD_DIR, f"*_{BOX}_*_FSD_10m.tif")):
    p = os.path.basename(f).split("_")
    fsd_idx.setdefault(int(p[0]), []).append(f)
years = sorted(fsd_idx)
print(f"{BOX} 10m mosaic years: {years}", flush=True)

def fsd_year(year, e, n):
    out = {k: np.nan for k in SKEEP}
    for f in fsd_idx.get(year, []):
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

# --- per-year rows for TREE ---
rows = []
for pid, r in meta.iterrows():
    e, n, dom = float(r.easting), float(r.northing), r.domainID
    for year in years:
        vals = fsd_year(year, e, n)
        if all(pd.isna(v) for v in vals.values()):
            continue
        rec = {'plotID': pid, 'siteID': SITE, 'year': int(year), 'domain': dom}
        rec.update(vals)
        rows.append(rec)
py = pd.DataFrame(rows)
print(f"TREE per-year structure rows: {len(py)} | plots covered: {py.plotID.nunique()} | years: {sorted(py.year.unique())}", flush=True)

# --- merge VI (same plot-year) into per_year rows ---
vi = pd.read_csv(os.path.join(D, "plot_vi_neon_brdf.csv"))
vic = vi[["plotID","year","NDVI_mean","EVI_mean","ARVI_mean","SAVI_mean"]].rename(
    columns={"NDVI_mean":"NDVI","EVI_mean":"EVI","ARVI_mean":"ARVI","SAVI_mean":"SAVI"})
py_out = py.merge(vic, on=["plotID","year"], how="left").reindex(columns=PY_COLS)

# --- pool to mean/sd/trend/nyears (identical to build_new_pooled_26) ---
def slope(y, x):
    ok = np.isfinite(y)
    return float(np.polyfit(x[ok], y[ok], 1)[0]) if ok.sum() >= 3 else np.nan
prows = []
for pid, g in py.groupby('plotID'):
    g = g.sort_values('year'); yrs = g['year'].values.astype(float)
    rec = {'plotID': pid}
    for m in STRUCT:
        v = g[m].values.astype(float); ok = np.isfinite(v)
        rec[f"{m}_mean"] = float(np.nanmean(v)) if ok.sum() else np.nan
        rec[f"{m}_sd"] = float(np.nanstd(v)) if ok.sum() >= 2 else np.nan
        rec[f"{m}_trend"] = slope(v, yrs) if ok.sum() >= 3 else np.nan
        rec[f"{m}_nyears"] = int(ok.sum())
    prows.append(rec)
pooled = pd.DataFrame(prows).set_index('plotID')
print(f"TREE pooled structure: {len(pooled)} plots", flush=True)

# --- update lidar_pooled_predictors_26.csv: fill TREE structure columns ---
LID = os.path.join(D, "lidar_pooled_predictors_26.csv")
lid = pd.read_csv(LID)
struct_cols = [c for c in lid.columns if any(c == f"{m}_{s}" for m in STRUCT for s in ["mean","sd","trend","nyears"])]
assert set(struct_cols) == set(pooled.columns), f"col mismatch: {set(pooled.columns) ^ set(struct_cols)}"
shutil.copy2(LID, LID + ".bak2")
tree_mask = lid.siteID == SITE
before = lid.loc[tree_mask, "Canopy_Ht_mean"].notna().sum()
for pid, prow in pooled.iterrows():
    m = tree_mask & (lid.plotID == pid)
    for c in struct_cols:
        lid.loc[m, c] = prow[c]
after = lid.loc[tree_mask, "Canopy_Ht_mean"].notna().sum()
lid.to_csv(LID, index=False)
print(f"lidar_26 TREE Canopy_Ht_mean non-null: {before} -> {after} (backup {LID}.bak2)", flush=True)

# --- append TREE rows to per_year_v2_26.csv (for scripts using per-year, e.g. 158) ---
PY26 = os.path.join(D, "per_year_v2_26.csv")
if os.path.exists(PY26):
    p26 = pd.read_csv(PY26)
    p26 = p26[p26.siteID != SITE]  # drop any stale TREE
    add = py_out.reindex(columns=p26.columns)
    shutil.copy2(PY26, PY26 + ".bak2")
    pd.concat([p26, add], ignore_index=True).to_csv(PY26, index=False)
    print(f"per_year_v2_26: appended {len(add)} TREE rows (backup {PY26}.bak2)", flush=True)

# save TREE per_year for the record
py_out.to_csv("scripts_pipeline/_pipeline_state/per_year_tree_10m.csv", index=False)
print("DONE", flush=True)
