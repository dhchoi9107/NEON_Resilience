"""
Aggregate NEW-site 1m FSD plot rasters -> 13 plot-year structure metrics, merge with VI,
producing per_year rows for the 7 new sites in the per_year_v2 format.
Out: scripts_pipeline/_pipeline_state/per_year_newsites.csv
FSD raster: E:/neon_lidar/structural_diversity_1m_plots/{site}/{year}_{site}_{plot}_FSD_1m.tif
"""
import os, glob, re, numpy as np, pandas as pd, rasterio
NEW = ["DELA", "LENO", "UKFS", "YELL", "BONA", "DEJU", "HEAL"]
FSD_BANDS = ['rumple','top_rugosity','mean_max_canopy_ht','max_canopy_ht','deepgap_fraction','meanH','vert_sd','vertCV',
 'mean_sd','sd_sd','GC','GFP','VCI','q25','q50','q75','q95','HeightRatio','FHD','LAI','LAI_subcanopy']
SKEEP = {'Canopy_Ht':'mean_max_canopy_ht','Max_Ht':'max_canopy_ht','Rumple':'rumple','Rugosity':'top_rugosity',
 'Deep_Gap':'deepgap_fraction','Vert_SD':'vert_sd','Vert_CV':'vertCV','Gini':'GC','VCI':'VCI','FHD':'FHD',
 'LAI':'LAI','Q95':'q95','Ht_Ratio':'HeightRatio'}
BASE = "E:/neon_lidar/structural_diversity_1m_plots"
bidx = {b: i for i, b in enumerate(FSD_BANDS)}

rows = []
for site in NEW:
    d = os.path.join(BASE, site)
    if not os.path.isdir(d):
        print(f"  {site}: no FSD dir"); continue
    n = 0
    for f in glob.glob(os.path.join(d, "*_FSD_1m.tif")):
        m = re.match(r"(\d{4})_([A-Z]{4})_(\w+)_FSD_1m\.tif", os.path.basename(f))
        if not m: continue
        year, s, plotnum = int(m.group(1)), m.group(2), m.group(3)
        plotID = f"{s}_{plotnum}"
        try:
            with rasterio.open(f) as src:
                arr = src.read()  # (bands, H, W)
        except Exception as e:
            print(f"  read fail {f}: {e}"); continue
        rec = {"plotID": plotID, "siteID": s, "year": year}
        with np.errstate(all="ignore"):
            for col, band in SKEEP.items():
                rec[col] = float(np.nanmean(arr[bidx[band]])) if band in bidx else np.nan
        rows.append(rec); n += 1
    print(f"  {site}: {n} plot-year FSD rasters aggregated", flush=True)

fsd = pd.DataFrame(rows)
# merge VI (same plot-year)
vi = pd.read_csv("NEON_v2/data/plot_vi_neon_brdf.csv")
vic = vi[["plotID", "year", "NDVI_mean", "EVI_mean", "ARVI_mean", "SAVI_mean"]].rename(
    columns={"NDVI_mean": "NDVI", "EVI_mean": "EVI", "ARVI_mean": "ARVI", "SAVI_mean": "SAVI"})
out = fsd.merge(vic, on=["plotID", "year"], how="left")
OUT = "scripts_pipeline/_pipeline_state/per_year_newsites.csv"
out.to_csv(OUT, index=False)
print(f"\nNEW per_year: {len(out)} plot-years, {out.siteID.nunique()} sites, {out.plotID.nunique()} plots")
print(f"  structure non-null: Canopy_Ht {out.Canopy_Ht.notna().sum()}, VCI {out.VCI.notna().sum()}, LAI {out.LAI.notna().sum()}")
print(f"  VI matched: NDVI {out.NDVI.notna().sum()}/{len(out)}")
print(f"  -> {OUT}")
print(out.groupby('siteID').agg(plot_years=('plotID','size'), plots=('plotID','nunique'), yrs=('year','nunique')).to_string())
