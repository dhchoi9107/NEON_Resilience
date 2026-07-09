"""
Extract Global 30m Forest Age (Zhang et al. 2025) at NEON plots.
Input: E:/neon_lidar/stand_age/{NF,PF}_NorthAmerica.7z  (natural / planted)
Steps: extract -> glob GeoTIFF tiles -> assign each plot to overlapping tile(s)
       (reprojecting plot lon/lat to each tile CRS) -> sample age.
Age = years since last spectral change point. NF preferred; PF if NF nodata.
Output: data/plot_stand_age_30m.csv  (plotID, stand_age_30m, forest_type_30m)
"""
import os, glob, numpy as np, pandas as pd, rasterio, py7zr, time
from rasterio.warp import transform as warp_transform
AGE_DIR=r"E:\neon_lidar\stand_age"
P=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
plots=pd.read_csv(os.path.join(P,"plot_lonlat.csv"))

def ensure_extracted(archive, subdir):
    outdir=os.path.join(AGE_DIR, subdir)
    if os.path.isdir(outdir) and glob.glob(os.path.join(outdir,"**","*.tif"),recursive=True):
        return outdir
    os.makedirs(outdir, exist_ok=True)
    print(f"extracting {archive} ...", flush=True); t=time.time()
    with py7zr.SevenZipFile(os.path.join(AGE_DIR,archive),'r') as z: z.extractall(outdir)
    print(f"  extracted {time.time()-t:.0f}s", flush=True)
    return outdir

def sample_set(outdir, label):
    tifs=glob.glob(os.path.join(outdir,"**","*.tif"),recursive=True)
    print(f"[{label}] {len(tifs)} tiles", flush=True)
    age=np.full(len(plots), np.nan)
    lon=plots.lon.values; lat=plots.lat.values
    for tp in tifs:
        try: src=rasterio.open(tp)
        except Exception: continue
        with src:
            xs,ys=warp_transform("EPSG:4326", src.crs, lon.tolist(), lat.tolist())
            xs=np.array(xs); ys=np.array(ys)
            l,b,r,t=src.bounds
            m=(xs>=l)&(xs<=r)&(ys>=b)&(ys<=t)&~np.isfinite(age) if False else (xs>=l)&(xs<=r)&(ys>=b)&(ys<=t)
            idx=np.where(m)[0]
            if not len(idx): continue
            nod=src.nodata
            for i,v in zip(idx, src.sample(list(zip(xs[idx],ys[idx])))):
                val=v[0]
                if nod is not None and val==nod: continue
                if val is None or (isinstance(val,(int,float)) and val<0): continue
                if np.isnan(age[i]): age[i]=float(val)
    return age

nf=ensure_extracted("NF_NorthAmerica.7z","NF"); age_nf=sample_set(nf,"NF")
pf=ensure_extracted("PF_NorthAmerica.7z","PF"); age_pf=sample_set(pf,"PF")
age=np.where(np.isfinite(age_nf), age_nf, age_pf)
ftype=np.where(np.isfinite(age_nf),"natural", np.where(np.isfinite(age_pf),"planted","none"))
plots['stand_age_30m']=age; plots['forest_type_30m']=ftype
out=os.path.join(P,"plot_stand_age_30m.csv"); plots.to_csv(out,index=False)
n=np.isfinite(age).sum()
print(f"\nDONE -> {out} | {n}/{len(plots)} plots with 30m age (NF {np.isfinite(age_nf).sum()}, PF-only {np.isfinite(age_pf).sum()})")
print(plots.groupby('siteID')['stand_age_30m'].median().round(0).to_string())
