"""
Obj2/Obj3 foundation — Sentinel-2 L2A intra-annual NDVI + NBR time series per NEON plot.
Source: Microsoft Planetary Computer STAC (no auth). 2016-2025.

Efficiency: process per ~8 km plot CLUSTER (so the read window stays small/fast even for
large sites). Outlier plots (georef errors >30 km from site median) dropped. For each scene
read bands ONCE over the cluster window, SCL-mask, sample every plot's 4x4 footprint.
Writes incremental per-site CSV (resumable: skips sites already done).

Bands: B04(red,10m), B08(nir,10m), B12(swir2,20m), SCL(20m). NDVI[Obj3 DHI], NBR[Obj2 disturb].
Output: NEON_v2/data/sentinel/s2_<SITE>.csv (plotID,date,year,doy,ndvi,nbr,npix)
"""
import os, sys, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from pyproj import Transformer
import pystac_client, planetary_computer, rasterio, rasterio.warp
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
from concurrent.futures import ThreadPoolExecutor
WORKERS=6; CELL=15000  # cluster cell size (m)
OUT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\sentinel"
os.makedirs(OUT,exist_ok=True)
YEARS=(2016,2025); CLOUD=30
SITES_ARG=sys.argv[1:] if len(sys.argv)>1 else None

ppy=pd.read_csv('E:/neon_lidar/vegetation_structure/vst_perplotperyear.csv',low_memory=False)
co=ppy.dropna(subset=['easting','northing','utmZone']).drop_duplicates('plotID')
cat=pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace)

def _rd(href,bbox,S,rs):
    with rasterio.open(href) as r:
        b=rasterio.warp.transform_bounds("EPSG:4326",r.crs,*bbox)
        return r.read(1,window=from_bounds(*b,r.transform),out_shape=S,resampling=rs).astype(np.float32)

def process_cluster(gc,tr):
    plon,plat=tr.transform(gc.easting.values,gc.northing.values); pids=gc.plotID.values
    bbox=[plon.min()-0.002,plat.min()-0.002,plon.max()+0.002,plat.max()+0.002]
    items=list(cat.search(collections=["sentinel-2-l2a"],bbox=bbox,
        datetime=f"{YEARS[0]}-01-01/{YEARS[1]}-12-31",query={"eo:cloud_cover":{"lt":CLOUD}}).items())
    def scene_work(it):
        out=[]
        try:
            a=it.assets
            with rasterio.open(a['B04'].href) as r:
                crs=r.crs; b=rasterio.warp.transform_bounds("EPSG:4326",crs,*bbox)
                win=from_bounds(*b,r.transform); tw=int(win.width); th=int(win.height)
                if tw<1 or th<1 or tw*th>3_000_000: return out
                S=(th,tw)
                red=r.read(1,window=win,out_shape=S,resampling=Resampling.bilinear).astype(np.float32)
                wt=rasterio.windows.transform(win,r.transform)
            nir=_rd(a['B08'].href,bbox,S,Resampling.bilinear)
            sw=_rd(a['B12'].href,bbox,S,Resampling.bilinear)
            scl=_rd(a['SCL'].href,bbox,S,Resampling.nearest)
            good=np.isin(scl,[4,5,6,7])
            ndvi=np.where(good,(nir-red)/(nir+red+1e-9),np.nan)
            nbr=np.where(good,(nir-sw)/(nir+sw+1e-9),np.nan)
            inv=~wt; xs,ys=rasterio.warp.transform("EPSG:4326",crs,list(plon),list(plat))
            dt=it.datetime
            for pid,xx,yy in zip(pids,xs,ys):
                col,row=inv*(xx,yy); col=int(col); row=int(row)
                r0,r1=max(0,row-2),min(th,row+2); c0,c1=max(0,col-2),min(tw,col+2)
                if r1<=r0 or c1<=c0: continue
                sub=ndvi[r0:r1,c0:c1]; n=np.isfinite(sub).sum()
                if n<1: continue
                out.append((pid,dt.date().isoformat(),dt.year,int(dt.timetuple().tm_yday),
                            round(float(np.nanmedian(sub)),4),
                            round(float(np.nanmedian(nbr[r0:r1,c0:c1])),4),int(n)))
        except Exception: pass
        return out
    rows=[]
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for res in ex.map(scene_work,items): rows.extend(res)
    return rows,len(items)

def process_site(site,g):
    fp=os.path.join(OUT,f"s2_{site}.csv")
    if os.path.exists(fp): print(f"{site}: done, skip",flush=True); return
    g=g.copy()
    ex,ny=g.easting.values,g.northing.values; mx,my=np.median(ex),np.median(ny)
    g=g[(np.abs(ex-mx)<30000)&(np.abs(ny-my)<30000)]
    if len(g)<5: print(f"{site}: <5 plots, skip",flush=True); return
    tr=Transformer.from_crs(f"EPSG:326{int(str(g.utmZone.iloc[0])[:2]):02d}","EPSG:4326",always_xy=True)
    g['cx']=(g.easting//CELL).astype(int); g['cy']=(g.northing//CELL).astype(int)
    clusters=list(g.groupby(['cx','cy']))
    print(f"{site}: {len(g)} plots in {len(clusters)} clusters",flush=True)
    allrows=[]
    for ci,((cx,cy),gc) in enumerate(clusters):
        rows,nsc=process_cluster(gc,tr); allrows+=rows
        print(f"  {site} cluster {ci+1}/{len(clusters)} ({len(gc)}p,{nsc}sc) -> {len(rows)} obs, total {len(allrows)}",flush=True)
    pd.DataFrame(allrows,columns=['plotID','date','year','doy','ndvi','nbr','npix']).to_csv(fp,index=False)
    print(f"{site}: wrote {len(allrows)} obs -> {fp}",flush=True)

sites=sorted(co.siteID.unique()) if not SITES_ARG else SITES_ARG
for s in sites:
    g=co[co.siteID==s]
    if len(g)>=5: process_site(s,g)
print("ALL DONE",flush=True)
