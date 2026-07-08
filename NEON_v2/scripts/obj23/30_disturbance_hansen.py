"""
Obj 2 (authoritative) — Hansen Global Forest Change v1.11 (2000-2023) tree-cover loss per plot.
Reads lossyear granules remotely via /vsicurl (only the 40x40m plot window) -> no full download.
lossyear: 0=no loss, 1..23 = loss in 2001..2023.
Output: plot_hansen_loss.csv (plotID, hansen_loss, hansen_lossyear, hansen_loss_frac, treecover2000)
"""
import os, warnings, numpy as np, pandas as pd; warnings.filterwarnings('ignore')
import rasterio
from rasterio.windows import from_bounds
from pyproj import Transformer
BASE="https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11"
OUT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_hansen_loss.csv"

ppy=pd.read_csv('E:/neon_lidar/vegetation_structure/vst_perplotperyear.csv',low_memory=False)
co=ppy.dropna(subset=['easting','northing','utmZone']).drop_duplicates('plotID')
# UTM -> lon/lat per plot
def lonlat(r):
    z=int(str(r.utmZone)[:2]); tr=Transformer.from_crs(f"EPSG:326{z:02d}","EPSG:4326",always_xy=True)
    return tr.transform(r.easting,r.northing)
co[['lon','lat']]=co.apply(lambda r: pd.Series(lonlat(r)),axis=1)

def granule(lat,lon):
    latt=int(np.ceil(lat/10.0)*10); lont=int(np.floor(lon/10.0)*10)
    return f"{latt}N_{abs(lont):03d}W"
co['gran']=[granule(la,lo) for la,lo in zip(co.lat,co.lon)]
print("granules:",sorted(co.gran.unique()))

half=0.0003  # ~33 m
rows=[]
for gr,g in co.groupby('gran'):
    ly=f"/vsicurl/{BASE}/Hansen_GFC-2023-v1.11_lossyear_{gr}.tif"
    tc=f"/vsicurl/{BASE}/Hansen_GFC-2023-v1.11_treecover2000_{gr}.tif"
    try:
        rly=rasterio.open(ly); rtc=rasterio.open(tc)
    except Exception as e:
        print(f"{gr}: open err {str(e)[:60]}"); continue
    print(f"{gr}: {len(g)} plots",flush=True)
    for _,p in g.iterrows():
        try:
            w=from_bounds(p.lon-half,p.lat-half,p.lon+half,p.lat+half,rly.transform)
            a=rly.read(1,window=w).astype(int)
            wt=from_bounds(p.lon-half,p.lat-half,p.lon+half,p.lat+half,rtc.transform)
            t=rtc.read(1,window=wt).astype(float)
            loss=a[a>0]
            rows.append((p.plotID, int((a>0).any()),
                         int(2000+loss.max()) if loss.size else 0,
                         round(float((a>0).mean()),3), round(float(np.nanmean(t)),1)))
        except Exception:
            rows.append((p.plotID,0,0,0.0,np.nan))
    rly.close(); rtc.close()
out=pd.DataFrame(rows,columns=['plotID','hansen_loss','hansen_lossyear','hansen_loss_frac','treecover2000'])
out.to_csv(OUT,index=False)
print(f"DONE -> {OUT}\nplots {len(out)} | with loss: {out.hansen_loss.sum()} | lossyear dist:")
print(out[out.hansen_loss==1]['hansen_lossyear'].value_counts().sort_index().to_string())
