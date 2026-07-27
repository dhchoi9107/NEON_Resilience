# Land-use heterogeneity for 26 sites via Microsoft Planetary Computer (ESA WorldCover). No GEE.
import os, numpy as np, pandas as pd, rasterio, rasterio.warp
from rasterio.windows import from_bounds
import pystac_client, planetary_computer
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
pl=pd.read_csv(os.path.join(D,"plot_lonlat_26.csv"))
cat=pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=planetary_computer.sign_inplace)
HALF=0.007
def het(lon,lat):
    bbox=[lon-HALF,lat-HALF,lon+HALF,lat+HALF]
    it=list(cat.search(collections=["esa-worldcover"],bbox=bbox).items())
    if not it: return None
    with rasterio.open(it[0].assets['map'].href) as r:
        b=rasterio.warp.transform_bounds("EPSG:4326",r.crs,*bbox)
        arr=r.read(1,window=from_bounds(*b,r.transform))
    a=arr[arr>0]
    if a.size<10: return None
    u,c=np.unique(a,return_counts=True); p=c/c.sum()
    edges=((arr[:,1:]!=arr[:,:-1]).sum()+(arr[1:,:]!=arr[:-1,:]).sum())/(2*arr.size)
    return dict(lc_shannon=float(-(p*np.log(p)).sum()),lc_nclass=int(len(u)),
                lc_dominant=float(p.max()),lc_edge=float(edges),lc_forest_frac=float((a==10).mean()))
rows=[]
for i,r in enumerate(pl.itertuples()):
    try:
        h=het(r.lon,r.lat)
        if h: h['plotID']=r.plotID; h['siteID']=r.siteID; rows.append(h)
    except Exception as e:
        if i%200==0: print(f"  err {r.plotID}: {type(e).__name__}",flush=True)
    if i%150==0: print(f"  {i}/{len(pl)}",flush=True)
lu=pd.DataFrame(rows); lu.to_csv(os.path.join(D,"plot_landuse_het_26.csv"),index=False)
new=['DELA','LENO','UKFS','YELL','BONA','DEJU','HEAL']
print(f"DONE plot_landuse_het_26.csv: {len(lu)} plots | new-site cov {lu[lu.siteID.isin(new)].siteID.nunique()}/7")
