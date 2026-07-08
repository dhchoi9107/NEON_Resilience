"""
Obj 2 (authoritative) — MTBS fire perimeters per plot (1984-present, US fire severity program).
Spatial-join NEON plot points to MTBS burned-area polygons.
Output: plot_mtbs_fire.csv (plotID, fire, fire_year, fire_type, n_fires)
"""
import os, glob, zipfile, warnings, numpy as np, pandas as pd; warnings.filterwarnings('ignore')
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Transformer
DD=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\disturbance"
OUT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_mtbs_fire.csv"

zp=os.path.join(DD,"mtbs_perims.zip")
if not glob.glob(os.path.join(DD,"*.shp")):
    print("unzipping MTBS...");
    with zipfile.ZipFile(zp) as z: z.extractall(DD)
shp=glob.glob(os.path.join(DD,"*.shp"))[0]
print("shapefile:",os.path.basename(shp))
fires=gpd.read_file(shp)
print("MTBS perimeters:",len(fires),"| cols:",[c for c in fires.columns][:12])
fires=fires.to_crs("EPSG:4326")
# date/type columns (handle naming)
dcol=next((c for c in fires.columns if c.lower() in ('ig_date','igdate','date')),None)
tcol=next((c for c in fires.columns if 'incid_type' in c.lower() or c.lower()=='type'),None)
fires['fyear']=pd.to_datetime(fires[dcol],errors='coerce').dt.year if dcol else np.nan

# plot points
ppy=pd.read_csv('E:/neon_lidar/vegetation_structure/vst_perplotperyear.csv',low_memory=False)
co=ppy.dropna(subset=['easting','northing','utmZone']).drop_duplicates('plotID')
def ll(r):
    z=int(str(r.utmZone)[:2]); tr=Transformer.from_crs(f"EPSG:326{z:02d}","EPSG:4326",always_xy=True)
    return tr.transform(r.easting,r.northing)
co[['lon','lat']]=co.apply(lambda r: pd.Series(ll(r)),axis=1)
pts=gpd.GeoDataFrame(co[['plotID']].copy(),geometry=[Point(x,y) for x,y in zip(co.lon,co.lat)],crs="EPSG:4326")

j=gpd.sjoin(pts,fires[['geometry','fyear']+([tcol] if tcol else [])],how='left',predicate='within')
rows=[]
for pid,g in j.groupby('plotID'):
    gg=g.dropna(subset=['fyear'])
    if len(gg):
        recent=gg.loc[gg['fyear'].idxmax()]
        rows.append((pid,1,int(gg['fyear'].max()),str(recent[tcol]) if tcol else 'fire',len(gg)))
    else:
        rows.append((pid,0,0,'',0))
out=pd.DataFrame(rows,columns=['plotID','fire','fire_year','fire_type','n_fires'])
out.to_csv(OUT,index=False)
print(f"DONE -> {OUT}\nplots {len(out)} | burned: {out.fire.sum()}")
print(out[out.fire==1]['fire_year'].value_counts().sort_index().to_string())
