"""
Obj 2 (self-detect, GEE) — annual growing-season NBR per plot, server-side.
Detect disturbance = largest year-over-year NBR drop (LandTrendr-lite).
Usage: python 61_gee_nbr_disturb.py geedankook
Output: plot_disturb_s2.csv (plotID,s2_dist,s2_dist_year,s2_dist_mag,s2_recovery,nbr_trend)
"""
import sys, numpy as np, pandas as pd, ee
PROJECT=sys.argv[1] if len(sys.argv)>1 else "geedankook"
ee.Initialize(project=PROJECT); print("EE init",PROJECT)
DROP=-0.08
pl=pd.read_csv(r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_lonlat.csv")
plots=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r.lon,r.lat]).buffer(20),{'plotID':r.plotID}) for r in pl.itertuples()])
S2=ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"); CS=ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")

def nbr_year(y):
    # growing season Jun-Sep median NBR (peak canopy)
    s=ee.Date.fromYMD(y,6,1); e=ee.Date.fromYMD(y,9,30)
    col=(S2.filterDate(s,e).filterBounds(plots).linkCollection(CS,['cs_cdf'])
           .map(lambda im: im.updateMask(im.select('cs_cdf').gte(0.6))
                 .normalizedDifference(['B8','B12']).rename('nbr')))
    return col.median().rename('nbr')

rows=[]
for y in range(2016,2026):
    fc=nbr_year(y).reduceRegions(plots,ee.Reducer.mean(),20).getInfo()
    for f in fc['features']:
        rows.append((f['properties']['plotID'],y,f['properties'].get('mean')))
    print(f"  {y}: {len(fc['features'])}",flush=True)
ts=pd.DataFrame(rows,columns=['plotID','year','nbr']).dropna(subset=['nbr'])
ts.to_csv(r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_nbr_annual.csv",index=False)
print("saved annual NBR series:",ts.shape)

out=[]
for pid,g in ts.groupby('plotID'):
    g=g.sort_values('year'); yrs=g['year'].values; v=g['nbr'].values
    if len(g)<4: out.append((pid,0,0,np.nan,np.nan,np.nan)); continue
    d=np.diff(v); i=int(np.argmin(d)); mag=float(d[i]); dyear=int(yrs[i+1])
    post=g[g['year']>=dyear]; rec=float(np.polyfit(post['year'],post['nbr'],1)[0]) if len(post)>=3 else np.nan
    trend=float(np.polyfit(yrs,v,1)[0])
    out.append((pid,int(mag<DROP),dyear if mag<DROP else 0,round(mag,3),
                round(rec,4) if rec==rec else np.nan,round(trend,4)))
o=pd.DataFrame(out,columns=['plotID','s2_dist','s2_dist_year','s2_dist_mag','s2_recovery','nbr_trend'])
o.to_csv(r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_disturb_s2.csv",index=False)
print(f"DONE plots {len(o)} | disturbance flagged {o.s2_dist.sum()}")
print(o[o.s2_dist==1]['s2_dist_year'].value_counts().sort_index().to_string())
