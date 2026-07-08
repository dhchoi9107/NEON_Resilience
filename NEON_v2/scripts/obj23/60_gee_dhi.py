"""
Obj 3 (GEE) — Sentinel-2 Dynamic Habitat Indices, computed SERVER-SIDE on Earth Engine.
Pulls only the small per-plot-year result table (fast; no scene downloads).
Usage:  python 60_gee_dhi.py <EE_PROJECT_ID>
DHI per plot-year (from monthly NDVI medians): cumulative(mean) / minimum / variation(CV).
Output: NEON_v2/data/plot_dhi_gee.csv  (plot-level mean + trend per component)
"""
import sys, os, numpy as np, pandas as pd, ee
PROJECT=sys.argv[1] if len(sys.argv)>1 else None
ee.Initialize(project=PROJECT)
print("EE initialized, project:",PROJECT)

pl=pd.read_csv(r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_lonlat.csv")
feats=[ee.Feature(ee.Geometry.Point([r.lon,r.lat]).buffer(20),{'plotID':r.plotID}) for r in pl.itertuples()]
plots=ee.FeatureCollection(feats)
print("plots:",len(feats))

S2=ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
CS=ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
def ndvi_monthly(year):
    """Reduce the year's cloud-masked NDVI collection directly (robust to empty winter months).
       cumulative=mean, minimum=10th percentile (noise-robust low), variation=stdDev/mean."""
    start=ee.Date.fromYMD(year,1,1); end=start.advance(1,'year')
    col=(S2.filterDate(start,end).filterBounds(plots)
           .linkCollection(CS,['cs_cdf'])
           .map(lambda im: im.updateMask(im.select('cs_cdf').gte(0.6))
                 .normalizedDifference(['B8','B4']).rename('ndvi')))
    cum=col.mean().rename('dhi_cum')
    mn=col.reduce(ee.Reducer.percentile([10])).rename('dhi_min')
    sd=col.reduce(ee.Reducer.stdDev()).rename('sd')
    var=sd.divide(cum.max(ee.Image(1e-6))).rename('dhi_var')
    return cum.addBands(mn).addBands(var)

rows=[]
for y in range(2016,2026):
    img=ndvi_monthly(y)
    fc=img.reduceRegions(plots,ee.Reducer.mean(),20).getInfo()
    for f in fc['features']:
        p=f['properties']
        rows.append((p['plotID'],y,p.get('dhi_cum'),p.get('dhi_min'),p.get('dhi_var')))
    print(f"  {y}: pulled {len(fc['features'])} plots",flush=True)

py=pd.DataFrame(rows,columns=['plotID','year','dhi_cum','dhi_min','dhi_var']).dropna(subset=['dhi_cum'])
def slope(s,x): return np.polyfit(x,s,1)[0] if len(x)>=3 else np.nan
out=[]
for pid,g in py.groupby('plotID'):
    g=g.sort_values('year'); yr=g['year'].values; rec={'plotID':pid,'dhi_nyears':len(g)}
    for c in ['dhi_cum','dhi_min','dhi_var']:
        rec[f'{c}_mean']=float(g[c].mean()); rec[f'{c}_trend']=float(slope(g[c].values,yr))
    out.append(rec)
out=pd.DataFrame(out)
out.to_csv(r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_dhi_sentinel.csv",index=False)
print(f"DONE plots {len(out)}")
print(out[['dhi_cum_mean','dhi_min_mean','dhi_var_mean']].describe().round(3).to_string())
