"""
Validate DHI(cumulative) as a productivity proxy against MODIS GPP (MOD17A3HGF, annual GPP).
Extract MODIS annual GPP at NEON plots (2016-2023, server-side GEE), average, aggregate to site,
and correlate site-mean DHI_cum vs site-mean MODIS GPP.
Usage: python 149_modis_gpp.py geedankook
Output: data/plot_modis_gpp.csv, results/dhi_gpp_validation.csv
"""
import sys, os, numpy as np, pandas as pd, ee
PROJECT=sys.argv[1] if len(sys.argv)>1 else 'geedankook'
ee.Initialize(project=PROJECT); print("EE init:",PROJECT)
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
R=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\results"

pl=pd.read_csv(os.path.join(D,"plot_lonlat.csv"))
feats=[ee.Feature(ee.Geometry.Point([r.lon,r.lat]),{'plotID':r.plotID}) for r in pl.itertuples()]
plots=ee.FeatureCollection(feats); print("plots:",len(feats))

GPP=ee.ImageCollection("MODIS/061/MOD17A3HGF")   # annual GPP, band Gpp, scale 1e-4 kgC/m2/yr
def gpp_year(y):
    img=GPP.filterDate(f"{y}-01-01",f"{y}-12-31").first().select('Gpp')
    return img.updateMask(img.lt(30000)).multiply(0.0001).rename('gpp')  # mask fill, scale

rows=[]
for y in range(2016,2024):
    fc=gpp_year(y).reduceRegions(plots,ee.Reducer.mean(),500).getInfo()
    for f in fc['features']:
        p=f['properties']; rows.append((p['plotID'],y,p.get('mean',p.get('gpp'))))
    print(f"  {y} done",flush=True)
g=pd.DataFrame(rows,columns=['plotID','year','gpp'])
gm=g.groupby('plotID')['gpp'].mean().reset_index().rename(columns={'gpp':'modis_gpp'})
gm.to_csv(os.path.join(D,"plot_modis_gpp.csv"),index=False)
print("saved plot_modis_gpp.csv | plots with GPP:",gm.modis_gpp.notna().sum())

# correlate with DHI at plot and site level
import scipy.stats as st
dhi=pd.read_csv(os.path.join(D,"plot_dhi_sentinel.csv"))[['plotID','dhi_cum_mean']]
m=gm.merge(dhi,on='plotID').merge(pl[['plotID','siteID']],on='plotID').dropna()
rp=st.pearsonr(m.dhi_cum_mean,m.modis_gpp)
site=m.groupby('siteID').agg(dhi=('dhi_cum_mean','median'),gpp=('modis_gpp','median')).dropna()
rs=st.pearsonr(site.dhi,site.gpp)
print(f"\n=== DHI_cum vs MODIS GPP ===")
print(f"plot-level  r={rp[0]:+.3f} p={rp[1]:.2g} (n={len(m)})")
print(f"site-level  r={rs[0]:+.3f} p={rs[1]:.2g} (n={len(site)})")
pd.DataFrame([dict(level='plot',r=rp[0],p=rp[1],n=len(m)),
              dict(level='site',r=rs[0],p=rs[1],n=len(site))]).to_csv(os.path.join(R,"dhi_gpp_validation.csv"),index=False)
print("\n사이트별 DHI / MODIS GPP:"); print(site.round(3).sort_values('gpp').to_string())
