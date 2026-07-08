"""
Obj 2 gap fill — LAND-USE HETEROGENEITY (proposal: 'disturbance events AND land-use heterogeneity').
Per plot (750 m buffer) from ESA WorldCover 10 m:
  compositional: LC Shannon, n_classes, dominant fraction
  configurational/fragmentation: edge density (class-boundary pixel fraction)
+ site-level configurational surface metrics (heterogeneity_all.csv: Sa/Sq at 500 m).
Then: (1) diversity ~ heterogeneity; (2) does heterogeneity MODERATE the RS<->taxonomic coupling?
Output: plot_landuse_het.csv + FINAL_v2_full.csv + results/obj2_heterogeneity.csv
"""
import os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, statsmodels.formula.api as smf
import pystac_client, planetary_computer, rasterio, rasterio.warp
from rasterio.windows import from_bounds
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
R=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\results"
def z(s): return (s-s.mean())/s.std()

pl=pd.read_csv(os.path.join(D,"plot_lonlat.csv"))
cat=pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=planetary_computer.sign_inplace)
HALF=0.007
def het(lon,lat):
    bbox=[lon-HALF,lat-HALF,lon+HALF,lat+HALF]
    it=list(cat.search(collections=["esa-worldcover"],bbox=bbox).items())
    if not it: return None
    with rasterio.open(it[0].assets['map'].href) as r:
        b=rasterio.warp.transform_bounds("EPSG:4326",r.crs,*bbox)
        a=r.read(1,window=from_bounds(*b,r.transform))
    a=a[a>0]
    if a.size<10: return None
    # reload as 2D for edge density
    with rasterio.open(it[0].assets['map'].href) as r:
        b=rasterio.warp.transform_bounds("EPSG:4326",r.crs,*bbox)
        arr=r.read(1,window=from_bounds(*b,r.transform))
    u,c=np.unique(a,return_counts=True); p=c/c.sum()
    edges=((arr[:,1:]!=arr[:,:-1]).sum()+(arr[1:,:]!=arr[:-1,:]).sum())/(2*arr.size)
    return dict(lc_shannon=float(-(p*np.log(p)).sum()),lc_nclass=int(len(u)),
                lc_dominant=float(p.max()),lc_edge=float(edges),lc_forest_frac=float((a==10).mean()))

rows=[]
for i,r in enumerate(pl.itertuples()):
    h=het(r.lon,r.lat)
    if h: h['plotID']=r.plotID; rows.append(h)
    if i%150==0: print(f"  {i}/{len(pl)}",flush=True)
lu=pd.DataFrame(rows); lu.to_csv(os.path.join(D,"plot_landuse_het.csv"),index=False)
print(f"land-use het: {len(lu)} plots | LC Shannon 중앙 {lu.lc_shannon.median():.2f}, edge 중앙 {lu.lc_edge.median():.3f}")

# site-level configurational heterogeneity (500m surface metrics)
sh=pd.read_csv("E:/neon_lidar/env_heterogeneity/heterogeneity_all.csv")
sh=sh.groupby('siteID')[['Sa_500m','Sq_500m','Sa_1000m']].mean().reset_index()

base=pd.read_csv(os.path.join(D,"FINAL_v2_specdiv.csv"))
m=base.merge(lu,on='plotID',how='left').merge(sh,on='siteID',how='left')
m.to_csv(os.path.join(D,"FINAL_v2_full.csv"),index=False)
mc=m[m['sample_coverage']>=0.9].copy()
HET=['lc_shannon','lc_nclass','lc_edge','lc_forest_frac','Sa_500m','Sa_1000m']
print(f"\ncov>=0.9 {len(mc)} | with land-use het {mc['lc_shannon'].notna().sum()}")

RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'Turnover','LCBD_nestedness_rare':'Nestedness'}
for r in RESP: mc[f'z_{r}']=z(mc[r])
rows=[]
# (1) diversity ~ heterogeneity
for resp,rl in RESP.items():
    for hcol in HET:
        d=mc[[f'z_{resp}',hcol,'siteID','domain']].dropna()
        if len(d)<40 or d['siteID'].nunique()<3: continue
        d=d.copy(); d['zx']=z(d[hcol])
        try:
            mm=smf.mixedlm(f"z_{resp} ~ zx + C(domain)",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
            if mm.converged: rows.append(dict(test='div~het',response=rl,term=hcol,beta=mm.fe_params['zx'],p=mm.pvalues['zx']))
        except Exception: pass
# (2) heterogeneity moderates RS<->taxonomic coupling
for resp,rl in RESP.items():
    for rs in ['SAVI_mean','FHD_mean','specdiv_rao_q']:
        for hcol in ['lc_shannon','lc_edge']:
            d=mc[[f'z_{resp}',rs,hcol,'siteID','domain']].dropna()
            if len(d)<50: continue
            d=d.copy(); d['zx']=z(d[rs]); d['zh']=z(d[hcol])
            try:
                mm=smf.mixedlm(f"z_{resp} ~ zx*zh + C(domain)",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
                if mm.converged and 'zx:zh' in mm.fe_params:
                    rows.append(dict(test='moderation',response=rl,term=f"{rs}×{hcol}",beta=mm.fe_params['zx:zh'],p=mm.pvalues['zx:zh']))
            except Exception: pass
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"obj2_heterogeneity.csv"),index=False)
print("\n=== Obj2: (1) 다양성 ~ 토지이용 이질성 (p<0.05) ===")
for _,x in res[(res.test=='div~het')&(res.p<0.05)].sort_values('p').iterrows():
    print(f"  {x['response']:10s} ~ {x['term']:14s} beta={x['beta']:+.3f} p={x['p']:.3g}")
print("\n=== Obj2: (2) 이질성이 RS↔종 관계 조절 (p<0.1) ===")
for _,x in res[(res.test=='moderation')&(res.p<0.1)].sort_values('p').iterrows():
    print(f"  {x['response']:10s} {x['term']:22s} interaction={x['beta']:+.3f} p={x['p']:.3g}")
print("DONE")
