"""
Third productivity proxy: PML-V2 GPP (Penman-Monteith-Leuning; different model from MOD17).
Triangulates DHI vs MODIS-GPP vs PML-GPP, and re-tests the diversity hump vs PML GPP.
Usage: python 150_pml_gpp.py geedankook
Output: data/plot_pml_gpp.csv, results/dhi_gpp_triangulation.csv, figures/L06_gpp_triangulation.png
"""
import sys, os, numpy as np, pandas as pd, ee, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
PROJECT=sys.argv[1] if len(sys.argv)>1 else 'geedankook'
ee.Initialize(project=PROJECT); print("EE init:",PROJECT)
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
R=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\results"
F=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\figures"
FT={'ABBY':'evergreen','WREF':'evergreen','TEAK':'evergreen','RMNP':'evergreen','SOAP':'evergreen','JERC':'evergreen','OSBS':'evergreen','TALL':'evergreen','BART':'deciduous','BLAN':'deciduous','SCBI':'deciduous','SERC':'deciduous','MLBS':'deciduous','ORNL':'deciduous','GRSM':'deciduous','HARV':'mixed','UNDE':'mixed','STEI':'mixed','TREE':'mixed'}

pl=pd.read_csv(os.path.join(D,"plot_lonlat.csv"))
plots=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r.lon,r.lat]),{'plotID':r.plotID}) for r in pl.itertuples()])
PML=ee.ImageCollection("CAS/IGSNRR/PML/V2_v017").select('GPP')   # 8-day GPP gC/m2/d, 500m
rows=[]
for y in range(2016,2021):   # PML v017 through ~2020
    img=PML.filterDate(f"{y}-01-01",f"{y}-12-31").mean()  # mean 8-day GPP ~ annual productivity index
    fc=img.reduceRegions(plots,ee.Reducer.mean(),500).getInfo()
    for f in fc['features']:
        p=f['properties']; rows.append((p['plotID'],y,p.get('mean',p.get('GPP'))))
    print(f"  {y} done",flush=True)
g=pd.DataFrame(rows,columns=['plotID','year','pml']).groupby('plotID')['pml'].mean().reset_index().rename(columns={'pml':'pml_gpp'})
g.to_csv(os.path.join(D,"plot_pml_gpp.csv"),index=False); print("saved plot_pml_gpp.csv | n:",g.pml_gpp.notna().sum())

# merge all
df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv")).merge(g,on='plotID',how='left')
df=df.merge(pd.read_csv(os.path.join(D,"plot_modis_gpp.csv")),on='plotID',how='left')
df=df.merge(pd.read_csv(os.path.join(D,"plot_dhi_sentinel.csv"))[['plotID','dhi_cum_mean']],on='plotID',how='left')
d=df[(df.sample_coverage>=0.9)].copy(); d['ft']=d.siteID.map(FT)

# correlations among 3 proxies (site-level)
sm=d.groupby('siteID').agg(DHI=('dhi_cum_mean','median'),MODIS=('modis_gpp','median'),PML=('pml_gpp','median')).dropna()
print("\n=== 3 프록시 사이트 상관 ===")
for a,b in [('DHI','MODIS'),('DHI','PML'),('MODIS','PML')]:
    r=st.pearsonr(sm[a],sm[b]); print(f"  {a}-{b}: r={r[0]:+.2f} p={r[1]:.2g}")

# hump test vs each proxy (quad sign, +forest_type)
print("\n=== Hill q1 혹형(2차) 부호: 프록시별 (+forest_type 통제) ===")
res=[]
for x,lab in [('dhi_cum_mean','DHI'),('modis_gpp','MODIS GPP'),('pml_gpp','PML GPP')]:
    dd=d[[x,'Hill_q1','ft']].dropna().copy(); dd['zx']=(dd[x]-dd[x].mean())/dd[x].std()
    m=smf.ols("Hill_q1 ~ zx + I(zx**2) + C(ft)",dd).fit()
    k=[c for c in m.params.index if '** 2' in c][0]
    shape='혹형∩' if m.params[k]<0 else 'U자∪/단조'
    print(f"  {lab:10s} 2차 β={m.params[k]:+.3f} p={m.pvalues[k]:.2g}  -> {shape}")
    res.append(dict(proxy=lab,quad_beta=m.params[k],quad_p=m.pvalues[k],shape=shape))
pd.DataFrame(res).to_csv(os.path.join(R,"dhi_gpp_triangulation.csv"),index=False)

# figure: Hill q1 vs each proxy
fig,ax=plt.subplots(1,3,figsize=(18,5.2))
for a,(x,lab) in zip(ax,[('dhi_cum_mean','DHI누적 (녹색도)'),('modis_gpp','MODIS GPP'),('pml_gpp','PML GPP')]):
    dd=d[[x,'Hill_q1','siteID']].dropna()
    a.scatter(dd[x],dd.Hill_q1,s=12,alpha=.3,color='#00695c')
    gg=dd.groupby('siteID').agg(x=(x,'median'),y=('Hill_q1','mean')); a.scatter(gg.x,gg.y,s=45,color='#c62828',zorder=3)
    xx=np.linspace(dd[x].min(),dd[x].max(),50); cf=np.polyfit(dd[x],dd.Hill_q1,2); a.plot(xx,np.polyval(cf,xx),'k-',lw=2.4)
    a.set_xlabel(lab); a.set_ylabel('Hill q1'); a.set_title(lab); a.grid(alpha=.2)
fig.suptitle("L06. 다양성 혹형 삼각검증 — DHI만 혹형(∩), 실제 생산성(MODIS·PML)은 단조증가",fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(F,"L06_gpp_triangulation.png"),dpi=120,bbox_inches='tight'); plt.close()
print("saved L06_gpp_triangulation.png")
