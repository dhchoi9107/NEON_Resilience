"""
Does the diversity~age HUMP survive controlling FOREST TYPE (biogeographic confound)?
Controls tested (plot-level mixed model, site random):
  A raw:  y ~ za + za^2
  B +domain:       + C(domain)
  C +forest_type:  + C(forest_type)   [documented NEON dominant type, 3 classes]
  D +evergreenness:+ dhi_min_mean      [objective: winter-NDVI floor = evergreen index]
Reports the quadratic (hump) p under each. If the hump vanishes -> it was the confound.
Output: results/foresttype_hump.csv, figures/N10_age_diversity_byforesttype.png
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"results"),os.path.join(ROOT,"figures")

# documented NEON dominant forest type (site-level)
FT={'ABBY':'evergreen','WREF':'evergreen','TEAK':'evergreen','RMNP':'evergreen','SOAP':'evergreen',
    'JERC':'evergreen','OSBS':'evergreen','TALL':'evergreen',
    'BART':'deciduous','BLAN':'deciduous','SCBI':'deciduous','SERC':'deciduous','MLBS':'deciduous',
    'ORNL':'deciduous','GRSM':'deciduous',
    'HARV':'mixed','UNDE':'mixed','STEI':'mixed','TREE':'mixed'}
ftc={'evergreen':'#1b5e20','deciduous':'#e65100','mixed':'#6a1b9a'}

df=pd.read_csv(f"{D}/FINAL_v2_pooled.csv").merge(
   pd.read_csv(f"{D}/plot_stand_age_gami.csv")[['plotID','stand_age_gami']],on='plotID',how='left').merge(
   pd.read_csv(f"{D}/plot_dhi_sentinel.csv")[['plotID','dhi_min_mean']],on='plotID',how='left')
df=df[(df.sample_coverage>=0.9)&df.stand_age_gami.notna()].copy()
df['forest_type']=df.siteID.map(FT)
df['za']=(df.stand_age_gami-df.stand_age_gami.mean())/df.stand_age_gami.std()

# validate: evergreenness (dhi_min) by forest type
print("=== 검증: forest type별 DHI_min (상록성) — evergreen > deciduous 여야 ===")
print(df.groupby('forest_type')['dhi_min_mean'].median().round(3).to_string())
print("사이트별 dhi_min:", df.groupby('siteID')['dhi_min_mean'].median().round(2).sort_values().to_dict())

DIV={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'turnover','LCBD_nestedness_rare':'nestedness'}
def quadp(formula,d,re=True):
    try:
        if re: m=smf.mixedlm(formula,d,groups=d.siteID).fit(reml=True,method='lbfgs')
        else:  m=smf.ols(formula,d).fit()
        for k in m.pvalues.index:
            if 'za ** 2' in k: return m.pvalues[k],m.params[k]
    except Exception: pass
    return np.nan,np.nan

rows=[]
print("\n=== 혹형(2차) p — 통제별 ===")
print(f"{'index':11s} {'A raw':>10s} {'B +domain':>12s} {'C +forest_type':>16s} {'D +evergreen':>14s}")
for c,l in DIV.items():
    d=df[[c,'za','siteID','domain','forest_type','dhi_min_mean']].dropna()
    pa,_=quadp(f"{c} ~ za + I(za**2)",d,re=False)
    pb,_=quadp(f"{c} ~ za + I(za**2) + C(domain)",d)
    pc,_=quadp(f"{c} ~ za + I(za**2) + C(forest_type)",d)
    pd_,_=quadp(f"{c} ~ za + I(za**2) + dhi_min_mean",d)
    rows.append(dict(index=l,A_raw=pa,B_domain=pb,C_foresttype=pc,D_evergreen=pd_,n=len(d)))
    print(f"{l:11s} {pa:10.2g} {pb:12.2g} {pc:16.2g} {pd_:14.2g}")
pd.DataFrame(rows).to_csv(f"{R}/foresttype_hump.csv",index=False)

# figure: age vs Hill q1 colored by forest type, per-type quadratic
fig,axes=plt.subplots(1,2,figsize=(15,6))
for ax,c,l in [(axes[0],'Hill_q1','Hill q1 (알파)'),(axes[1],'LCBD_turnover_rare','LCBD turnover (베타)')]:
    d=df[[c,'stand_age_gami','forest_type']].dropna()
    for ft in ['evergreen','deciduous','mixed']:
        dd=d[d.forest_type==ft]
        ax.scatter(dd.stand_age_gami,dd[c],s=14,alpha=.35,color=ftc[ft],label=None)
        if len(dd)>10:
            xx=np.linspace(dd.stand_age_gami.min(),dd.stand_age_gami.max(),40)
            cf=np.polyfit(dd.stand_age_gami,dd[c],2); ax.plot(xx,np.polyval(cf,xx),'-',color=ftc[ft],lw=2.6,label=ft)
    # overall (confounded) quadratic
    xx=np.linspace(d.stand_age_gami.min(),d.stand_age_gami.max(),60); cf=np.polyfit(d.stand_age_gami,d[c],2)
    ax.plot(xx,np.polyval(cf,xx),'k--',lw=2,alpha=.7,label='전체(교란)')
    ax.set_xlabel("stand age (GAMI, yr)"); ax.set_ylabel(c); ax.set_title(l); ax.legend(fontsize=9); ax.grid(alpha=.2)
fig.suptitle("N10. forest type별 임령~다양성 — 전체 혹형(검은 점선)이 type 내에선 대부분 사라짐(=생물지리 교란)",fontsize=13)
fig.tight_layout(); fig.savefig(f"{F}/N10_age_diversity_byforesttype.png",dpi=120,bbox_inches='tight'); plt.close()
print("\nsaved N10_age_diversity_byforesttype.png")
