"""
(B) Stand age (succession) analyses. Age from GAMI 100m (2020) and/or Global 30m.
  (i)   diversity ~ age  (linear + quadratic -> hump?)
  (ii)  age MODERATES the RS<->diversity coupling:  div ~ z(RS)*z(age)
  (iii) does structural COMPLEXITY GROWTH (trend) depend on age? trend ~ z(age)
        (negative => younger stands complexify faster = growth hypothesis)
Outputs: results/stand_age_{models,moderation,trendage}.csv, figures/N01_stand_age.png
Re-runs automatically with whatever plot_stand_age_*.csv are present (gami / 30m).
"""
import os, numpy as np, pandas as pd, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"results"),os.path.join(ROOT,"figures")
def z(s): return (s-s.mean())/s.std()

df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
AGES={}
for src,col in [("gami","stand_age_gami"),("30m","stand_age_30m")]:
    fp=os.path.join(D,f"plot_stand_age_{src}.csv")
    if os.path.exists(fp):
        a=pd.read_csv(fp)
        if col in a: df=df.merge(a[['plotID',col]],on='plotID',how='left'); AGES[src]=col
print("age sources:",AGES)
df=df[df['sample_coverage']>=0.9].copy()

RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'Turnover','LCBD_nestedness_rare':'Nestedness'}
RS=['SAVI_mean','EVI_mean','NDVI_mean','Deep_Gap_mean','VCI_mean','Gini_mean','LAI_mean','FHD_mean']
TREND=['Rugosity_trend','FHD_trend','VCI_trend','Gini_trend','Deep_Gap_trend','Max_Ht_trend','LAI_trend','Canopy_Ht_trend']

def mm(d,formula,grp='siteID'):
    try: m=smf.mixedlm(formula,d,groups=d[grp]).fit(reml=True,method='lbfgs')
    except Exception: return None
    return m if m.converged else None

for src,acol in AGES.items():
    print(f"\n########## AGE SOURCE = {src} ({acol}) ##########")
    d0=df[df[acol].notna()].copy(); d0['zage']=z(d0[acol])
    print(f"plots with age & cov>=0.9: {len(d0)}  age med={d0[acol].median():.0f} [{d0[acol].min():.0f}-{d0[acol].max():.0f}]")

    # (i) diversity ~ age linear + quadratic
    rows=[]
    for r,rl in RESP.items():
        d=d0[[r,acol,'siteID','domain']].dropna().copy()
        if len(d)<40: continue
        d['zy']=z(d[r]); d['za']=z(d[acol])
        m=mm(d,"zy ~ za + I(za**2) + C(domain)")
        if m is None: continue
        bl,pl=m.fe_params['za'],m.pvalues['za']
        bq,pq=m.fe_params['I(za ** 2)'],m.pvalues['I(za ** 2)']
        rows.append(dict(src=src,response=rl,beta_lin=bl,p_lin=pl,beta_quad=bq,p_quad=pq,n=len(d),
                         hump=(bq<0 and pq<0.05)))
    md=pd.DataFrame(rows); md.to_csv(os.path.join(R,f"stand_age_models_{src}.csv"),index=False)
    print("\n(i) diversity ~ age (linear + quad):")
    for _,x in md.iterrows():
        tag=' <== HUMP' if x.hump else (' (U)' if x.beta_quad>0 and x.p_quad<0.05 else '')
        print(f"  {x.response:11s} lin={x.beta_lin:+.3f}(p={x.p_lin:.3g}) quad={x.beta_quad:+.3f}(p={x.p_quad:.3g}){tag}")

    # (ii) age moderates RS<->diversity
    rows=[]
    for r,rl in RESP.items():
        for rs in RS:
            if rs not in d0: continue
            d=d0[[r,rs,acol,'siteID','domain']].dropna().copy()
            if len(d)<40: continue
            d['zy']=z(d[r]); d['zrs']=z(d[rs]); d['za']=z(d[acol])
            m=mm(d,"zy ~ zrs*za + C(domain)")
            if m is None: continue
            ik=[k for k in m.fe_params.index if ':' in k]
            if not ik: continue
            rows.append(dict(src=src,response=rl,rs=rs,beta_int=m.fe_params[ik[0]],p_int=m.pvalues[ik[0]],n=len(d)))
    mo=pd.DataFrame(rows); mo.to_csv(os.path.join(R,f"stand_age_moderation_{src}.csv"),index=False)
    sig=mo[mo.p_int<0.05].sort_values('p_int')
    print(f"\n(ii) age moderates RS<->diversity: {len(sig)}/{len(mo)} sig")
    for _,x in sig.iterrows():
        print(f"  {x.response:11s} ~ {x.rs:14s} x age  int={x.beta_int:+.3f} p={x.p_int:.3g}")

    # (iii) structural complexity GROWTH (trend) ~ age
    rows=[]
    for tr in TREND:
        if tr not in d0: continue
        d=d0[[tr,acol,'siteID','domain']].dropna().copy()
        if len(d)<40: continue
        d['zt']=z(d[tr]); d['za']=z(d[acol])
        m=mm(d,"zt ~ za + C(domain)")
        if m is None: continue
        rows.append(dict(src=src,trend=tr,beta_age=m.fe_params['za'],p=m.pvalues['za'],n=len(d)))
    mt=pd.DataFrame(rows); mt.to_csv(os.path.join(R,f"stand_age_trendage_{src}.csv"),index=False)
    sig=mt[mt.p<0.05].sort_values('beta_age')
    print(f"\n(iii) complexity growth (trend) ~ age: {len(sig)}/{len(mt)} sig  (neg = younger complexify faster)")
    for _,x in sig.iterrows():
        print(f"  {x.trend:16s} beta_age={x.beta_age:+.3f} p={x.p:.3g}")

# ===== figure (primary source = gami if present else first) =====
src=list(AGES)[0]; acol=AGES[src]
d0=df[df[acol].notna()].copy()
fig,axes=plt.subplots(1,3,figsize=(18,5.2))
# panel1: Hill q1 ~ age with quadratic fit
ax=axes[0]; d=d0[['Hill_q1',acol]].dropna()
ax.scatter(d[acol],d['Hill_q1'],s=12,alpha=.35,color='#2e7d32')
xx=np.linspace(d[acol].min(),d[acol].max(),100)
cf=np.polyfit(d[acol],d['Hill_q1'],2); ax.plot(xx,np.polyval(cf,xx),'r-',lw=2.3)
ax.set_xlabel(f"stand age ({src}, yr)"); ax.set_ylabel("Hill q1"); ax.set_title("(i) 알파 다양성 ~ 임령")
# panel2: site median age vs site mean Hill q1
ax=axes[1]; g=d0.groupby('siteID').agg(age=(acol,'median'),h=('Hill_q1','mean')).dropna()
ax.scatter(g['age'],g['h'],s=45,color='#1565c0')
for s,r in g.iterrows(): ax.annotate(s,(r['age'],r['h']),fontsize=7,alpha=.7)
ax.set_xlabel("사이트 중앙 임령"); ax.set_ylabel("사이트 평균 Hill q1"); ax.set_title("(사이트 수준)")
# panel3: complexity growth (FHD_trend) ~ age
ax=axes[2]
if 'FHD_trend' in d0:
    d=d0[['FHD_trend',acol]].dropna()
    ax.scatter(d[acol],d['FHD_trend'],s=12,alpha=.35,color='#6a1b9a')
    sl,ic=np.polyfit(d[acol],d['FHD_trend'],1); ax.plot(xx,ic+sl*xx,'r-',lw=2.3)
    ax.axhline(0,color='k',ls=':',lw=1)
    ax.set_xlabel(f"stand age ({src}, yr)"); ax.set_ylabel("FHD_trend (복잡도 생장률)"); ax.set_title("(iii) 복잡화 속도 ~ 임령")
fig.suptitle(f"N01. Stand age (succession) — Obj: 천이단계와 다양성·복잡화 (age={src})",fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(F,"N01_stand_age.png"),dpi=120,bbox_inches='tight'); plt.close()
print("\nsaved N01_stand_age.png ; DONE")
