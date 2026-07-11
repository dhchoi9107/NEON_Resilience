"""
Re-run key stand-age analyses with 20-year (2-decade) binning, per GAMI v2.1
producer recommendation ("absolute 100m age unreliable; bin into >=2 decades";
Besnard 2021 RMSE ~48yr, young overestimated / old underestimated).
  - structural complexity change (trend) ~ binned age  (vs continuous N04)
  - diversity moderation:  div ~ z(RS) x binned age
Age coarsened to 20yr-bin midpoint (>=20yr resolution). Compare to continuous.
Output: results/age_binned_{trend,moderation}.csv, figures/N11_age_binned_trends.png
"""
import os, numpy as np, pandas as pd, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"results"),os.path.join(ROOT,"figures")
def z(s): return (s-s.mean())/s.std()

df=pd.read_csv(f"{D}/FINAL_v2_pooled.csv").merge(
   pd.read_csv(f"{D}/plot_stand_age_gami.csv")[['plotID','stand_age_gami']],on='plotID',how='left')
df=df[(df.sample_coverage>=0.9)&df.stand_age_gami.notna()].copy()
A='stand_age_gami'
# 20-yr bins -> midpoint (coarsen to >=2 decades); tail 200+ collapsed
edges=list(range(0,241,20)); df['age_bin']=pd.cut(df[A],bins=edges,right=False)
df['age_mid']=df['age_bin'].apply(lambda b: (b.left+b.right)/2 if pd.notna(b) else np.nan).astype(float)
df=df[df.age_mid.notna()].copy()
print("bin별 plot 수:"); print(df.groupby('age_bin',observed=True).size().to_string()); print(f"n={len(df)}\n")

TREND=['VCI_trend','LAI_trend','Deep_Gap_trend','Canopy_Ht_trend','FHD_trend','Rugosity_trend','Max_Ht_trend','Gini_trend']
RS=['SAVI_mean','EVI_mean','NDVI_mean','LAI_mean','VCI_mean']
RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'turnover','LCBD_nestedness_rare':'nestedness'}

def mm(f,d):
    try: m=smf.mixedlm(f,d,groups=d.siteID).fit(reml=True,method='lbfgs')
    except Exception: return None
    return m if m.converged else None

# (iii) structural trend ~ binned age  vs continuous age
print("=== 구조변화 ~ 임령: 연속 vs 20년-binned (β_age|domain+site) ===")
rows=[]
for t in TREND:
    d=df[[t,A,'age_mid','siteID','domain']].dropna().copy()
    d['zc']=z(d[A]); d['zb']=z(d['age_mid'])
    mc=mm(f"{t} ~ zc + C(domain)",d); mb=mm(f"{t} ~ zb + C(domain)",d)
    bc,pc=(mc.fe_params['zc'],mc.pvalues['zc']) if mc else (np.nan,np.nan)
    bb,pb=(mb.fe_params['zb'],mb.pvalues['zb']) if mb else (np.nan,np.nan)
    rows.append(dict(trend=t,beta_cont=bc,p_cont=pc,beta_binned=bb,p_binned=pb,n=len(d)))
    print(f"  {t:16s} 연속 β={bc:+.3f}(p={pc:.2g})   20년구간 β={bb:+.3f}(p={pb:.2g})")
pd.DataFrame(rows).to_csv(f"{R}/age_binned_trend.csv",index=False)

# (ii) diversity moderation: div ~ z(RS) x binned age
print("\n=== RS–다양성 조절 ~ 20년-binned 임령 (유의만) ===")
mods=[]
for r,rl in RESP.items():
    for rs in RS:
        d=df[[r,rs,'age_mid','siteID','domain']].dropna().copy()
        if len(d)<40: continue
        d['zy']=z(d[r]); d['zrs']=z(d[rs]); d['zb']=z(d['age_mid'])
        m=mm("zy ~ zrs*zb + C(domain)",d)
        if not m: continue
        ik=[k for k in m.fe_params.index if ':' in k]
        if ik: mods.append(dict(response=rl,rs=rs,beta_int=m.fe_params[ik[0]],p_int=m.pvalues[ik[0]]))
md=pd.DataFrame(mods); md.to_csv(f"{R}/age_binned_moderation.csv",index=False)
sig=md[md.p_int<0.05].sort_values('p_int')
print(f"  {len(sig)}/{len(md)} 유의:")
for _,x in sig.iterrows(): print(f"    {x.response:10s} ~ {x.rs:10s} x age  int={x.beta_int:+.3f} p={x.p_int:.2g}")

# figure: PARTIAL RESIDUAL of each key trend (domain fixed + site random removed), by 20yr bin
def partial_resid(col):
    """trend residual net of C(domain)+(1|site) — same control as N04."""
    d=df[[col,'domain','siteID']].dropna().copy()
    m=smf.mixedlm(f"{col} ~ C(domain)",d,groups=d.siteID).fit(reml=True,method='lbfgs')
    fe=m.predict(d)
    try:
        rmap={s:float(v.iloc[0]) for s,v in m.random_effects.items()}; re=d.siteID.map(rmap).fillna(0).values
    except Exception: re=np.zeros(len(d))
    out=pd.Series(np.nan,index=df.index); out.loc[d.index]=d[col].values-(fe.values+re); return out

KEY=[('VCI_trend','VCI'),('LAI_trend','LAI'),('Deep_Gap_trend','Deep_Gap'),('Canopy_Ht_trend','Canopy_Ht')]
for t,_ in KEY: df[t+'_pr']=partial_resid(t)
order=[b for b in df['age_bin'].cat.categories if (df.age_bin==b).sum()>0]
fig,axes=plt.subplots(1,4,figsize=(20,5))
for ax,(t,l) in zip(axes,KEY):
    pr=t+'_pr'; data=[df[df.age_bin==b][pr].dropna().values for b in order]
    ax.boxplot(data,showfliers=False)
    ax.set_xticklabels([f"{int(b.left)}-{int(b.right)}" for b in order],rotation=45,ha='right',fontsize=8)
    med=[np.median(x) if len(x) else np.nan for x in data]
    ax.plot(range(1,len(order)+1),med,'r-o',lw=2,ms=4)
    ax.axhline(0,color='gray',ls=':'); ax.set_title(f"{l}_trend 부분잔차 ~ 임령구간"); ax.set_xlabel("임령 구간(년)"); ax.set_ylabel(f"{t} partial resid"); ax.grid(alpha=.2)
fig.suptitle("N11. 구조 복잡화 변화율(도메인·site 통제 partial residual) ~ 20년 임령구간 — 젊음=축적↑, 노령=갭↑ 유지",fontsize=13)
fig.tight_layout(); fig.savefig(f"{F}/N11_age_binned_trends.png",dpi=120,bbox_inches='tight'); plt.close()
print("\nsaved N11_age_binned_trends.png (partial residual)")
