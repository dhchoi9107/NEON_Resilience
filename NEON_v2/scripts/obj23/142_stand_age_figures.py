"""
Stand-age figures:
  N02  GAMI vs Global-30m age cross-validation (30m saturates at ~40 = recency)
  N03  stand age (GAMI) vs structural complexity CHANGE RATE (per-year trend)
Bivariate (raw) scatter; domain-controlled mixed models are in 140_stand_age.py.
Outputs: figures/N02_age_crossval.png, N03_age_vs_structural_change.png, results/stand_age_crossval.csv
"""
import os, pandas as pd, numpy as np, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
def z(s): return (s-s.mean())/s.std()
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")

# ===== N02: cross-validation =====
g=pd.read_csv(f"{D}/plot_stand_age_gami.csv")[['plotID','siteID','stand_age_gami']]
m=pd.read_csv(f"{D}/plot_stand_age_30m.csv")[['plotID','stand_age_30m']]
d=g.merge(m,on='plotID').dropna()
young=d[d.stand_age_gami<40]
stats=dict(n=len(d),sat_pct=round(100*(d.stand_age_30m>=39).mean()),
           rho_all=round(st.spearmanr(d.stand_age_gami,d.stand_age_30m).correlation,3),
           rho_young=round(st.spearmanr(young.stand_age_gami,young.stand_age_30m).correlation,3),n_young=len(young))
pd.DataFrame([stats]).to_csv(f"{R}/stand_age_crossval.csv",index=False); print(stats)
fig,ax=plt.subplots(1,2,figsize=(12,5))
ax[0].scatter(d.stand_age_gami,d.stand_age_30m,s=14,alpha=.35,color='#00695c')
ax[0].axhline(40,color='r',ls='--',lw=1.2,label='30m 포화(~40, Landsat 레코드)')
ax[0].set_xlabel("GAMI 100m 임령 (yr, 진짜)"); ax[0].set_ylabel("Global 30m (yr, 검열≤40)")
ax[0].set_title(f"GAMI vs 30m — 노령림서 30m 포화\nSpearman 전체={stats['rho_all']}, 젊은(GAMI<40)={stats['rho_young']}"); ax[0].legend(fontsize=8)
gs=d.groupby('siteID').agg(gami=('stand_age_gami','median'),m30=('stand_age_30m','median'))
ax[1].scatter(gs.gami,gs.m30,s=45,color='#c62828')
for s,r in gs.iterrows(): ax[1].annotate(s,(r.gami,r.m30),fontsize=7,alpha=.7)
ax[1].axhline(40,color='r',ls='--',lw=1); ax[1].set_xlabel("사이트 중앙 GAMI 임령"); ax[1].set_ylabel("사이트 중앙 30m")
ax[1].set_title("사이트 수준 — 노령 사이트(WREF/TEAK/RMNP)만 30m 포화")
fig.suptitle("N02. Stand age 교차검증 — GAMI(진짜 임령) vs Global 30m(Landsat 검열)",fontsize=12)
fig.tight_layout(); fig.savefig(f"{F}/N02_age_crossval.png",dpi=120,bbox_inches='tight'); plt.close(); print("saved N02")

# ===== N03: stand age vs structural complexity change rate =====
df=pd.read_csv(f"{D}/FINAL_v2_pooled.csv").merge(
   pd.read_csv(f"{D}/plot_stand_age_gami.csv")[['plotID','stand_age_gami']],on='plotID',how='left')
df=df[(df.sample_coverage>=0.9)&df.stand_age_gami.notna()].copy(); A='stand_age_gami'
sites=sorted(df.siteID.unique())
cmap=plt.get_cmap('tab20'); scol={s:cmap(i%20) for i,s in enumerate(sites)}
METR=[('VCI_trend','VCI (수직복잡도)'),('FHD_trend','FHD (엽층다양성)'),('LAI_trend','LAI (엽면적)'),
      ('Rugosity_trend','Rugosity (표면거칠기)'),('Gini_trend','Gini (구조불평등)'),
      ('Deep_Gap_trend','Deep_Gap (갭비율)'),('Max_Ht_trend','Max_Ht (최고수고)'),('Canopy_Ht_trend','Canopy_Ht (평균수고)')]
fig,axes=plt.subplots(2,4,figsize=(20,9.5))
for ax,(col,lab) in zip(axes.ravel(),METR):
    dd=df[[A,col,'siteID']].dropna(); x=dd[A].values; y=dd[col].values
    ax.scatter(x,y,s=16,alpha=.6,c=[scol[s] for s in dd.siteID],edgecolors='none')
    sl,ic,r,p,se=st.linregress(x,y); sig=p<0.05
    xx=np.linspace(x.min(),x.max(),50); ax.plot(xx,ic+sl*xx,'-',color='k',lw=2.6)  # overall fit
    ax.axhline(0,color='k',ls=':',lw=1)
    star='***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
    ax.set_title(f"{lab}\nr={r:+.2f}{star}  (기울기 {sl:+.2e}/yr)",fontsize=10.5,color=('#4a148c' if sig else '#616161'))
    ax.set_xlabel("stand age (GAMI, yr)"); ax.set_ylabel(col); ax.grid(alpha=.25)
from matplotlib.lines import Line2D
handles=[Line2D([0],[0],marker='o',ls='',mfc=scol[s],mec='none',ms=8,label=s) for s in sites]
fig.legend(handles=handles,loc='lower center',ncol=10,fontsize=9,frameon=False,
           title="site (색)  ·  검은선=전체 회귀",bbox_to_anchor=(0.5,-0.02))
fig.suptitle("N03. 임분연령 vs 구조적 복잡도 변화속도(trend) — 점 색=site, 젊을수록 축적 빠름(VCI/LAI −), 노령=갭동태(Deep_Gap +)",fontsize=13)
fig.tight_layout(rect=[0,0.05,1,0.97]); fig.savefig(f"{F}/N03_age_vs_structural_change.png",dpi=120,bbox_inches='tight'); plt.close()
print("saved N03 (colored by site)")

# ===== N04: PARTIAL regression — pure age effect controlling for domain + site random =====
# Mixed model  z(trend) ~ z(age) + C(domain) + (1|site).  Partial residual for age:
#   pr = y - [C(domain) fit + site random effect]  (= beta_age*z_age + residual)
# Plotting pr vs age isolates the age effect net of domain & between-site differences.
fig,axes=plt.subplots(2,4,figsize=(20,9.5))
for ax,(col,lab) in zip(axes.ravel(),METR):
    d=df[[A,col,'siteID','domain']].dropna().copy()
    d['za']=z(d[A]); d['zt']=z(d[col])
    try:
        m=smf.mixedlm("zt ~ za + C(domain)",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
    except Exception:
        ax.set_title(f"{lab}\n(fit 실패)"); continue
    beta=m.fe_params['za']; p=m.pvalues['za']
    fe_pred=m.predict(d)                                   # X*beta (fixed: za + domain)
    try:
        remap={s:float(v.iloc[0]) for s,v in m.random_effects.items()}
        re=d['siteID'].map(remap).fillna(0.0).values
    except Exception:
        re=np.zeros(len(d))                                # singular site variance -> no RE
    pr=(d['zt'].values-(fe_pred.values+re))+beta*d['za'].values   # partial resid for age
    ax.scatter(d[A],pr,s=16,alpha=.6,c=[scol[s] for s in d.siteID],edgecolors='none')
    xx=np.linspace(d[A].min(),d[A].max(),50)
    ax.plot(xx, beta*(xx-d[A].mean())/d[A].std(),'-',color='k',lw=2.6)   # partial slope
    ax.axhline(0,color='k',ls=':',lw=1)
    star='***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
    ax.set_title(f"{lab}\nβ_age|(domain+site)={beta:+.3f}{star}  p={p:.2g}",fontsize=10.5,color=('#1a237e' if p<0.05 else '#616161'))
    ax.set_xlabel("stand age (GAMI, yr)"); ax.set_ylabel(f"{col} partial resid"); ax.grid(alpha=.25)
fig.legend(handles=handles,loc='lower center',ncol=10,fontsize=9,frameon=False,
           title="site (색)  ·  검은선=부분 임령효과(β_age)",bbox_to_anchor=(0.5,-0.02))
fig.suptitle("N04. 순수 임령효과 (부분회귀) — 도메인·site 랜덤효과 통제 후 임령 partial residual",fontsize=13)
fig.tight_layout(rect=[0,0.05,1,0.97]); fig.savefig(f"{F}/N04_age_partial_regression.png",dpi=120,bbox_inches='tight'); plt.close()
print("saved N04 (partial regression)")
