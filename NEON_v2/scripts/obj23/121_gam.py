"""
GAM (data-driven smooth) test of nonlinearity. For each key predictor: control domain
(residualize y on C(domain)), fit LinearGAM s(x) with lambda gridsearch, report EDF
(>1 = nonlinear), smooth-term p-value, and plot the fitted curve + 95% CI (shows hump/
threshold/saturation directly). Output: results/gam.csv + figure L02.
"""
import os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from pygam import LinearGAM, s
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
R=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\results"
F=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\figures"
def z(s_): return (s_-s_.mean())/s_.std()
m=pd.read_csv(D+r"\FINAL_v2_full.csv").merge(pd.read_csv(D+r"\plot_dhi_sentinel.csv"),on='plotID',how='left')
m=m.merge(pd.read_csv(D+r"\plot_disturbance_neon.csv")[['plotID','severity','disturbed']],on='plotID',how='left')
mc=m[m['sample_coverage']>=0.9].copy()

PAIRS=[('Hill_q1','dhi_cum_mean','생산성(DHI누적)',None),('Hill_q1','dhi_var_mean','생산성 계절변동',None),
       ('LCBD_turnover_rare','lc_edge','파편화(edge)',None),('Hill_q1','lc_forest_frac','산림비율',None),
       ('Hill_q1','VCI_mean','구조 VCI',None),('Hill_q1','specdiv_rao_q','분광다양성 RaoQ',None),
       ('Hill_q1','severity','교란강도',1),('LCBD_turnover_rare','severity','교란강도→Turnover',1)]
rows=[]; figd=[]
for resp,pred,lab,donly in PAIRS:
    d=(mc[mc.disturbed==1] if donly else mc)[[resp,pred,'domain']].dropna().copy()
    if len(d)<40: continue
    # domain 통제: y를 domain에 회귀 후 잔차
    d['y']=z(d[resp]);
    d['yr']=smf.ols("y ~ C(domain)",d).fit().resid
    X=z(d[pred]).values.reshape(-1,1); y=d['yr'].values
    gam=LinearGAM(s(0,n_splines=10)).gridsearch(X,y,progress=False)
    edf=gam.statistics_['edof']; pval=gam.statistics_['p_values'][0]
    XX=np.linspace(X.min(),X.max(),100).reshape(-1,1)
    pdep,ci=gam.partial_dependence(term=0,X=XX,width=.95)
    # 형태: 곡선 극점이 내부 & 봉우리면 hump
    imax=int(np.argmax(pdep)); imin=int(np.argmin(pdep))
    hump = 5<imax<95 and pdep[imax]>pdep[0] and pdep[imax]>pdep[-1]
    dip  = 5<imin<95 and pdep[imin]<pdep[0] and pdep[imin]<pdep[-1]
    shape=('단봉hump' if hump else 'U자dip' if dip else ('포화/단조'))
    rows.append(dict(resp=resp,pred=pred,n=len(d),EDF=round(edf,2),p=round(pval,4),shape=shape))
    figd.append((lab,XX.ravel(),pdep,ci,d['pred_z'] if False else X.ravel(),y,edf,pval,shape))
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"gam.csv"),index=False)
print(f"{'관계':32s}| n  | EDF  | p      | 형태")
for _,x in res.iterrows():
    print(f"{x['resp']+'~'+x['pred']:32s}| {x['n']:3d}| {x['EDF']:4.2f} | {x['p']:.4f} | {x['shape']}  {'(EDF>1.5=비선형)' if x['EDF']>1.5 else ''}")

fig,axes=plt.subplots(2,4,figsize=(20,9))
for ax,(lab,xx,pdep,ci,X,y,edf,pval,shape) in zip(axes.ravel(),figd):
    ax.scatter(X,y,s=10,alpha=.2,color='#90a4ae')
    ax.plot(xx,pdep,'r-',lw=2.5,label=f'GAM smooth (EDF={edf:.1f})')
    ax.fill_between(xx,ci[:,0],ci[:,1],color='r',alpha=.15)
    ax.set_title(f"{lab}\n({shape}, EDF={edf:.2f}, p={pval:.3f})",
                 color='#b71c1c' if 'hump' in shape or 'dip' in shape else '#455a64')
    ax.set_xlabel("z(predictor)"); ax.set_ylabel("domain-보정 다양성"); ax.grid(alpha=.3); ax.legend(fontsize=8)
for ax in axes.ravel()[len(figd):]: ax.axis('off')
fig.suptitle("L02. GAM 스무딩 — 데이터가 그린 곡선(EDF>1.5=비선형). 생산성·파편화 단봉, 나머지 단조",fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(F,"L02_gam.png"),dpi=120,bbox_inches='tight'); plt.close()
print("\nsaved L02_gam.png")
