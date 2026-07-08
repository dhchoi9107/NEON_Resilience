"""
Nonlinearity test — are RS/productivity/disturbance <-> diversity relationships UNIMODAL (hump),
not linear? Linear/mean models miss classic ecological shapes (productivity-diversity hump,
intermediate disturbance). Test quadratic term (OLS + domain FE, cluster-robust SE by site) and
visualize decile-binned means. Output: results/nonlinear.csv + figure L01.
"""
import os, numpy as np, pandas as pd, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
R=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\results"
F=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\figures"
def z(s): return (s-s.mean())/s.std()
m=pd.read_csv(D+r"\FINAL_v2_full.csv").merge(pd.read_csv(D+r"\plot_dhi_sentinel.csv"),on='plotID',how='left')
m=m.merge(pd.read_csv(D+r"\plot_disturbance_neon.csv")[['plotID','severity','disturbed']],on='plotID',how='left')
mc=m[m['sample_coverage']>=0.9].copy()

def fit(resp,pred,d):
    d=d[[resp,pred,'siteID','domain']].dropna().copy()
    if len(d)<40: return None
    d['x']=z(d[pred]); d['y']=z(d[resp]); d['x2']=d['x']**2
    try:
        mm=smf.ols("y ~ x + x2 + C(domain)",d).fit(cov_type='cluster',cov_kwds={'groups':d['siteID']})
        b2,p2=mm.params['x2'],mm.pvalues['x2']
    except Exception: return None
    # vertex (극점) 위치: -b1/(2b2), 데이터범위 내면 진짜 단봉/U
    b1=mm.params['x']; vert=-b1/(2*b2) if b2!=0 else np.nan
    inside=(-2<vert<2)
    shape=('단봉(hump)' if (p2<0.05 and b2<0 and inside) else 'U자' if (p2<0.05 and b2>0 and inside) else '선형/단조')
    return dict(resp=resp,pred=pred,n=len(d),x2_beta=b2,x2_p=p2,vertex=vert,shape=shape,d=d)

PAIRS=[('Hill_q1','dhi_cum_mean','생산성(DHI누적)'),('Hill_q1','VCI_mean','구조 VCI'),
       ('LCBD_turnover_rare','lc_edge','파편화'),('Hill_q1','lc_forest_frac','산림비율'),
       ('Hill_q1','specdiv_rao_q','분광다양성 RaoQ'),('Hill_q1','FHD_mean','구조 FHD')]
rows=[]; figdata=[]
for resp,pred,lab in PAIRS:
    r=fit(resp,pred,mc)
    if r: rows.append({k:r[k] for k in ['resp','pred','n','x2_beta','x2_p','vertex','shape']}); figdata.append((r,lab))
# 교란 강도 (중간교란가설) — 교란 plot, 단순 모델
dz=mc[mc.disturbed==1]
for resp,lab in [('Hill_q1','교란강도→Hillq1'),('LCBD_turnover_rare','교란강도→Turnover')]:
    r=fit(resp,'severity',dz)
    if r: rows.append({k:r[k] for k in ['resp','pred','n','x2_beta','x2_p','vertex','shape']}); figdata.append((r,lab))
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"nonlinear.csv"),index=False)
print(f"{'관계':34s}| n  | 2차계수(p)       | 극점 | 형태")
for _,x in res.iterrows():
    print(f"{x['resp']+'~'+x['pred']:34s}| {x['n']:3d}| {x['x2_beta']:+.3f}(p={x['x2_p']:.3f})| {x['vertex']:+.1f} | {x['shape']}")

# ---- figure: decile-binned means + quadratic fit ----
n=len(figdata); fig,axes=plt.subplots(2,4,figsize=(20,9))
for ax,(r,lab) in zip(axes.ravel(),figdata):
    d=r['d']; d=d.copy(); d['q']=pd.qcut(d['x'],10,labels=False,duplicates='drop')
    g=d.groupby('q')['y'].agg(['mean','sem']); xq=d.groupby('q')['x'].mean()
    ax.errorbar(xq,g['mean'],yerr=g['sem'],fmt='o',color='#37474f',ms=6,capsize=3,label='10분위 평균±SE')
    xx=np.linspace(d['x'].min(),d['x'].max(),50)
    mm=smf.ols("y ~ x + x2",d.assign(x2=d['x']**2)).fit()
    ax.plot(xx,mm.params['Intercept']+mm.params['x']*xx+mm.params['x2']*xx**2,'r-',lw=2.3,label='2차 적합')
    sig='*' if r['x2_p']<0.05 else ''
    ax.set_title(f"{lab}\n({r['shape']}, 2차 p={r['x2_p']:.3f}{sig})",
                 color='#b71c1c' if '단봉' in r['shape'] or 'U자' in r['shape'] else '#455a64')
    ax.set_xlabel(f"z({r['pred']})"); ax.set_ylabel(f"z({r['resp']})"); ax.grid(alpha=.3); ax.legend(fontsize=7)
for ax in axes.ravel()[len(figdata):]: ax.axis('off')
fig.suptitle("L01. 비선형(단봉) 검정 — 선형이 놓친 hump: 생산성·구조·파편화-다양성 관계",fontsize=14)
fig.tight_layout(); fig.savefig(os.path.join(F,"L01_nonlinear.png"),dpi=120,bbox_inches='tight'); plt.close()
print("\nsaved L01_nonlinear.png")
