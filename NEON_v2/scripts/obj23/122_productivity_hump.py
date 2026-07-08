"""
Obj 3 headline figure — hump-shaped productivity-diversity relationship (Sentinel-2 DHI).
Clean panels: Hill q1 & q2 vs DHI cumulative productivity, GAM smooth + 95% CI + decile means.
Domain-controlled. Marks the productivity level of peak diversity.
Output: figure L03_productivity_hump.png + results row.
"""
import os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from pygam import LinearGAM, s
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"; F=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\figures"
def z(s_): return (s_-s_.mean())/s_.std()
m=pd.read_csv(D+r"\FINAL_v2_full.csv").merge(pd.read_csv(D+r"\plot_dhi_sentinel.csv"),on='plotID',how='left')
mc=m[m['sample_coverage']>=0.9].copy()

fig,axes=plt.subplots(1,2,figsize=(14,5.6))
for ax,(resp,rl) in zip(axes,[('Hill_q1','Hill q1 (알파 다양성)'),('Hill_q2','Hill q2 (알파 다양성)')]):
    d=mc[[resp,'dhi_cum_mean','domain']].dropna().copy()
    d['y']=z(d[resp]); d['yr']=smf.ols("y ~ C(domain)",d).fit().resid
    X=z(d['dhi_cum_mean']).values.reshape(-1,1); y=d['yr'].values
    gam=LinearGAM(s(0,n_splines=8)).gridsearch(X,y,progress=False)
    edf=gam.statistics_['edof']; p=gam.statistics_['p_values'][0]
    XX=np.linspace(X.min(),X.max(),100).reshape(-1,1); pdep,ci=gam.partial_dependence(term=0,X=XX,width=.95)
    ax.scatter(X.ravel(),y,s=12,alpha=.18,color='#78909c')
    # decile means
    d['q']=pd.qcut(z(d['dhi_cum_mean']),10,labels=False,duplicates='drop')
    g=d.groupby('q').apply(lambda t:pd.Series({'x':z(d['dhi_cum_mean']).loc[t.index].mean(),'m':t['yr'].mean(),'se':t['yr'].sem()}))
    ax.errorbar(g['x'],g['m'],yerr=g['se'],fmt='o',color='#263238',ms=6,capsize=3,zorder=5,label='10분위 평균±SE')
    ax.plot(XX.ravel(),pdep,'r-',lw=2.6,label=f'GAM smooth (EDF={edf:.1f}, p={p:.3f})')
    ax.fill_between(XX.ravel(),ci[:,0],ci[:,1],color='r',alpha=.15)
    peak=XX.ravel()[int(np.argmax(pdep))]
    ax.axvline(peak,color='green',ls='--',lw=1.5,label=f'다양성 최대 지점 (z={peak:+.1f})')
    ax.set_xlabel("Sentinel-2 DHI 누적생산성 (z)"); ax.set_ylabel(f"{rl} (도메인 보정)"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    ax.set_title(f"{rl} ~ 생산성 — 단봉(hump)")
fig.suptitle("L03 (Obj3). 생산성–다양성 혹형 — 중간 생산성에서 종다양성 최대 (Sentinel-2 DHI, GAM)",fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(F,"L03_productivity_hump.png"),dpi=125,bbox_inches='tight'); plt.close()
print("saved L03_productivity_hump.png")
