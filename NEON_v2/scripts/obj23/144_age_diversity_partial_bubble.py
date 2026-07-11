"""
N08: PARTIAL (added-variable) version of N07. All three variables residualized on
domain(fixed)+site(random):  x = age resid,  y = VCI_trend resid,  and the
DIVERSITY is also residualized -> size = |resid| (magnitude), color = sign
(red = more diverse than site/domain predicts, blue = less). One panel per index.
Shows, net of site/domain, where in the pure age x complexity-change space plots
are more/less diverse than expected. Output: figures/N08_partial_diversity_bubbles.png
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures")

df=pd.read_csv(f"{D}/FINAL_v2_pooled.csv").merge(
   pd.read_csv(f"{D}/plot_stand_age_gami.csv")[['plotID','stand_age_gami']],on='plotID',how='left')
df=df[(df.sample_coverage>=0.9)&df.stand_age_gami.notna()].copy()
A='stand_age_gami'; Y='VCI_trend'
DIV=[('Hill_q1','Hill q1 (알파)'),('Hill_q2','Hill q2 (알파)'),
     ('LCBD_turnover_rare','LCBD turnover (베타)'),('LCBD_nestedness_rare','LCBD nestedness (베타)')]

def resid(d,col):
    """within-site residual: control siteID as fixed effect (absorbs domain, nested in site).
    Robust vs the singular domain-fixed + site-random combo (single-site domains D08/D10)."""
    dd=d.copy(); dd['yy']=dd[col].astype(float)
    return smf.ols("yy ~ C(siteID)",dd).fit().resid.values

def sizes(v):
    a=np.abs(v); hi=np.nanpercentile(a,95); s=np.clip(a/(hi+1e-12),0,1); return 15+300*s

fig,axes=plt.subplots(2,2,figsize=(16,13))
for ax,(dv,dl) in zip(axes.ravel(),DIV):
    d=df[[A,Y,dv,'siteID','domain']].dropna()
    xr=resid(d,A); yr=resid(d,Y); dr=resid(d,dv)          # all residualized on site (fixed)
    col=np.where(dr>0,'#c62828','#1565c0')                 # red=above expected, blue=below
    ax.scatter(xr,yr,s=sizes(dr),c=col,alpha=.5,edgecolors='k',linewidths=.3)
    sl,ic,r,p,se=st.linregress(xr,yr); xx=np.linspace(xr.min(),xr.max(),50)
    ax.plot(xx,ic+sl*xx,'k--',lw=2,alpha=.85)
    ax.axhline(0,color='gray',ls=':',lw=1); ax.axvline(0,color='gray',ls=':',lw=1)
    star='***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
    ax.set_title(f"점 크기 = |{dl} 잔차|,  색: 빨강=기대↑ 파랑=기대↓\n(부분회귀 VCI_trend~age: r={r:+.2f}{star})",fontsize=11)
    ax.set_xlabel("stand age 잔차 (site 통제(고정))"); ax.set_ylabel("VCI_trend 잔차")
    ax.grid(alpha=.2)
# legends
csz=[Line2D([0],[0],marker='o',ls='',mfc='gray',mec='k',ms=np.sqrt(s0)/1.3,label=l)
     for l,s0 in zip(['작음','중간','큼'],sizes(np.array([0.0,0.5,1.0]))*np.array([0.1,1,2]))]
ax=axes[0,0]
sign=[Line2D([0],[0],marker='o',ls='',mfc='#c62828',mec='k',ms=9,label='기대보다 다양(+)'),
      Line2D([0],[0],marker='o',ls='',mfc='#1565c0',mec='k',ms=9,label='기대보다 낮음(−)')]
ax.legend(handles=sign,loc='upper right',fontsize=8.5)
fig.suptitle("N08. 부분회귀(added-variable) — x·y·다양성 모두 site 잔차화(고정). 순수 임령×복잡화 공간의 다양성 편차",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(f"{F}/N08_partial_diversity_bubbles.png",dpi=120,bbox_inches='tight'); plt.close()
print("saved N08_partial_diversity_bubbles.png")
