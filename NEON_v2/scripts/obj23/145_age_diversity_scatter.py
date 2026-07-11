"""
N09: stand age (GAMI) vs each of the 4 taxonomic diversity indices (extends N01's
Hill-q1 panel to all four). Small points = plots (site-colored); large markers =
site means (the meaningful level, since age is ~75% between-site, ICC=0.75).
Annotates bivariate r and the mixed-model (domain fixed + site random) linear/quad p.
Output: figures/N09_age_vs_diversity.png
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures")
def z(s): return (s-s.mean())/s.std()

df=pd.read_csv(f"{D}/FINAL_v2_pooled.csv").merge(
   pd.read_csv(f"{D}/plot_stand_age_gami.csv")[['plotID','stand_age_gami']],on='plotID',how='left')
df=df[(df.sample_coverage>=0.9)&df.stand_age_gami.notna()].copy()
A='stand_age_gami'
DIV=[('Hill_q1','Hill q1 (알파)'),('Hill_q2','Hill q2 (알파)'),
     ('LCBD_turnover_rare','LCBD turnover (베타)'),('LCBD_nestedness_rare','LCBD nestedness (베타)')]
sites=sorted(df.siteID.unique()); cmap=plt.get_cmap('tab20'); scol={s:cmap(i%20) for i,s in enumerate(sites)}

def mixed_lin_quad(d,resp):
    dd=d[[A,resp,'siteID','domain']].dropna().copy()
    dd['zy']=z(dd[resp]); dd['za']=z(dd[A])
    try:
        m=smf.mixedlm("zy ~ za + I(za**2) + C(domain)",dd,groups=dd['siteID']).fit(reml=True,method='lbfgs')
        return m.pvalues.get('za',np.nan),m.pvalues.get('I(za ** 2)',np.nan)
    except Exception: return np.nan,np.nan

fig,axes=plt.subplots(2,2,figsize=(15,12))
for ax,(dv,dl) in zip(axes.ravel(),DIV):
    d=df[[A,dv,'siteID']].dropna()
    ax.scatter(d[A].values,d[dv].values,s=13,alpha=.30,c=[scol[s] for s in d.siteID.values],edgecolors='none')
    g=d.groupby('siteID').agg(age=(A,'median'),div=(dv,'mean'),n=(dv,'size')).reset_index()
    ax.scatter(g['age'].values,g['div'].values,s=np.clip(g['n'].values*5,50,320),
               c=[scol[s] for s in g['siteID'].values],edgecolors='k',lw=.8,zorder=4)
    # plot-level quadratic (black dashed)
    xx=np.linspace(d[A].min(),d[A].max(),100); cf=np.polyfit(d[A],d[dv],2); ax.plot(xx,np.polyval(cf,xx),'k--',lw=1.6,alpha=.7)
    # between-site QUADRATIC through site means (red curve) = meaningful level
    cfb=np.polyfit(g['age'],g['div'],2); ax.plot(xx,np.polyval(cfb,xx),'-',color='#b71c1c',lw=2.6,alpha=.9)
    gg=g.rename(columns={'div':'yy'}).copy(); gg['za']=(gg['age']-gg['age'].mean())/gg['age'].std()
    ob=smf.ols("yy ~ za + I(za**2)",gg).fit(); sq=ob.pvalues['I(za ** 2)']; sqb=ob.params['I(za ** 2)']
    pl,pq=mixed_lin_quad(df,dv)
    shape='혹형∩' if sqb<0 else 'U자∪'
    ax.set_title(f"{dl}\n사이트평균 2차 {shape} p={sq:.2g} · 혼합모델(domain통제) quad p={pq:.2g}",fontsize=10.5)
    ax.set_xlabel("stand age (GAMI, yr)"); ax.set_ylabel(dv); ax.grid(alpha=.2)
sh=[Line2D([0],[0],marker='o',ls='',mfc=scol[s],mec='none',ms=8,label=s) for s in sites]
fig.legend(handles=sh,loc='lower center',ncol=10,fontsize=8.5,frameon=False,
           title="site (작은점=plot, 큰점=사이트평균) · 빨강곡선=사이트평균 2차 · 검은점선=plot 2차",bbox_to_anchor=(0.5,-0.02))
fig.suptitle("N09. 임령 vs 종다양성 4지수 — 사이트평균(빨강곡선) 혹형∩(nestedness만 U∪), 단 domain통제시 소멸=생물지리 교란",fontsize=13)
fig.tight_layout(rect=[0,0.04,1,0.97]); fig.savefig(f"{F}/N09_age_vs_diversity.png",dpi=120,bbox_inches='tight'); plt.close()
print("saved N09_age_vs_diversity.png")
