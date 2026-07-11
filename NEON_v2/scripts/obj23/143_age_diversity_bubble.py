"""
N07: stand age (x) vs structural complexity change rate (y = VCI_trend, representative),
point SIZE = each of the 4 taxonomic diversity indices (one per panel), color = site.
Shows how species diversity is distributed over the age x complexity-change space.
Output: figures/N07_age_change_diversity_bubbles.png
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures")

df=pd.read_csv(f"{D}/FINAL_v2_pooled.csv").merge(
   pd.read_csv(f"{D}/plot_stand_age_gami.csv")[['plotID','stand_age_gami']],on='plotID',how='left')
df=df[(df.sample_coverage>=0.9)&df.stand_age_gami.notna()].copy()
A='stand_age_gami'; Y='VCI_trend'; Ylab="VCI_trend (수직복잡도 변화속도)"
DIV=[('Hill_q1','Hill q1 (알파)'),('Hill_q2','Hill q2 (알파)'),
     ('LCBD_turnover_rare','LCBD turnover (베타)'),('LCBD_nestedness_rare','LCBD nestedness (베타)')]
sites=sorted(df.siteID.unique()); cmap=plt.get_cmap('tab20'); scol={s:cmap(i%20) for i,s in enumerate(sites)}

def sizes(v):  # robust min-max (5-95 pct) -> marker area [20,320]
    lo,hi=np.nanpercentile(v,5),np.nanpercentile(v,95)
    s=np.clip((v-lo)/(hi-lo+1e-12),0,1); return 20+300*s

fig,axes=plt.subplots(2,2,figsize=(16,13))
for ax,(dv,dl) in zip(axes.ravel(),DIV):
    d=df[[A,Y,dv,'siteID']].dropna()
    sz=sizes(d[dv].values)
    ax.scatter(d[A],d[Y],s=sz,alpha=.55,c=[scol[s] for s in d.siteID],edgecolors='k',linewidths=.3)
    sl,ic,r,p,se=st.linregress(d[A],d[Y]); xx=np.linspace(d[A].min(),d[A].max(),50)
    ax.plot(xx,ic+sl*xx,'k--',lw=2,alpha=.8)
    ax.axhline(0,color='gray',ls=':',lw=1)
    ax.set_title(f"점 크기 = {dl}",fontsize=12)
    ax.set_xlabel("stand age (GAMI, yr)"); ax.set_ylabel(Ylab); ax.grid(alpha=.2)
    # per-panel size legend (small/mid/large of this index)
    lo,hi=np.nanpercentile(d[dv],5),np.nanpercentile(d[dv],95)
    qs=[lo,(lo+hi)/2,hi]
    hs=[Line2D([0],[0],marker='o',ls='',mfc='gray',mec='k',ms=np.sqrt(s0)/1.3,
        label=f"{q:.3g}") for q,s0 in zip(qs,sizes(np.array(qs)))]
    ax.legend(handles=hs,title=dl,loc='upper right',fontsize=8,labelspacing=1.1,borderpad=.8)
# shared site legend
sh=[Line2D([0],[0],marker='o',ls='',mfc=scol[s],mec='none',ms=8,label=s) for s in sites]
fig.legend(handles=sh,loc='lower center',ncol=10,fontsize=8.5,frameon=False,title="site (색)",bbox_to_anchor=(0.5,-0.02))
fig.suptitle("N07. 임령 × 구조 변화속도(VCI_trend), 점 크기 = 종다양성 4지수 — 다양성의 공간 분포",fontsize=14)
fig.tight_layout(rect=[0,0.04,1,0.97]); fig.savefig(f"{F}/N07_age_change_diversity_bubbles.png",dpi=120,bbox_inches='tight'); plt.close()
print("saved N07_age_change_diversity_bubbles.png")
