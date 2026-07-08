"""
Figures for the two gap-fills:
  K01 (Obj1): taxonomic diversity ~ RS-derived diversity (spectral + structural), forest plot.
  K02 (Obj2): diversity ~ land-use heterogeneity (scatter grid).
  K03 (Obj2): land-use heterogeneity MODERATES the RS<->taxonomic coupling (low vs high split).
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=120,bbox_inches='tight'); plt.close(fig); print(" saved",n)

# ===== K01: Obj1 forest — spectral vs structural diversity =====
o1=pd.read_csv(os.path.join(R,"obj1_specdiv.csv"))
resps=['Hill q1(alpha)','Hill q2(alpha)','LCBD turnover(beta)','LCBD nestedness(beta)']
catc={'spectral_div':'#2e7d32','structural_div':'#c62828'}
fig,axes=plt.subplots(1,4,figsize=(20,6),sharex=True)
for ax,rl in zip(axes,resps):
    d=o1[o1.response==rl].sort_values('beta')
    for i,(_,x) in enumerate(d.iterrows()):
        ax.scatter(x['beta'],i,c=catc[x['kind']],s=70,alpha=1 if x['p']<0.05 else .3)
    ax.set_yticks(range(len(d))); ax.set_yticklabels(d['predictor'],fontsize=9)
    ax.axvline(0,color='k',ls='--'); ax.set_title(rl,fontsize=11); ax.set_xlabel("beta (z)"); ax.grid(axis='x',alpha=.3)
axes[-1].legend(handles=[Patch(color=v,label={'spectral_div':'분광다양성','structural_div':'구조다양성'}[k]) for k,v in catc.items()],loc='lower right')
fig.suptitle("K01 (Obj1). 종다양성(알파 Hill + 베타 LCBD) ~ RS-유래 다양성 (분광 + 구조)",fontsize=13)
fig.tight_layout(); save(fig,"K01_obj1_rs_diversity.png")

# ===== K02: Obj2 diversity ~ land-use heterogeneity =====
m=pd.read_csv(os.path.join(D,"FINAL_v2_full.csv")); mc=m[m['sample_coverage']>=0.9]
pairs=[('lc_edge','파편화(edge density)'),('lc_forest_frac','산림 비율'),('lc_shannon','토지피복 Shannon')]
fig,axes=plt.subplots(2,3,figsize=(16,9))
for j,(resp,rl) in enumerate([('Hill_q1','Hill q1 (알파)'),('LCBD_turnover_rare','LCBD turnover (베타)')]):
    for k,(h,hl) in enumerate(pairs):
        ax=axes[j,k]; d=mc[[h,resp]].dropna()
        ax.scatter(d[h],d[resp],s=12,alpha=.3,color='#5d4037')
        sl,ic,r,p,se=st.linregress(d[h],d[resp]); xx=np.linspace(d[h].min(),d[h].max(),20)
        ax.plot(xx,ic+sl*xx,'r-',lw=2)
        star='***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
        ax.set_xlabel(hl); ax.set_ylabel(rl); ax.set_title(f"{rl} ~ {hl}  r={r:+.2f}{star}",fontsize=10)
fig.suptitle("K02 (Obj2). 종다양성 ~ 토지이용 이질성 (파편화↑→다양성↓, 산림비율↑→다양성↑)",fontsize=13)
fig.tight_layout(); save(fig,"K02_obj2_heterogeneity.png")

# ===== K03: heterogeneity moderates RS<->taxonomic coupling =====
fig,axes=plt.subplots(1,2,figsize=(14,5.5))
def z(s): return (s-s.mean())/s.std()
for ax,(rs,resp,hcol,tt) in zip(axes,[
     ('SAVI_mean','LCBD_turnover_rare','lc_shannon','SAVI-Turnover'),
     ('SAVI_mean','LCBD_turnover_rare','lc_edge','SAVI-Turnover')]):
    hl=hcol; d=mc[[rs,resp,hcol]].dropna(); med=d[hcol].median()
    for lohi,c,lab in [('low','#1565c0','저이질성'),('high','#d32f2f','고이질성')]:
        dd=d[d[hcol]<med] if lohi=='low' else d[d[hcol]>=med]
        ax.scatter(dd[rs],dd[resp],s=12,alpha=.3,color=c)
        sl,ic,r,p,se=st.linregress(dd[rs],dd[resp]); xx=np.linspace(dd[rs].min(),dd[rs].max(),20)
        ax.plot(xx,ic+sl*xx,'-',color=c,lw=2.4,label=f"{lab}: 기울기 {sl:+.2f}")
    ax.set_xlabel(rs); ax.set_ylabel(resp); ax.legend(fontsize=9); ax.set_title(f"{tt}  (분할: {hl})")
fig.suptitle("K03 (Obj2). 토지이용 이질성이 RS↔종 관계를 조절 — 이질성 高/低에서 SAVI–Turnover 기울기 다름",fontsize=12)
fig.tight_layout(); save(fig,"K03_obj2_moderation.png")
print("DONE gap figures")
