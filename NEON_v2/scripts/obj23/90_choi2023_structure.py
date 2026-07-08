"""
Apply Choi et al. 2023 (J. Ecology) framework to our NEON data.
Choi 2023: multi-temporal LiDAR shows (a) PULSE (fire/defoliation) disturbances cause clearer
canopy structural change than PRESS (insect/disease, masked by compensatory growth); (b) INITIAL
canopy complexity (rugosity) BUFFERS the disturbance impact -> resilience.

Here: impact plots' pre/post (split at NEON disturbance year) ΔLiDAR structure + initial complexity.
Tests: (1) ΔStructure press vs pulse; (2) does initial rugosity buffer ΔCanopyHt & ΔRichness?
Output: results/choi2023.csv + figures J01-J03
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")
def z(s): return (s-s.mean())/s.std()
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=118,bbox_inches='tight'); plt.close(fig); print(" saved",n)

neon=pd.read_csv(os.path.join(D,"plot_disturbance_neon.csv"))
imp=neon[neon.disturbed==1][['plotID','dist_year','neon_dist_type','severity']].dropna(subset=['dist_year'])
# press vs pulse (Choi): press=insect/disease/slow mortality ; pulse=fire/wind/flood/harvest (acute)
PRESS={'insect','mortality','natural'}; PULSE={'fire','wind','flood','harvest'}
imp['regime']=imp['neon_dist_type'].map(lambda t:'press' if t in PRESS else ('pulse' if t in PULSE else 'other'))
splity=imp.set_index('plotID')['dist_year'].to_dict()

py=pd.read_csv(os.path.join(D,"per_year_v2.csv"))
LID=['Canopy_Ht','Rumple','Rugosity','Deep_Gap','LAI','Gini','FHD','Max_Ht']
rows=[]
for pid,g in py.groupby('plotID'):
    if pid not in splity: continue
    sy=splity[pid]; b=g[g.year<sy]; a=g[g.year>=sy]
    if len(b)<1 or len(a)<1: continue
    rec={'plotID':pid}
    for m in LID:
        if m in g: rec['d_'+m]=a[m].mean()-b[m].mean(); rec['pre_'+m]=b[m].mean()
    rows.append(rec)
ch=pd.DataFrame(rows).merge(imp,on='plotID',how='left')
print(f"impact plots with pre/post LiDAR: {len(ch)} | press {int((ch.regime=='press').sum())} pulse {int((ch.regime=='pulse').sum())}")

# initial complexity class (Choi: rugosity SD; high>+1, low<-1)
ch['init_complex']=ch['pre_Rugosity']; ch['zc']=z(ch['init_complex'])
ch['complex_class']=pd.cut(ch['zc'],[-99,-1,1,99],labels=['low','med','high'])

# ===== (1) ΔStructure press vs pulse =====
print("\n=== (1) ΔLiDAR structure: press vs pulse (Choi: pulse clearer) ===")
res=[]
for m in ['Canopy_Ht','Deep_Gap','Rugosity','LAI','Gini','FHD']:
    c=f'd_{m}'; pr=ch[ch.regime=='press'][c].dropna(); pu=ch[ch.regime=='pulse'][c].dropna()
    if len(pr)<5 or len(pu)<5: continue
    t,p=st.mannwhitneyu(pr,pu)
    print(f"  Δ{m:10s}: press median {pr.median():+.2f} | pulse median {pu.median():+.2f} | diff p={p:.3f}")
    res.append(dict(test='press_vs_pulse',metric=m,press_med=pr.median(),pulse_med=pu.median(),p=p))

# ===== (2) initial complexity buffers impact? =====
print("\n=== (2) 초기 복잡도(rugosity)가 교란 영향 완충? (Choi 핵심) ===")
ch['siteID']=ch.plotID.str[:4]
# diversity change (richness) merge
bac=pd.read_csv(os.path.join(D,"plot_baci_diversity_neon.csv"))[['plotID','drich','dq1']]
ch=ch.merge(bac,on='plotID',how='left')
for dv,lab in [('d_Canopy_Ht','ΔCanopyHt'),('d_FHD','ΔFHD'),('drich','ΔRichness'),('dq1','ΔHillq1')]:
    d=ch.dropna(subset=[dv,'zc'])
    if len(d)<20: continue
    r,p=st.pearsonr(d['zc'],d[dv])
    # also controlling severity
    try:
        mm=smf.mixedlm(f"{dv} ~ zc + severity",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
        b,pp=mm.fe_params['zc'],mm.pvalues['zc']
    except Exception: b,pp=np.nan,np.nan
    print(f"  {lab:11s} ~ init_complexity: r={r:+.2f} (p={p:.3f}) | lmm(+severity) beta={b:+.2f} p={pp:.3f}")
    res.append(dict(test='complexity_buffer',metric=lab,r=r,p=p,lmm_beta=b,lmm_p=pp))
pd.DataFrame(res).to_csv(os.path.join(R,"choi2023.csv"),index=False)

# ---- J01: ΔStructure press vs pulse ----
mets=['Canopy_Ht','Deep_Gap','Rugosity','LAI','Gini','FHD']
fig,axes=plt.subplots(2,3,figsize=(15,8))
for ax,m in zip(axes.ravel(),mets):
    c=f'd_{m}'; data=[ch[ch.regime=='press'][c].dropna(),ch[ch.regime=='pulse'][c].dropna()]
    ax.boxplot(data,labels=['press\n(insect/disease)','pulse\n(fire/wind)'],showfliers=False)
    ax.axhline(0,color='r',ls='--');ax.set_ylabel(f"Δ{m}");ax.set_title(f"Δ{m}")
fig.suptitle("J01. ΔLiDAR canopy structure (post−pre): press vs pulse disturbances [Choi 2023]",fontsize=13)
fig.tight_layout(); save(fig,"J01_press_pulse.png")

# ---- J02: initial complexity buffer ----
fig,axes=plt.subplots(1,3,figsize=(16,5))
for ax,(dv,lab) in zip(axes,[('d_Canopy_Ht','ΔCanopy height'),('d_FHD','ΔFHD'),('drich','ΔRichness')]):
    d=ch.dropna(subset=[dv,'init_complex'])
    ax.scatter(d['init_complex'],d[dv],s=22,alpha=.6,c=(d.regime=='pulse').map({True:'#d32f2f',False:'#1565c0'}))
    if len(d)>10:
        sl,ic,r,p,se=st.linregress(d['init_complex'],d[dv]);xx=np.linspace(d.init_complex.min(),d.init_complex.max(),20)
        ax.plot(xx,ic+sl*xx,'k-',lw=2,label=f"r={r:+.2f} p={p:.2g}")
    ax.axhline(0,color='gray',ls=':');ax.set_xlabel("initial canopy complexity (pre rugosity)")
    ax.set_ylabel(lab);ax.legend();ax.set_title(f"{lab} vs initial complexity")
fig.suptitle("J02. Does initial canopy complexity BUFFER disturbance impact? (red=pulse, blue=press) [Choi 2023]",fontsize=12)
fig.tight_layout(); save(fig,"J02_complexity_buffer.png")

# ---- J03: ΔStructure by complexity class x regime ----
fig,ax=plt.subplots(figsize=(9,5))
ch.boxplot(column='d_Canopy_Ht',by='complex_class',ax=ax)
ax.axhline(0,color='r',ls='--');ax.set_xlabel("initial complexity class");ax.set_ylabel("ΔCanopy height")
ax.set_title("J03. Canopy height change by initial complexity class");plt.suptitle("")
fig.tight_layout(); save(fig,"J03_complexity_class.png")
print("DONE Choi2023")
