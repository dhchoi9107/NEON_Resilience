"""
Obj 2 (continuous) — disturbance SEVERITY (intensity) + RECENCY (time-since) analysis.
  (A) main effects : diversity ~ z(severity)          [all plots; intact=0]
  (B) recovery     : diversity ~ z(recency)+z(severity)  [disturbed only -> does diversity recover with time?]
  (C) moderation   : diversity ~ z(RS) * z(severity)   [does intensity bend the RS<->taxonomic slope?]
  (D) moderation   : diversity ~ z(RS) * z(recency)    [disturbed; does the bend fade with time?]
Outputs: results/obj2_severity_recency.csv + figures H01-H04
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")
def z(s): return (s-s.mean())/s.std()
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=118,bbox_inches='tight'); plt.close(fig); print(" saved",n)

base=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
dist=pd.read_csv(os.path.join(D,"plot_disturbance_robust.csv"))
dhi=pd.read_csv(os.path.join(D,"plot_dhi_sentinel.csv"))
m=base.merge(dist,on='plotID',how='left').merge(dhi,on='plotID',how='left')
m=m[m['sample_coverage']>=0.9].copy()
for c in ['severity','sev_spec','sev_loss']: m[c]=m[c].fillna(0)
m['disturbed']=m['disturbed'].fillna(0)
print(f"plots {len(m)} | disturbed {int(m.disturbed.sum())} | severity>0 {(m.severity>0).sum()}")

RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'Turnover','LCBD_nestedness_rare':'Nestedness'}
for r in RESP: m[f"z_{r}"]=z(m[r])
RS=['SAVI_mean','EVI_mean','FHD_mean','LAI_mean','Deep_Gap_mean','dhi_cum_mean']
dz=m[m.disturbed==1].copy()
print(f"disturbed subset {len(dz)} (for recency/recovery)")

def fit(d,resp,formula):
    d=d.dropna(subset=[resp]+[c for c in ['zx','zsev','zrec'] if c in formula])
    if len(d)<40 or d['siteID'].nunique()<3: return None
    try: mm=smf.mixedlm(f"{resp} ~ {formula}",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
    except Exception: return None
    return mm if mm.converged else None

rows=[]
for resp,rl in RESP.items():
    zr=f"z_{resp}"
    # (A) severity main effect (all plots)
    d=m.copy(); d['zsev']=z(d['severity'])
    mm=fit(d,zr,"zsev + C(domain)")
    if mm is not None: rows.append(dict(test='A_severity',response=rl,term='severity',
        beta=mm.fe_params['zsev'],p=mm.pvalues['zsev'],n=len(d)))
    # (B) recovery: recency (disturbed)
    d=dz.copy(); d['zrec']=z(d['recency']); d['zsev']=z(d['severity'])
    mm=fit(d,zr,"zrec + zsev + C(domain)")
    if mm is not None:
        rows.append(dict(test='B_recency',response=rl,term='recency',beta=mm.fe_params['zrec'],p=mm.pvalues['zrec'],n=len(d)))
        rows.append(dict(test='B_recency',response=rl,term='severity',beta=mm.fe_params['zsev'],p=mm.pvalues['zsev'],n=len(d)))
    # (C) moderation by severity ; (D) by recency
    for rs in RS:
        d=m.copy(); d['zx']=z(d[rs]); d['zsev']=z(d['severity'])
        mm=fit(d,zr,"zx * zsev + C(domain)")
        if mm is not None and 'zx:zsev' in mm.fe_params:
            rows.append(dict(test='C_sev_x_RS',response=rl,term=rs,beta=mm.fe_params['zx:zsev'],p=mm.pvalues['zx:zsev'],n=len(d)))
        d=dz.copy(); d['zx']=z(d[rs]); d['zrec']=z(d['recency'])
        mm=fit(d,zr,"zx * zrec + C(domain)")
        if mm is not None and 'zx:zrec' in mm.fe_params:
            rows.append(dict(test='D_rec_x_RS',response=rl,term=rs,beta=mm.fe_params['zx:zrec'],p=mm.pvalues['zx:zrec'],n=len(d)))
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"obj2_severity_recency.csv"),index=False)
print("\n=== significant (p<0.05) ===")
for _,x in res[res.p<0.05].sort_values('p').iterrows():
    print(f"  [{x['test']:11s}] {x['response']:10s} {x['term']:14s} beta={x['beta']:+.3f} p={x['p']:.3g}")

# ---- H01: severity gradient ----
fig,axes=plt.subplots(1,4,figsize=(20,4.6))
for ax,(resp,rl) in zip(axes,RESP.items()):
    d=m[m.severity>0][['severity',resp]].dropna()
    ax.scatter(m['severity'],m[resp],s=10,alpha=.25,color='#8d6e63')
    if len(d)>10:
        sl,ic,r,p,se=st.linregress(d['severity'],d[resp]);xx=np.linspace(0,m['severity'].max(),20)
        ax.plot(xx,ic+sl*xx,'r-',lw=2,label=f"disturbed r={r:+.2f} p={p:.2g}")
    ax.set_xlabel("disturbance severity");ax.set_ylabel(rl);ax.legend(fontsize=8);ax.set_title(rl)
fig.suptitle("H01. Diversity vs disturbance SEVERITY (intensity)",fontsize=14); fig.tight_layout(); save(fig,"H01_severity.png")

# ---- H02: recovery vs recency, colored by severity ----
fig,axes=plt.subplots(1,4,figsize=(20,4.6))
for ax,(resp,rl) in zip(axes,RESP.items()):
    d=dz[['recency','severity',resp]].dropna()
    sc=ax.scatter(d['recency'],d[resp],s=18,c=d['severity'],cmap='YlOrRd',alpha=.7)
    if len(d)>10:
        sl,ic,r,p,se=st.linregress(d['recency'],d[resp]);xx=np.linspace(0,d['recency'].max(),20)
        ax.plot(xx,ic+sl*xx,'k-',lw=2,label=f"r={r:+.2f} p={p:.2g}")
    ax.set_xlabel("years since disturbance");ax.set_ylabel(rl);ax.legend(fontsize=8);ax.set_title(rl)
fig.colorbar(sc,ax=axes[-1],label="severity")
fig.suptitle("H02. Recovery: diversity vs RECENCY (disturbed plots; color=severity)",fontsize=14); fig.tight_layout(); save(fig,"H02_recency_recovery.png")

# ---- H03: moderation coefficients (severity x RS, recency x RS) ----
fig,axes=plt.subplots(1,2,figsize=(15,6))
for ax,(test,tit) in zip(axes,[('C_sev_x_RS','RS x SEVERITY'),('D_rec_x_RS','RS x RECENCY')]):
    d=res[res.test==test]
    if d.empty: ax.axis('off'); continue
    labs=[f"{r['response']}~{r['term'].replace('_mean','')}" for _,r in d.iterrows()]
    y=range(len(d))
    ax.scatter(d['beta'],y,c=['#d32f2f' if p<0.05 else '#90a4ae' for p in d['p']],s=60)
    ax.axvline(0,color='k',ls='--');ax.set_yticks(list(y));ax.set_yticklabels(labs,fontsize=8)
    ax.set_xlabel("interaction beta");ax.set_title(tit)
fig.suptitle("H03. Does disturbance intensity/recency bend the RS<->taxonomic slope?",fontsize=13); fig.tight_layout(); save(fig,"H03_moderation.png")
print("DONE severity/recency")
