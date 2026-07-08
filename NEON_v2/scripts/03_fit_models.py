"""
NEON_v2 STEP 3 — Mixed models: POOLED diversity ~ z(RS feature) + C(domain) + (1|site).
coverage>=0.9. Spectral predictors = NEON .002 BRDF only (8 years).
Mirrors prior model spec exactly; only the VI source changed (single, consistent).
Output: results/v2_coeff.csv, figures/Fig_v2_forest.png
"""
import os, numpy as np, pandas as pd
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D, R, F = os.path.join(ROOT,"data"), os.path.join(ROOT,"results"), os.path.join(ROOT,"figures")
def z(s): return (s-s.mean())/s.std()

df = pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
df = df[df['sample_coverage']>=0.9].copy()
print("plots (cov>=0.9):",len(df),"| sites:",df['siteID'].nunique())

STRUCT=['Canopy_Ht','Max_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Vert_CV','Gini','VCI','FHD','LAI','Q95','Ht_Ratio']
VI=['NDVI','EVI','ARVI','SAVI']   # NEON .002 BRDF, fPAR/PRI/LAI_opt excluded
PRED=STRUCT+VI; FEAT=['mean','sd','trend']
RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'Turnover(rare)','LCBD_nestedness_rare':'Nestedness(rare)'}
for r in RESP:
    if r in df: df[f"z_{r}"]=z(df[r])

def run(resp,col):
    if col not in df: return None
    d=df[[resp,col,'siteID','domain']].dropna()
    if len(d)<40 or d['siteID'].nunique()<3 or d['domain'].nunique()<2: return None
    d=d.copy(); d['zx']=z(d[col])
    try: m=smf.mixedlm(f"{resp} ~ zx + C(domain)",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
    except Exception: return None
    if not m.converged: return None
    b=m.fe_params['zx']; X=np.asarray(m.model.exog); fe=np.asarray(m.fe_params); yh=X@fe
    try: vs=max(float(m.cov_re.iloc[0,0]),0)
    except Exception: vs=0
    tot=np.var(yh)+vs+m.scale
    return dict(beta=b,se=m.bse_fe['zx'],p=m.pvalues['zx'],n=len(d),r2=b**2*np.var(d['zx'])/tot if tot>0 else 0)

rows=[]
for resp,rl in RESP.items():
    for p in PRED:
        for feat in FEAT:
            o=run(f"z_{resp}",f"{p}_{feat}")
            if o: rows.append(dict(response=rl,predictor=p,feature=feat,
                                   category='Structural' if p in STRUCT else 'VI',**o))
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"v2_coeff.csv"),index=False)
print("models:",len(res))

print("\n=== top |beta| per response x feature (p<0.05) ===")
for rl in RESP.values():
    for feat in FEAT:
        d=res[(res['response']==rl)&(res['feature']==feat)&(res['p']<0.05)]
        if len(d):
            top=d.loc[d['beta'].abs().idxmax()]
            print(f"  {rl:16s} {feat:5s}: {top['predictor']}={top['beta']:+.2f} (R2={top['r2']:.3f}) [{len(d)} sig]")
        else: print(f"  {rl:16s} {feat:5s}: (none sig)")

# forest figure (Hill q1, consistent predictor order across feature panels)
fig,axes=plt.subplots(1,3,figsize=(18,8),sharey=True)
catc={'Structural':'#c62828','VI':'#2e7d32'}
base=res[(res['response']=='Hill q1')&(res['feature']=='mean')]
order=base.sort_values('beta')['predictor'].tolist()
catmap=res.drop_duplicates('predictor').set_index('predictor')['category'].to_dict()
for ax,feat in zip(axes,FEAT):
    d=res[(res['response']=='Hill q1')&(res['feature']==feat)].set_index('predictor')
    for i,pname in enumerate(order):
        if pname not in d.index: continue
        b=d.loc[pname,'beta']; se=d.loc[pname,'se']; p=d.loc[pname,'p']
        ax.errorbar(b,i,xerr=1.96*se,fmt='none',ecolor='gray',zorder=2)
        ax.scatter(b,i,c=catc[catmap[pname]],s=70,alpha=1 if p<0.05 else .25,zorder=3)
    ax.set_yticks(range(len(order)));ax.set_yticklabels(order,fontsize=8)
    ax.axvline(0,color='k',ls='--',lw=1);ax.set_title(f"Hill q1 ~ RS {feat}");ax.set_xlabel("beta");ax.grid(axis='x',alpha=.3)
axes[-1].legend(handles=[Patch(color=v,label=k) for k,v in catc.items()],loc='lower right')
fig.suptitle("NEON_v2 — Pooled Hill q1 ~ per-year RS features (VI = NEON .002 BRDF only)",fontsize=14)
fig.tight_layout(); fig.savefig(os.path.join(F,"Fig_v2_forest.png"),dpi=120,bbox_inches='tight'); plt.close()
print("\nForest figure saved. DONE")
