"""
BACI step 2 — models + figures. Resolves the time-invariance problem:
  (A) BACI test:  ΔDiversity ~ impact (+severity)  -> did disturbance change taxonomic diversity
                  beyond control (background) change?
  (B) severity:   ΔDiversity ~ severity  (impact plots) -> dose-response
  (C) change<->change coupling: ΔDiversity ~ ΔRS (ΔNBR / Δstructure)
      -> does RS CHANGE track taxonomic CHANGE? (both vary in time; no pooled-response problem)
Outputs: results/baci.csv + figures I01-I03
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")
def z(s): return (s-s.mean())/s.std()
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=118,bbox_inches='tight'); plt.close(fig); print(" saved",n)

bac=pd.read_csv(os.path.join(D,"plot_baci_diversity.csv"))
bac['siteID']=bac['plotID'].str[:4]
DOM={"HARV":"D01","BART":"D01","SCBI":"D02","SERC":"D02","BLAN":"D02","GRSM":"D07","MLBS":"D07","ORNL":"D07",
 "JERC":"D03","OSBS":"D03","TALL":"D08","UNDE":"D05","STEI":"D05","TREE":"D05","WREF":"D16","ABBY":"D16",
 "SOAP":"D17","TEAK":"D17","RMNP":"D10"}
bac['domain']=bac['siteID'].map(DOM)

# ---- split year (same logic as 70) to build ΔRS ----
rob=pd.read_csv(os.path.join(D,"plot_disturbance_robust.csv"))
imp=rob[rob.dist_year.notna()][['plotID','dist_year']]
impmap=imp.set_index('plotID')['dist_year'].to_dict()
site_med=imp.assign(s=imp.plotID.str[:4]).groupby('s')['dist_year'].median().to_dict()
gmed=float(imp['dist_year'].median())
def split_year(pid): return impmap.get(pid, site_med.get(pid[:4],gmed))

def delta(df,val,ycol='year'):
    out={}
    for pid,g in df.groupby('plotID'):
        sy=split_year(pid); b=g[g[ycol]<sy][val]; a=g[g[ycol]>=sy][val]
        if len(b) and len(a): out[pid]=a.mean()-b.mean()
    return pd.Series(out,name='d'+val)

# ΔNBR (Sentinel annual)
nbr=pd.read_csv(os.path.join(D,"plot_nbr_annual.csv"))
dnbr=delta(nbr,'nbr')
# Δstructure (LiDAR per-year)
py=pd.read_csv(os.path.join(D,"per_year_v2.csv"))
dFHD=delta(py,'FHD'); dDG=delta(py,'Deep_Gap'); dLAI=delta(py,'LAI')
ch=pd.concat([dnbr,dFHD,dDG,dLAI],axis=1).reset_index().rename(columns={'index':'plotID'})
m=bac.merge(ch,on='plotID',how='left')
print("plots:",len(m),"| with ΔNBR:",m['dnbr'].notna().sum())

# ===== (A) BACI test =====
print("\n=== (A) BACI: ΔDiversity ~ impact (+severity) ===")
rows=[]
for dv,lab in [('dq1','ΔHill q1'),('dq2','ΔHill q2'),('drich','ΔRichness')]:
    d=m.dropna(subset=[dv,'impact','siteID','domain']).copy()
    try:
        mm=smf.mixedlm(f"{dv} ~ impact + C(domain)",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
        b=mm.fe_params['impact']; p=mm.pvalues['impact']
        print(f"  {lab:11s}: impact beta={b:+.3f} p={p:.3g} (impact n={int(d.impact.sum())}, control {int((d.impact==0).sum())})")
        rows.append(dict(test='BACI',response=lab,term='impact',beta=b,p=p))
    except Exception as e: print(f"  {lab}: fail {str(e)[:40]}")

# ===== (B) severity dose-response (impact only) =====
print("\n=== (B) ΔDiversity ~ severity (impact plots) ===")
di=m[m.impact==1].copy()
for dv,lab in [('dq1','ΔHill q1'),('dq2','ΔHill q2')]:
    d=di.dropna(subset=[dv,'severity','siteID'])
    if len(d)<25: continue
    d=d.copy(); d['zs']=z(d['severity'])
    try:
        mm=smf.mixedlm(f"{dv} ~ zs",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
        print(f"  {lab:11s}: severity beta={mm.fe_params['zs']:+.3f} p={mm.pvalues['zs']:.3g} n={len(d)}")
        rows.append(dict(test='severity',response=lab,term='severity',beta=mm.fe_params['zs'],p=mm.pvalues['zs']))
    except Exception: pass

# ===== (C) change<->change coupling =====
print("\n=== (C) ΔDiversity ~ ΔRS (change tracks change) ===")
for dv,lab in [('dq1','ΔHill q1'),('dq2','ΔHill q2')]:
    for rs,rl in [('dnbr','ΔNBR'),('dFHD','ΔFHD'),('dDeep_Gap','ΔDeepGap'),('dLAI','ΔLAI')]:
        d=m.dropna(subset=[dv,rs])
        if len(d)<30: continue
        r,p=st.pearsonr(d[rs],d[dv])
        rows.append(dict(test='coupling',response=lab,term=rl,beta=r,p=p))
        if p<0.1: print(f"  {lab:11s} ~ {rl:9s}: r={r:+.3f} p={p:.3g} n={len(d)}")
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"baci.csv"),index=False)

# ---- I01: BACI boxplot ----
fig,axes=plt.subplots(1,3,figsize=(15,5))
for ax,(dv,lab) in zip(axes,[('dq1','ΔHill q1'),('dq2','ΔHill q2'),('drich','ΔRichness')]):
    d=m.dropna(subset=[dv])
    ax.boxplot([d[d.impact==0][dv],d[d.impact==1][dv]],labels=['control','impact'],showfliers=False)
    ax.axhline(0,color='r',ls='--');ax.set_ylabel(lab);ax.set_title(lab)
fig.suptitle("I01. BACI — taxonomic diversity CHANGE (after−before): impact vs control",fontsize=13)
fig.tight_layout(); save(fig,"I01_baci_boxplot.png")

# ---- I02: severity dose-response ----
fig,axes=plt.subplots(1,2,figsize=(12,5))
for ax,(dv,lab) in zip(axes,[('dq1','ΔHill q1'),('dq2','ΔHill q2')]):
    d=di.dropna(subset=[dv,'severity'])
    ax.scatter(d['severity'],d[dv],s=20,alpha=.6,color='#c62828')
    if len(d)>10:
        sl,ic,r,p,se=st.linregress(d['severity'],d[dv]);xx=np.linspace(d.severity.min(),d.severity.max(),20)
        ax.plot(xx,ic+sl*xx,'k-',lw=2,label=f"r={r:+.2f} p={p:.2g}")
    ax.axhline(0,color='gray',ls=':');ax.set_xlabel("severity");ax.set_ylabel(lab);ax.legend();ax.set_title(lab)
fig.suptitle("I02. Dose-response: diversity change vs disturbance severity (impact plots)",fontsize=13)
fig.tight_layout(); save(fig,"I02_severity_dose.png")

# ---- I03: change<->change ----
fig,axes=plt.subplots(1,4,figsize=(20,4.6))
for ax,(rs,rl) in zip(axes,[('dnbr','ΔNBR'),('dFHD','ΔFHD'),('dDeep_Gap','ΔDeepGap'),('dLAI','ΔLAI')]):
    d=m.dropna(subset=['dq1',rs])
    ax.scatter(d[rs],d['dq1'],s=16,alpha=.4,c=d['impact'],cmap='coolwarm')
    if len(d)>10:
        sl,ic,r,p,se=st.linregress(d[rs],d['dq1']);xx=np.linspace(d[rs].min(),d[rs].max(),20)
        ax.plot(xx,ic+sl*xx,'k-',lw=2,label=f"r={r:+.2f} p={p:.2g}")
    ax.axhline(0,color='gray',ls=':');ax.axvline(0,color='gray',ls=':')
    ax.set_xlabel(rl);ax.set_ylabel("ΔHill q1");ax.legend(fontsize=8);ax.set_title(f"ΔHill q1 ~ {rl}")
fig.suptitle("I03. Change tracks change — ΔTaxonomic diversity ~ ΔRS (red=impact, blue=control)",fontsize=14)
fig.tight_layout(); save(fig,"I03_change_coupling.png")
print("DONE BACI")
