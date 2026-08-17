"""
Reviewer consistency fixes:
 (A) O4 interactions recomputed with the SAME framework predictors as the nested models
     (spectral = EVI only, retained structure, GPP) x context (stand age, disturbance, land use),
     cluster-robust SE, Benjamini-Hochberg per family. Replaces the SAVI-based O4 reporting.
 (B) Press vs pulse structural change recomputed at SITE level (site-median post-pre change),
     avoiding plot-level pseudoreplication.
Out: papers/UNIFIED/results/O0_framework/o4_interactions_consistent.csv, press_pulse_site.csv
"""
import os, warnings, numpy as np, pandas as pd, scipy.stats as st
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"; D=os.path.join(BASE,"data")
OUT=os.path.join(BASE,"papers","UNIFIED","results","O0_framework")

df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
for f in ["plot_modis_gpp.csv","plot_pml_gpp.csv"]: df=df.merge(pd.read_csv(os.path.join(D,f)),on="plotID",how="left")
df=df.merge(pd.read_csv(os.path.join(D,"plot_stand_age_gami.csv"))[["plotID","stand_age_gami"]],on="plotID",how="left")
df=df.merge(pd.read_csv(os.path.join(D,"plot_disturbance_robust.csv"))[["plotID","severity","recency"]],on="plotID",how="left")
df=df.merge(pd.read_csv(os.path.join(D,"plot_landuse_het.csv"))[["plotID","lc_edge","lc_shannon","lc_forest_frac"]],on="plotID",how="left")
df=df[df.sample_coverage>=0.90].copy()

RS=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean","EVI_mean","modis_gpp","pml_gpp"]   # framework set (EVI only)
RESP={"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"turnover","LCBD_nestedness_rare":"nestedness"}
LOG={"Hill_q1","Hill_q2"}
FAM={"stand age":["stand_age_gami"], "disturbance":["severity","recency"], "land use":["lc_edge","lc_shannon","lc_forest_frac"]}

def z(s): return (s-s.mean())/s.std()
rows=[]
for fam,ctxs in FAM.items():
    for ctx in ctxs:
        for rp,lab in RESP.items():
            for rs in RS:
                d=df[[rp,"domain","siteID",rs,ctx]].dropna().copy()
                if d.siteID.nunique()<5 or len(d)<40: continue
                d["y"]=np.log(d[rp]) if rp in LOG else d[rp].astype(float)
                d["zR"]=z(d[rs]); d["zC"]=z(d[ctx])
                try:
                    m=smf.ols("y ~ C(domain) + zR + zC + zR:zC",d).fit(cov_type="cluster",cov_kwds={"groups":d["siteID"]})
                    rows.append(dict(family=fam,context=ctx,response=lab,rs=rs,
                        beta_int=m.params.get("zR:zC",np.nan),p=m.pvalues.get("zR:zC",np.nan),n=len(d)))
                except Exception: pass
I=pd.DataFrame(rows).dropna(subset=["p"])
# BH within each family
I["q"]=np.nan
for fam in FAM:
    mask=I.family==fam
    I.loc[mask,"q"]=multipletests(I.loc[mask,"p"].values,method="fdr_bh")[1]
I["sig"]=I["q"]<0.05
I.to_csv(os.path.join(OUT,"o4_interactions_consistent.csv"),index=False)
print("=== O4 interactions (framework predictors incl EVI only) ===")
print(f"total tests={len(I)} | survive BH per family: {int(I.sig.sum())}")
for fam in FAM:
    s=I[(I.family==fam)&I.sig]
    print(f"\n[{fam}] {len(s)}/{(I.family==fam).sum()} sig")
    print(s[["context","response","rs","beta_int","p","q"]].round(4).to_string(index=False) if len(s) else "  (none)")

# ---------- (B) press vs pulse at SITE level ----------
neon=pd.read_csv(os.path.join(D,"plot_disturbance_neon.csv"))
imp=neon[neon.disturbed==1][["plotID","dist_year","neon_dist_type"]].dropna(subset=["dist_year"])
PRESS={"insect","mortality","natural"}; PULSE={"fire","wind","flood","harvest"}
imp["regime"]=imp.neon_dist_type.map(lambda t:"press" if t in PRESS else ("pulse" if t in PULSE else "other"))
splity=imp.set_index("plotID")["dist_year"].to_dict()
py=pd.read_csv(os.path.join(D,"per_year_v2.csv")); rows2=[]
for pid,g in py.groupby("plotID"):
    if pid not in splity: continue
    sy=splity[pid]; b=g[g.year<sy]; a=g[g.year>=sy]
    if len(b)<1 or len(a)<1 or "LAI" not in g: continue
    rows2.append(dict(plotID=pid,siteID=pid[:4],dLAI=a.LAI.mean()-b.LAI.mean()))
ch=pd.DataFrame(rows2).merge(imp[["plotID","regime"]],on="plotID",how="left")
site=ch.groupby(["siteID","regime"]).dLAI.median().reset_index()   # site-median change per regime
pr=site[site.regime=="press"].dLAI; pu=site[site.regime=="pulse"].dLAI
res=dict(press_n_sites=len(pr),pulse_n_sites=len(pu),press_med=pr.median(),pulse_med=pu.median())
if len(pr)>=3 and len(pu)>=3:
    res["mw_p"]=st.mannwhitneyu(pr,pu).pvalue
pd.DataFrame([res]).to_csv(os.path.join(OUT,"press_pulse_site.csv"),index=False)
print("\n=== press vs pulse ΔLAI (post-pre) at SITE level ===")
print(res)
