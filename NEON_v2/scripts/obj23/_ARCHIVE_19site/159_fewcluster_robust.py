"""
Reviewer statistical hardening (19 clusters):
 (A) WILD CLUSTER BOOTSTRAP (Rademacher, restricted) for the dynamics & function block joint tests
     -> few-cluster-valid p-values vs the asymptotic cluster-robust Wald.
 (B) turnover GPP suppression diagnostics: MODIS-only, PML-only, joint coefs, r, VIF.
 (C) simple slopes for the two confirmatory interactions (GPP x age, structure x disturbance).
Out: papers/UNIFIED/results/O0_framework/wildboot_blocks.csv, gpp_suppression.csv, simple_slopes.csv
"""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
np.random.seed(7)
BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"; D=os.path.join(BASE,"data")
OUT=os.path.join(BASE,"papers","UNIFIED","results","O0_framework")
df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
for f in ["plot_modis_gpp.csv","plot_pml_gpp.csv"]: df=df.merge(pd.read_csv(os.path.join(D,f)),on="plotID",how="left")
df=df.merge(pd.read_csv(os.path.join(D,"plot_stand_age_gami.csv"))[["plotID","stand_age_gami"]],on="plotID",how="left")
df=df.merge(pd.read_csv(os.path.join(D,"plot_disturbance_robust.csv"))[["plotID","severity"]],on="plotID",how="left")
df=df[df.sample_coverage>=0.90].copy()
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];G=["modis_gpp","pml_gpp"]
keep=S+P+Dn+G
RESP={"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"}
LOG={"Hill_q1","Hill_q2"}
def z(s): return (s-s.mean())/s.std()

# ---------- (A) wild cluster bootstrap ----------
def wald_F(res,terms): return float(res.wald_test(",".join(f"{t}=0" for t in terms),use_f=True,scalar=True).statistic)
def wcb(d,base,block,B=999):
    g=d["siteID"]; cl=g.unique()
    ff="y ~ "+" + ".join(base+block); rf="y ~ "+" + ".join(base)
    full=smf.ols(ff,d).fit(cov_type="cluster",cov_kwds={"groups":g}); Fobs=wald_F(full,block)
    r0=smf.ols(rf,d).fit(); fit0=r0.fittedvalues.values; res0=r0.resid.values
    Fs=np.empty(B)
    for b in range(B):
        w=dict(zip(cl,np.where(np.random.rand(len(cl))<0.5,1.0,-1.0)))
        d2=d.copy(); d2["y"]=fit0+g.map(w).values*res0
        Fs[b]=wald_F(smf.ols(ff,d2).fit(cov_type="cluster",cov_kwds={"groups":g}),block)
    return Fobs,(np.sum(Fs>=Fobs)+1)/(B+1)
rowsA=[]
for rp,lab in RESP.items():
    d=df[[rp,"domain","siteID"]+keep].dropna().copy()
    for c in keep: d[c]=z(d[c])
    d["y"]=np.log(d[rp]) if rp in LOG else d[rp].astype(float); d["y"]=z(d["y"])
    F1,p1=wcb(d,["C(domain)"]+S+P,Dn)               # dynamics beyond state
    F2,p2=wcb(d,["C(domain)"]+S+P+Dn,G)             # function beyond state+dynamics
    rowsA.append(dict(response=lab,wald_F_dyn=F1,p_dyn_wildboot=p1,wald_F_fun=F2,p_fun_wildboot=p2))
    print(f"[WCB] {lab:16s} dyn F={F1:.2f} p={p1:.3f} | fun F={F2:.2f} p={p2:.3f}",flush=True)
pd.DataFrame(rowsA).to_csv(os.path.join(OUT,"wildboot_blocks.csv"),index=False)

# ---------- (B) turnover GPP suppression ----------
d=df[["LCBD_turnover_rare","domain","siteID"]+keep].dropna().copy()
for c in keep: d[c]=z(d[c]);
d["y"]=d["LCBD_turnover_rare"]
base="y ~ C(domain) + "+" + ".join(S+P+Dn)
mM=smf.ols(base+" + modis_gpp",d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
mP=smf.ols(base+" + pml_gpp",d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
mB=smf.ols(base+" + modis_gpp + pml_gpp",d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
vif_g=variance_inflation_factor(d[["modis_gpp","pml_gpp"]].assign(c=1).values,0)
sup=dict(r_modis_pml=d.modis_gpp.corr(d.pml_gpp),
    modis_only_beta=mM.params["modis_gpp"],modis_only_p=mM.pvalues["modis_gpp"],
    pml_only_beta=mP.params["pml_gpp"],pml_only_p=mP.pvalues["pml_gpp"],
    modis_joint=mB.params["modis_gpp"],pml_joint=mB.params["pml_gpp"],vif_between_gpp=vif_g)
pd.DataFrame([sup]).to_csv(os.path.join(OUT,"gpp_suppression.csv"),index=False)
print("\n[SUPPRESSION turnover]",{k:round(v,3) for k,v in sup.items()})

# ---------- (C) simple slopes ----------
def simple_slopes(rp,rs,ctx,vals=(-1,0,1)):
    d=df[[rp,"domain","siteID",rs,ctx]].dropna().copy()
    d["zR"]=z(d[rs]); d["zC"]=z(d[ctx]); d["y"]=np.log(d[rp]) if rp in LOG else d[rp].astype(float); d["y"]=z(d["y"])
    m=smf.ols("y ~ C(domain) + zR + zC + zR:zC",d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
    b1,b3=m.params["zR"],m.params["zR:zC"]; V=m.cov_params()
    out=[]
    for v in vals:
        sl=b1+b3*v; se=np.sqrt(V.loc["zR","zR"]+v*v*V.loc["zR:zC","zR:zC"]+2*v*V.loc["zR","zR:zC"])
        out.append(dict(rp=RESP[rp],rs=rs,ctx=ctx,ctx_z=v,slope=sl,se=se,p=2*(1-__import__("scipy.stats",fromlist=["norm"]).norm.cdf(abs(sl/se)))))
    return out
ss=[]
ss+=simple_slopes("Hill_q1","modis_gpp","stand_age_gami")
ss+=simple_slopes("Hill_q1","VCI_mean","severity")
S_=pd.DataFrame(ss); S_.to_csv(os.path.join(OUT,"simple_slopes.csv"),index=False)
print("\n[SIMPLE SLOPES]"); print(S_.round(3).to_string(index=False))
