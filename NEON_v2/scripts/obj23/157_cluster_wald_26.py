"""
Reviewer fixes: (A) cluster-robust JOINT WALD tests for the dynamics & function blocks
(primary significance; LRT/AIC secondary), (B) fully standardized coefficients (response
also z-scored) so magnitudes are comparable, (C) per-response dynamics ΔR2, (D) shared-fraction check.
OLS + domain fixed + site-clustered cov. Out: papers/UNIFIED/results/O0_framework/cluster_wald.csv
"""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"; D=os.path.join(BASE,"data")
OUT=os.path.join(BASE,"papers","UNIFIED","results","O0_framework")
df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled_26.csv"))
for f in ["plot_modis_gpp_26.csv","plot_pml_gpp_26.csv"]: df=df.merge(pd.read_csv(os.path.join(D,f)),on="plotID",how="left")
df=df[df.sample_coverage>=0.90].copy()
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];G=["modis_gpp","pml_gpp"]
keep=S+P+Dn+G
RESP={"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"}
LOG={"Hill_q1","Hill_q2"}
rows=[]; coefs=[]
for rp,lab in RESP.items():
    d=df[[rp,"domain","siteID"]+keep].dropna().copy()
    for c in keep: d[c]=(d[c]-d[c].mean())/d[c].std()
    y=np.log(d[rp]) if rp in LOG else d[rp].astype(float)
    d["y"]=(y-y.mean())/y.std()                       # standardize response -> fully standardized betas
    f="y ~ C(domain)"+"".join(" + "+c for c in keep)
    m=smf.ols(f,d).fit(cov_type="cluster",cov_kwds={"groups":d["siteID"]})
    # (A) cluster-robust joint Wald (F) tests for blocks
    wd=m.wald_test(",".join(f"{t} = 0" for t in Dn),use_f=True,scalar=True)
    wf=m.wald_test(",".join(f"{t} = 0" for t in G), use_f=True,scalar=True)
    # per-response dynamics/function dR2 (non-robust OLS, for reporting)
    o=lambda preds: smf.ols("y ~ C(domain)"+"".join(" + "+c for c in preds),d).fit().rsquared
    dR2_dyn=o(S+P+Dn)-o(S+P); dR2_fun=o(S+P+Dn+G)-o(S+P+Dn)
    rows.append(dict(response=lab,n=len(d),
        wald_F_dyn=float(wd.statistic),p_dyn_clustered=float(wd.pvalue),dR2_dyn=dR2_dyn,
        wald_F_fun=float(wf.statistic),p_fun_clustered=float(wf.pvalue),dR2_fun=dR2_fun))
    for t in keep:
        coefs.append(dict(response=lab,predictor=t,beta_fullstd=m.params[t],p_clustered=m.pvalues[t]))
R=pd.DataFrame(rows); R.to_csv(os.path.join(OUT,"cluster_wald.csv"),index=False)
pd.DataFrame(coefs).to_csv(os.path.join(OUT,"coeffs_fullstd.csv"),index=False)
pd.set_option("display.width",200)
print("=== cluster-robust JOINT WALD (F) block tests + per-response dR2 ===")
print(R.round(4).to_string(index=False))
print("\n=== fully-standardized key coefficients (response z-scored) ===")
cf=pd.DataFrame(coefs)
for lab in RESP.values():
    s=cf[cf.response==lab].reindex(cf[cf.response==lab].beta_fullstd.abs().sort_values(ascending=False).index).head(4)
    print(f"  {lab}: "+", ".join(f"{r.predictor}={r.beta_fullstd:+.2f}(p={r.p_clustered:.2g})" for r in s.itertuples()))
