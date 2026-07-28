# Restructured framework: NO productivity block (fragile/version-dependent, coarse scale).
# M0 domain -> M1 state(structure+spectral) -> M2 +dynamics. Variance partition into
# structure/spectral/dynamics/shared. Dynamics block: cluster-robust Wald + wild cluster bootstrap.
import os, warnings, numpy as np, pandas as pd; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"; D=os.path.join(BASE,"data")
OUT=os.path.join(BASE,"papers","UNIFIED","results","O0_framework")
df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled_26.csv"))
df=df[df.sample_coverage>=0.90].copy()
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"]; P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"]
keep=S+P+Dn
RESP={"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"}
LOG={"Hill_q1","Hill_q2"}
def wcb(d,base_terms,block,B=999,seed=11):
    rng=np.random.RandomState(seed)
    f_full="y ~ "+" + ".join(base_terms+block); f_red="y ~ "+" + ".join(base_terms)
    m=smf.ols(f_full,d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
    Fobs=float(m.wald_test(",".join(f"{t}=0" for t in block),use_f=True,scalar=True).statistic)
    r0=smf.ols(f_red,d).fit(); fit0=r0.fittedvalues; res0=r0.resid; cl=d.siteID.values; uniq=np.unique(cl)
    cnt=0
    for _ in range(B):
        w=pd.Series(np.where(rng.rand(len(uniq))<0.5,1,-1),index=uniq)
        d=d.assign(y=fit0+res0*d.siteID.map(w).values)
        mb=smf.ols(f_full,d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
        Fb=float(mb.wald_test(",".join(f"{t}=0" for t in block),use_f=True,scalar=True).statistic)
        if Fb>=Fobs: cnt+=1
    return Fobs,(cnt+1)/(B+1)
nested=[]; vp=[]; cw=[]; wb=[]
for rp,lab in RESP.items():
    d=df[[rp,"domain","siteID"]+keep].dropna().copy()
    for c in keep: d[c]=(d[c]-d[c].mean())/d[c].std()
    y=np.log(d[rp]) if rp in LOG else d[rp].astype(float); d["y"]=(y-y.mean())/y.std()
    o=lambda pr: smf.ols("y ~ C(domain)"+"".join(" + "+c for c in pr),d).fit().rsquared
    r_m0=0.0; r_m1=o(S+P)-0; r_m2=o(S+P+Dn)  # beyond domain
    R2dom=smf.ols("y ~ C(domain)",d).fit().rsquared
    nested.append(dict(response=lab,n=len(d),R2_domain=R2dom,R2_M1_beyond=o(S+P),R2_M2_beyond=o(S+P+Dn)))
    # unique blocks in M2 (beyond domain)
    full=o(S+P+Dn)
    u_str=full-o(P+Dn); u_spe=full-o(S+Dn); u_dyn=full-o(S+P)
    shared=full-(u_str+u_spe+u_dyn)
    vp.append(dict(response=lab,R2_RS=full,unique_structure=u_str,unique_spectral=u_spe,unique_dynamics=u_dyn,shared_RS=shared))
    # dynamics block: cluster-robust Wald + dR2
    m=smf.ols("y ~ C(domain)"+"".join(" + "+c for c in keep),d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
    wd=m.wald_test(",".join(f"{t}=0" for t in Dn),use_f=True,scalar=True)
    cw.append(dict(response=lab,n=len(d),wald_F_dyn=float(wd.statistic),p_dyn_clustered=float(wd.pvalue),dR2_dyn=o(S+P+Dn)-o(S+P)))
    F,p=wcb(d,["C(domain)"]+S+P,Dn)
    wb.append(dict(response=lab,wald_F_dyn=F,p_dyn_wildboot=p))
pd.DataFrame(nested).to_csv(os.path.join(OUT,"nested_models_sd.csv"),index=False)
pd.DataFrame(vp).to_csv(os.path.join(OUT,"variance_partition_sd.csv"),index=False)
pd.DataFrame(cw).to_csv(os.path.join(OUT,"cluster_wald_sd.csv"),index=False)
pd.DataFrame(wb).to_csv(os.path.join(OUT,"wildboot_sd.csv"),index=False)
print("=== RESTRUCTURED FRAMEWORK (state + dynamics, NO productivity) ===")
print(pd.DataFrame(vp).round(4).to_string(index=False))
print()
m1=pd.DataFrame(cw).merge(pd.DataFrame(wb),on="response",suffixes=("_wald","_wb"))
print(m1[["response","n","dR2_dyn","p_dyn_clustered","p_dyn_wildboot"]].round(4).to_string(index=False))
